#include "EquilibriumThermoModel.H"
#include "SpeciesRegistry.H"
#include "TransportProps.H"

#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <new>
#include <type_traits>
#include <vector>

namespace
{

std::atomic<bool> countAllocations{false};
std::atomic<std::size_t> allocationCount{0};

void recordAllocation() noexcept
{
    if(countAllocations.load(std::memory_order_relaxed))
        allocationCount.fetch_add(1,std::memory_order_relaxed);
}

void* allocate(std::size_t size)
{
    void* memory=std::malloc(size==0 ? 1 : size);
    if(!memory) throw std::bad_alloc();
    recordAllocation();
    return memory;
}

} // namespace

void* operator new(std::size_t size)
{
    return allocate(size);
}

void* operator new[](std::size_t size)
{
    return allocate(size);
}

void operator delete(void* memory) noexcept
{
    std::free(memory);
}

void operator delete[](void* memory) noexcept
{
    std::free(memory);
}

void operator delete(void* memory,std::size_t) noexcept
{
    std::free(memory);
}

void operator delete[](void* memory,std::size_t) noexcept
{
    std::free(memory);
}

using namespace Foam;

int main()
{
    static_assert(
        std::is_trivially_copyable<TransportProps::Result>::value,
        "direct-Q transport results must not own dynamic storage");

    SpeciesRegistry registry;
    const std::vector<std::string> names={"O2","N2"};
    const std::vector<double> massFractions={0.7,0.3};
    const auto model=thermo::makeLegacySinglePhaseThermoModel(names);
    const auto initialized=model->initializePTY(
        1.0e6,500.0,massFractions);
    if(!initialized.success())
    {
        std::cerr<<"legacy transport allocation fixture failed to initialize\n";
        return 1;
    }

    PrimitiveVars primitive(static_cast<int>(names.size()));
    primitive.p=initialized.state->p;
    primitive.T=initialized.state->T;
    primitive.Y=massFractions;

    const TransportProps::Result warm=
        TransportProps::evaluateCommittedLegacyGas(
            names,primitive,*initialized.state,*model,registry);
    if(!warm.success())
    {
        std::cerr<<"legacy transport allocation fixture failed to warm up\n";
        return 1;
    }

    allocationCount.store(0,std::memory_order_relaxed);
    countAllocations.store(true,std::memory_order_relaxed);
    bool success=true;
    double checksum=0.0;
    for(int repetition=0;repetition<512;++repetition)
    {
        const TransportProps::Result evaluated=
            TransportProps::evaluateCommittedLegacyGas(
                names,primitive,*initialized.state,*model,registry);
        success=success&&evaluated.success();
        checksum+=evaluated.point.mu;
    }
    countAllocations.store(false,std::memory_order_relaxed);

    const std::size_t observed=
        allocationCount.load(std::memory_order_relaxed);
    if(!success || !(checksum>0.0) || observed!=0)
    {
        std::cerr<<"direct-Q transport allocated "<<observed
                 <<" times after warm-up\n";
        return 1;
    }

    std::cout<<"TRANSPORT_NO_ALLOC_PASS\n";
    return 0;
}
