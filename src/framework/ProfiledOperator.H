#ifndef ProfiledOperator_H
#define ProfiledOperator_H

#include "Operator.H"
#include <chrono>
#include <memory>

namespace Foam
{

class ProfiledOperator : public Operator
{
    std::unique_ptr<Operator> inner_;
    mutable double applySeconds_=0.0,timestepSeconds_=0.0,lastStableDt_=1e30;
    mutable long applyCalls_=0,timestepCalls_=0;
public:
    explicit ProfiledOperator(std::unique_ptr<Operator> inner):inner_(std::move(inner)) {}
    Operator* inner() const { return inner_.get(); }
    bool ownsPersistentState() const override
    { return inner_->ownsPersistentState(); }
    void beginStepTransaction(const Domain& d) const override
    { inner_->beginStepTransaction(d); }
    void acceptStepTransaction() const noexcept override
    { inner_->acceptStepTransaction(); }
    void rollbackStepTransaction() const noexcept override
    { inner_->rollbackStepTransaction(); }
    void averageWithStepStartState() const override
    { inner_->averageWithStepStartState(); }
    void apply(Domain& d,double dt) const override
    {
        const auto start=std::chrono::steady_clock::now(); inner_->apply(d,dt);
        applySeconds_+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        ++applyCalls_;
    }
    double maxStableDt(const Domain& d,double cfl) const override
    {
        const auto start=std::chrono::steady_clock::now();
        const double result=inner_->maxStableDt(d,cfl);
        timestepSeconds_+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        ++timestepCalls_; lastStableDt_=result; return result;
    }
    void localStableDt(
        const Domain& d,double cfl,std::vector<double>& dt) const override
    {
        const auto start=std::chrono::steady_clock::now();
        inner_->localStableDt(d,cfl,dt);
        timestepSeconds_+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-start).count();
        ++timestepCalls_;
        lastStableDt_=dt.empty()?1e30:
            *std::min_element(dt.begin(),dt.end());
    }
    void addStructuredImplicitStencil(
        const Domain& d,StructuredImplicitStencil& stencil) const override
    {
        const auto start=std::chrono::steady_clock::now();
        inner_->addStructuredImplicitStencil(d,stencil);
        timestepSeconds_+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-start).count();
        ++timestepCalls_;
    }
    bool supportsStructuredBlockCoupling() const override
    { return inner_->supportsStructuredBlockCoupling(); }
    bool supportsStructuredDiagonalBlock() const override
    { return inner_->supportsStructuredDiagonalBlock(); }
    bool addStructuredImplicitDiagonalBlock(
        const Domain& d,int cell,
        StructuredImplicitBlock& contribution,
        double& scalarContribution) const override
    {
        return inner_->addStructuredImplicitDiagonalBlock(
            d,cell,contribution,scalarContribution);
    }
    bool addStructuredImplicitNeighborContribution(
        const Domain& d,int cell,int neighbor,int face,
        const ConservedVars& neighborCorrection,
        ConservedVars& contribution) const override
    {
        return inner_->addStructuredImplicitNeighborContribution(
            d,cell,neighbor,face,neighborCorrection,contribution);
    }
    std::string name() const override { return inner_->name(); }
    double applySeconds() const { return applySeconds_; }
    double timestepSeconds() const { return timestepSeconds_; }
    double lastStableDt() const { return lastStableDt_; }
    long applyCalls() const { return applyCalls_; }
    long timestepCalls() const { return timestepCalls_; }
};

} // End namespace Foam
#endif
