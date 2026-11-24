#include "RestartManifestCodec.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <new>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace Foam
{
namespace
{

constexpr std::array<unsigned char,8> diskMagic = {{
    'R','F','A','M','R','M','2','\0'
}};

constexpr std::uint64_t maximumDocumentBytes = 32ULL*1024ULL*1024ULL*1024ULL;
constexpr std::uint64_t maximumStringBytes = 4ULL*1024ULL*1024ULL;
constexpr std::uint64_t maximumPayloadBytes = 16ULL*1024ULL*1024ULL*1024ULL;
constexpr std::uint32_t maximumCollectionEntries = 16U*1024U*1024U;
constexpr std::size_t diskHeaderBytes = 8U+4U+8U;

enum class WireType : std::uint8_t
{
    U8 = 1,
    U32 = 2,
    U64 = 3,
    F64 = 4,
    BOOL = 5,
    STRING = 6,
    BYTES = 7,
    RECORD = 8,
    RECORD_LIST = 9
};

bool validWireType(std::uint8_t raw)
{
    return raw>=static_cast<std::uint8_t>(WireType::U8)
        &&raw<=static_cast<std::uint8_t>(WireType::RECORD_LIST);
}

struct DecodeError
{
    RestartCodecStatus status = RestartCodecStatus::SUCCESS;
    std::string detail;

    bool fail(RestartCodecStatus next,const std::string& message)
    {
        if(status==RestartCodecStatus::SUCCESS)
        {
            status=next;
            detail=message;
        }
        return false;
    }
};

class Reader
{
    std::string_view bytes_;
    std::size_t position_ = 0;

public:
    explicit Reader(std::string_view bytes):bytes_(bytes) {}

    std::size_t remaining() const
    {
        return bytes_.size()-position_;
    }

    bool readU8(std::uint8_t& value)
    {
        if(remaining()<1U) return false;
        value=static_cast<std::uint8_t>(
            static_cast<unsigned char>(bytes_[position_++]));
        return true;
    }

    bool readU16(std::uint16_t& value)
    {
        if(remaining()<2U) return false;
        value=0U;
        for(unsigned int shift=0U;shift<16U;shift+=8U)
        {
            value|=static_cast<std::uint16_t>(
                static_cast<unsigned char>(bytes_[position_++]))<<shift;
        }
        return true;
    }

    bool readU32(std::uint32_t& value)
    {
        if(remaining()<4U) return false;
        value=0U;
        for(unsigned int shift=0U;shift<32U;shift+=8U)
        {
            value|=static_cast<std::uint32_t>(
                static_cast<unsigned char>(bytes_[position_++]))<<shift;
        }
        return true;
    }

    bool readU64(std::uint64_t& value)
    {
        if(remaining()<8U) return false;
        value=0U;
        for(unsigned int shift=0U;shift<64U;shift+=8U)
        {
            value|=static_cast<std::uint64_t>(
                static_cast<unsigned char>(bytes_[position_++]))<<shift;
        }
        return true;
    }

    bool readView(std::uint64_t length,std::string_view& value)
    {
        if(length>static_cast<std::uint64_t>(remaining())
           ||length>static_cast<std::uint64_t>(
                std::numeric_limits<std::size_t>::max()))
            return false;
        const std::size_t count=static_cast<std::size_t>(length);
        value=bytes_.substr(position_,count);
        position_+=count;
        return true;
    }
};

void requireAppendCapacity(const std::string& output,std::uint64_t count)
{
    if(count>maximumDocumentBytes
       ||static_cast<std::uint64_t>(output.size())
            >maximumDocumentBytes-count
       ||count>static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max()-output.size()))
        throw std::length_error("restart codec output exceeds resource limit");
}

void appendRaw(std::string& output,std::string_view bytes)
{
    requireAppendCapacity(output,bytes.size());
    output.append(bytes.data(),bytes.size());
}

void appendU8(std::string& output,std::uint8_t value)
{
    requireAppendCapacity(output,1U);
    output.push_back(static_cast<char>(value));
}

void appendU16(std::string& output,std::uint16_t value)
{
    requireAppendCapacity(output,2U);
    for(unsigned int shift=0U;shift<16U;shift+=8U)
        output.push_back(static_cast<char>((value>>shift)&0xffU));
}

void appendU32(std::string& output,std::uint32_t value)
{
    requireAppendCapacity(output,4U);
    for(unsigned int shift=0U;shift<32U;shift+=8U)
        output.push_back(static_cast<char>((value>>shift)&0xffU));
}

void appendU64(std::string& output,std::uint64_t value)
{
    requireAppendCapacity(output,8U);
    for(unsigned int shift=0U;shift<64U;shift+=8U)
        output.push_back(static_cast<char>((value>>shift)&0xffU));
}

std::string scalarU8(std::uint8_t value)
{
    std::string bytes;
    appendU8(bytes,value);
    return bytes;
}

std::string scalarU32(std::uint32_t value)
{
    std::string bytes;
    appendU32(bytes,value);
    return bytes;
}

std::string scalarU64(std::uint64_t value)
{
    std::string bytes;
    appendU64(bytes,value);
    return bytes;
}

std::string scalarF64(double value)
{
    if(!std::numeric_limits<double>::is_iec559
       ||sizeof(double)!=sizeof(std::uint64_t))
        throw std::logic_error(
            "restart codec host does not support IEEE binary64");
    std::uint64_t bits=0U;
    std::memcpy(&bits,&value,sizeof(bits));
    return scalarU64(bits);
}

std::string scalarBool(bool value)
{
    return scalarU8(value?1U:0U);
}

struct EncodedField
{
    std::uint16_t id = 0U;
    WireType type = WireType::BYTES;
    std::string data;
};

std::string encodedRecord(std::vector<EncodedField> fields)
{
    if(fields.size()>std::numeric_limits<std::uint16_t>::max())
        throw std::length_error("restart codec record has too many fields");
    std::string bytes;
    appendU16(bytes,static_cast<std::uint16_t>(fields.size()));
    std::uint16_t previous=0U;
    for(std::size_t i=0;i<fields.size();++i)
    {
        const EncodedField& field=fields[i];
        if(field.id==0U||(i>0U&&field.id<=previous))
            throw std::logic_error("restart codec field order is not canonical");
        if(field.type==WireType::STRING
           &&field.data.size()>maximumStringBytes)
            throw std::length_error(
                "restart codec string exceeds resource limit");
        if(field.type==WireType::BYTES
           &&field.data.size()>maximumPayloadBytes)
            throw std::length_error(
                "restart codec payload exceeds resource limit");
        previous=field.id;
        appendU16(bytes,field.id);
        appendU8(bytes,static_cast<std::uint8_t>(field.type));
        appendU8(bytes,0U);
        appendU64(bytes,field.data.size());
        appendRaw(bytes,field.data);
    }
    return bytes;
}

std::string encodedRecordList(const std::vector<std::string>& records)
{
    if(records.size()>std::numeric_limits<std::uint32_t>::max()
       ||records.size()>maximumCollectionEntries)
        throw std::length_error("restart codec list exceeds resource limit");
    std::string bytes;
    appendU32(bytes,static_cast<std::uint32_t>(records.size()));
    for(const std::string& record:records)
    {
        appendU64(bytes,record.size());
        appendRaw(bytes,record);
    }
    return bytes;
}

struct FieldView
{
    std::uint16_t id = 0U;
    WireType type = WireType::BYTES;
    std::string_view data;
};

struct RecordView
{
    std::vector<FieldView> fields;
};

bool parseRecord(
    std::string_view bytes,
    RecordView& record,
    DecodeError& error,
    const std::string& context)
{
    Reader reader(bytes);
    std::uint16_t count=0U;
    if(!reader.readU16(count))
        return error.fail(RestartCodecStatus::TRUNCATED_INPUT,
            context+": missing record field count");
    if(count>256U)
        return error.fail(RestartCodecStatus::RESOURCE_LIMIT,
            context+": record field count exceeds codec limit");
    record.fields.clear();
    record.fields.reserve(count);
    std::vector<std::uint16_t> seen;
    seen.reserve(count);
    std::uint16_t previous=0U;
    for(std::uint16_t i=0U;i<count;++i)
    {
        std::uint16_t id=0U;
        std::uint8_t rawType=0U;
        std::uint8_t reserved=0U;
        std::uint64_t length=0U;
        if(!reader.readU16(id)||!reader.readU8(rawType)
           ||!reader.readU8(reserved)||!reader.readU64(length))
            return error.fail(RestartCodecStatus::TRUNCATED_INPUT,
                context+": truncated field header");
        if(id==0U)
            return error.fail(RestartCodecStatus::UNKNOWN_FIELD,
                context+": field id zero is reserved");
        if(std::find(seen.begin(),seen.end(),id)!=seen.end())
            return error.fail(RestartCodecStatus::DUPLICATE_FIELD,
                context+": duplicate field id "+std::to_string(id));
        if(i>0U&&id<=previous)
            return error.fail(RestartCodecStatus::NONCANONICAL_FIELD_ORDER,
                context+": fields are not strictly ordered");
        if(!validWireType(rawType))
            return error.fail(RestartCodecStatus::TYPE_MISMATCH,
                context+": unknown wire type");
        if(reserved!=0U)
            return error.fail(RestartCodecStatus::MALFORMED_RECORD,
                context+": reserved field byte is nonzero");
        std::string_view data;
        if(!reader.readView(length,data))
            return error.fail(RestartCodecStatus::TRUNCATED_INPUT,
                context+": field payload exceeds record");
        seen.push_back(id);
        previous=id;
        record.fields.push_back(
            {id,static_cast<WireType>(rawType),data});
    }
    if(reader.remaining()!=0U)
        return error.fail(RestartCodecStatus::TRAILING_DATA,
            context+": record has an undecoded suffix");
    return true;
}

struct ExpectedField
{
    std::uint16_t id;
    WireType type;
};

bool exactFields(
    const RecordView& record,
    std::initializer_list<ExpectedField> expected,
    DecodeError& error,
    const std::string& context)
{
    for(const FieldView& field:record.fields)
    {
        const auto found=std::find_if(
            expected.begin(),expected.end(),[&](const ExpectedField& item)
            { return item.id==field.id; });
        if(found==expected.end())
            return error.fail(RestartCodecStatus::UNKNOWN_FIELD,
                context+": unknown field id "+std::to_string(field.id));
        if(found->type!=field.type)
            return error.fail(RestartCodecStatus::TYPE_MISMATCH,
                context+": field "+std::to_string(field.id)
                    +" has the wrong wire type");
    }
    for(const ExpectedField& item:expected)
    {
        const bool found=std::any_of(
            record.fields.begin(),record.fields.end(),[&](const FieldView& field)
            { return field.id==item.id; });
        if(!found)
            return error.fail(RestartCodecStatus::MISSING_FIELD,
                context+": missing field id "+std::to_string(item.id));
    }
    return true;
}

const FieldView& field(const RecordView& record,std::uint16_t id)
{
    const auto found=std::find_if(
        record.fields.begin(),record.fields.end(),[&](const FieldView& item)
        { return item.id==id; });
    if(found==record.fields.end())
        throw std::logic_error("restart codec field lookup failed");
    return *found;
}

bool decodeU8(
    const FieldView& input,
    std::uint8_t& value,
    DecodeError& error,
    const std::string& context)
{
    if(input.data.size()!=1U)
        return error.fail(RestartCodecStatus::TYPE_MISMATCH,
            context+": u8 field length is not one");
    value=static_cast<std::uint8_t>(
        static_cast<unsigned char>(input.data[0]));
    return true;
}

bool decodeU32(
    const FieldView& input,
    std::uint32_t& value,
    DecodeError& error,
    const std::string& context)
{
    Reader reader(input.data);
    if(input.data.size()!=4U||!reader.readU32(value))
        return error.fail(RestartCodecStatus::TYPE_MISMATCH,
            context+": u32 field length is not four");
    return true;
}

bool decodeU64(
    const FieldView& input,
    std::uint64_t& value,
    DecodeError& error,
    const std::string& context)
{
    Reader reader(input.data);
    if(input.data.size()!=8U||!reader.readU64(value))
        return error.fail(RestartCodecStatus::TYPE_MISMATCH,
            context+": u64 field length is not eight");
    return true;
}

bool decodeF64(
    const FieldView& input,
    double& value,
    DecodeError& error,
    const std::string& context)
{
    if(!std::numeric_limits<double>::is_iec559||sizeof(double)!=8U)
        return error.fail(RestartCodecStatus::RANGE_ERROR,
            context+": host does not support IEEE binary64");
    std::uint64_t bits=0U;
    if(!decodeU64(input,bits,error,context)) return false;
    std::memcpy(&value,&bits,sizeof(value));
    return true;
}

bool decodeBool(
    const FieldView& input,
    bool& value,
    DecodeError& error,
    const std::string& context)
{
    std::uint8_t raw=0U;
    if(!decodeU8(input,raw,error,context)) return false;
    if(raw>1U)
        return error.fail(RestartCodecStatus::RANGE_ERROR,
            context+": boolean value is outside {0,1}");
    value=raw==1U;
    return true;
}

bool decodeString(
    const FieldView& input,
    std::string& value,
    DecodeError& error,
    const std::string& context)
{
    if(input.data.size()>maximumStringBytes)
        return error.fail(RestartCodecStatus::RESOURCE_LIMIT,
            context+": string exceeds codec limit");
    value.assign(input.data.data(),input.data.size());
    return true;
}

bool decodeBytes(
    const FieldView& input,
    std::string& value,
    DecodeError& error,
    const std::string& context)
{
    if(input.data.size()>maximumPayloadBytes)
        return error.fail(RestartCodecStatus::RESOURCE_LIMIT,
            context+": payload exceeds codec limit");
    value.assign(input.data.data(),input.data.size());
    return true;
}

bool parseRecordList(
    const FieldView& input,
    std::vector<std::string_view>& records,
    DecodeError& error,
    const std::string& context)
{
    Reader reader(input.data);
    std::uint32_t count=0U;
    if(!reader.readU32(count))
        return error.fail(RestartCodecStatus::TRUNCATED_INPUT,
            context+": missing list count");
    if(count>maximumCollectionEntries)
        return error.fail(RestartCodecStatus::RESOURCE_LIMIT,
            context+": list count exceeds codec limit");
    if(static_cast<std::uint64_t>(count)
        >static_cast<std::uint64_t>(reader.remaining())/10U)
        return error.fail(RestartCodecStatus::TRUNCATED_INPUT,
            context+": list count cannot fit in encoded bytes");
    records.clear();
    records.reserve(count);
    for(std::uint32_t i=0U;i<count;++i)
    {
        std::uint64_t length=0U;
        std::string_view record;
        if(!reader.readU64(length)||!reader.readView(length,record))
            return error.fail(RestartCodecStatus::TRUNCATED_INPUT,
                context+": truncated list element");
        records.push_back(record);
    }
    if(reader.remaining()!=0U)
        return error.fail(RestartCodecStatus::TRAILING_DATA,
            context+": list has an undecoded suffix");
    return true;
}

std::string encodeGenerationState(const RestartGenerationState& state)
{
    return encodedRecord({
        {1U,WireType::F64,scalarF64(state.physicalTime)},
        {2U,WireType::U64,scalarU64(state.acceptedIteration)},
        {3U,WireType::BOOL,scalarBool(state.previousAcceptedDtPresent)},
        {4U,WireType::F64,scalarF64(state.previousAcceptedDt)}
    });
}

std::string encodeIdentities(const RestartIdentitySet& identities)
{
    return encodedRecord({
        {1U,WireType::STRING,identities.configurationSha256},
        {2U,WireType::STRING,identities.buildSha256},
        {3U,WireType::STRING,identities.meshSha256},
        {4U,WireType::STRING,identities.thermoModelId},
        {5U,WireType::STRING,identities.kijManifestSha256},
        {6U,WireType::STRING,identities.mechanismSha256},
        {7U,WireType::STRING,identities.chemistryIntegratorSha256},
        {8U,WireType::STRING,identities.transportSha256},
        {9U,WireType::STRING,identities.topologyEnvelopeSha256}
    });
}

std::string encodeOwner(const RestartOwnerRecord& owner)
{
    return encodedRecord({
        {1U,WireType::U8,scalarU8(static_cast<std::uint8_t>(owner.owner))},
        {2U,WireType::U8,scalarU8(
            static_cast<std::uint8_t>(owner.disposition))},
        {3U,WireType::U32,scalarU32(owner.publicationOrder)},
        {4U,WireType::STRING,owner.generationId},
        {5U,WireType::STRING,owner.identitySha256},
        {6U,WireType::U64,scalarU64(owner.payloadBytes)},
        {7U,WireType::STRING,owner.payloadSha256}
    });
}

std::string encodeSeed(const RestartRecoverySeed& seed)
{
    return encodedRecord({
        {1U,WireType::U64,scalarU64(seed.cell)},
        {2U,WireType::F64,scalarF64(seed.pressure)},
        {3U,WireType::F64,scalarF64(seed.temperature)},
        {4U,WireType::U8,scalarU8(
            static_cast<std::uint8_t>(seed.expectedPhase))},
        {5U,WireType::STRING,seed.phaseIdentitySha256}
    });
}

std::string encodePayload(const RestartPayloadCandidate& payload)
{
    return encodedRecord({
        {1U,WireType::U8,scalarU8(static_cast<std::uint8_t>(payload.owner))},
        {2U,WireType::STRING,payload.generationId},
        {3U,WireType::BYTES,payload.bytes}
    });
}

bool decodeGenerationState(
    const FieldView& input,
    RestartGenerationState& state,
    DecodeError& error)
{
    RecordView record;
    if(!parseRecord(input.data,record,error,"generation-state")
       ||!exactFields(record,{
            {1U,WireType::F64},
            {2U,WireType::U64},
            {3U,WireType::BOOL},
            {4U,WireType::F64}},error,"generation-state"))
        return false;
    if(!decodeF64(field(record,1U),state.physicalTime,error,
            "generation-state.physical-time")
       ||!decodeU64(field(record,2U),state.acceptedIteration,error,
            "generation-state.accepted-iteration")
       ||!decodeBool(field(record,3U),state.previousAcceptedDtPresent,error,
            "generation-state.previous-dt-present")
       ||!decodeF64(field(record,4U),state.previousAcceptedDt,error,
            "generation-state.previous-dt"))
        return false;
    if(!std::isfinite(state.physicalTime)||state.physicalTime<0.0
       ||(state.physicalTime==0.0&&std::signbit(state.physicalTime))
       ||state.acceptedIteration>static_cast<std::uint64_t>(
            std::numeric_limits<int>::max())
       ||(state.previousAcceptedDtPresent
            &&(!std::isfinite(state.previousAcceptedDt)
               ||state.previousAcceptedDt<=0.0))
       ||(!state.previousAcceptedDtPresent
            &&(state.previousAcceptedDt!=0.0
               ||std::signbit(state.previousAcceptedDt))))
        return error.fail(RestartCodecStatus::RANGE_ERROR,
            "generation-state value is outside the canonical range");
    return true;
}

bool decodeIdentities(
    const FieldView& input,
    RestartIdentitySet& identities,
    DecodeError& error)
{
    RecordView record;
    if(!parseRecord(input.data,record,error,"identities")
       ||!exactFields(record,{
            {1U,WireType::STRING},{2U,WireType::STRING},
            {3U,WireType::STRING},{4U,WireType::STRING},
            {5U,WireType::STRING},{6U,WireType::STRING},
            {7U,WireType::STRING},{8U,WireType::STRING},
            {9U,WireType::STRING}},error,"identities"))
        return false;
    return decodeString(field(record,1U),identities.configurationSha256,error,
            "identities.configuration")
        &&decodeString(field(record,2U),identities.buildSha256,error,
            "identities.build")
        &&decodeString(field(record,3U),identities.meshSha256,error,
            "identities.mesh")
        &&decodeString(field(record,4U),identities.thermoModelId,error,
            "identities.thermo-model")
        &&decodeString(field(record,5U),identities.kijManifestSha256,error,
            "identities.kij")
        &&decodeString(field(record,6U),identities.mechanismSha256,error,
            "identities.mechanism")
        &&decodeString(field(record,7U),
            identities.chemistryIntegratorSha256,error,
            "identities.chemistry-integrator")
        &&decodeString(field(record,8U),identities.transportSha256,error,
            "identities.transport")
        &&decodeString(field(record,9U),identities.topologyEnvelopeSha256,error,
            "identities.topology-envelope");
}

bool decodeOwner(
    std::string_view bytes,
    RestartOwnerRecord& owner,
    DecodeError& error,
    std::size_t index)
{
    const std::string context="owner["+std::to_string(index)+"]";
    RecordView record;
    if(!parseRecord(bytes,record,error,context)
       ||!exactFields(record,{
            {1U,WireType::U8},{2U,WireType::U8},{3U,WireType::U32},
            {4U,WireType::STRING},{5U,WireType::STRING},
            {6U,WireType::U64},{7U,WireType::STRING}},error,context))
        return false;
    std::uint8_t rawOwner=0U;
    std::uint8_t rawDisposition=0U;
    if(!decodeU8(field(record,1U),rawOwner,error,context+".kind")
       ||!decodeU8(field(record,2U),rawDisposition,error,
            context+".disposition")
       ||!decodeU32(field(record,3U),owner.publicationOrder,error,
            context+".publication-order")
       ||!decodeString(field(record,4U),owner.generationId,error,
            context+".generation")
       ||!decodeString(field(record,5U),owner.identitySha256,error,
            context+".identity")
       ||!decodeU64(field(record,6U),owner.payloadBytes,error,
            context+".payload-bytes")
       ||!decodeString(field(record,7U),owner.payloadSha256,error,
            context+".payload-sha256"))
        return false;
    if(rawOwner>=static_cast<std::uint8_t>(RestartOwner::COUNT)
       ||rawDisposition>static_cast<std::uint8_t>(
            RestartOwnerDisposition::STATE_PAYLOAD))
        return error.fail(RestartCodecStatus::RANGE_ERROR,
            context+": enum value is outside the declared domain");
    owner.owner=static_cast<RestartOwner>(rawOwner);
    owner.disposition=static_cast<RestartOwnerDisposition>(rawDisposition);
    return true;
}

bool decodeSeed(
    std::string_view bytes,
    RestartRecoverySeed& seed,
    DecodeError& error,
    std::size_t index)
{
    const std::string context="seed["+std::to_string(index)+"]";
    RecordView record;
    if(!parseRecord(bytes,record,error,context)
       ||!exactFields(record,{
            {1U,WireType::U64},{2U,WireType::F64},{3U,WireType::F64},
            {4U,WireType::U8},{5U,WireType::STRING}},error,context))
        return false;
    std::uint8_t rawPhase=0U;
    if(!decodeU64(field(record,1U),seed.cell,error,context+".cell")
       ||!decodeF64(field(record,2U),seed.pressure,error,context+".pressure")
       ||!decodeF64(field(record,3U),seed.temperature,error,
            context+".temperature")
       ||!decodeU8(field(record,4U),rawPhase,error,context+".phase")
       ||!decodeString(field(record,5U),seed.phaseIdentitySha256,error,
            context+".phase-identity"))
        return false;
    if(rawPhase>static_cast<std::uint8_t>(
            RestartPhaseAssertion::THREE_PHASE_REQUIRED))
        return error.fail(RestartCodecStatus::RANGE_ERROR,
            context+": phase enum is outside the declared domain");
    seed.expectedPhase=static_cast<RestartPhaseAssertion>(rawPhase);
    return true;
}

bool decodePayload(
    std::string_view bytes,
    RestartPayloadCandidate& payload,
    DecodeError& error,
    std::size_t index)
{
    const std::string context="payload["+std::to_string(index)+"]";
    RecordView record;
    if(!parseRecord(bytes,record,error,context)
       ||!exactFields(record,{
            {1U,WireType::U8},{2U,WireType::STRING},
            {3U,WireType::BYTES}},error,context))
        return false;
    std::uint8_t rawOwner=0U;
    if(!decodeU8(field(record,1U),rawOwner,error,context+".owner")
       ||!decodeString(field(record,2U),payload.generationId,error,
            context+".generation")
       ||!decodeBytes(field(record,3U),payload.bytes,error,
            context+".bytes"))
        return false;
    if(rawOwner>=static_cast<std::uint8_t>(RestartOwner::COUNT))
        return error.fail(RestartCodecStatus::RANGE_ERROR,
            context+": owner enum is outside the declared domain");
    payload.owner=static_cast<RestartOwner>(rawOwner);
    payload.consumedBytes=payload.bytes.size();
    return true;
}

RestartCodecEncodeResult encodeFailure(
    RestartCodecStatus status,
    const std::string& detail)
{
    RestartCodecEncodeResult result;
    result.status=status;
    result.detail=detail;
    return result;
}

RestartCodecDecodeResult decodeFailure(
    RestartCodecStatus status,
    const std::string& detail)
{
    RestartCodecDecodeResult result;
    result.status=status;
    result.detail=detail;
    return result;
}

} // namespace

const char* restartCodecStatusName(RestartCodecStatus status)
{
    switch(status)
    {
        case RestartCodecStatus::SUCCESS: return "SUCCESS";
        case RestartCodecStatus::INVALID_INPUT: return "INVALID_INPUT";
        case RestartCodecStatus::INVALID_MAGIC: return "INVALID_MAGIC";
        case RestartCodecStatus::UNSUPPORTED_CODEC_VERSION:
            return "UNSUPPORTED_CODEC_VERSION";
        case RestartCodecStatus::TRUNCATED_INPUT: return "TRUNCATED_INPUT";
        case RestartCodecStatus::TRAILING_DATA: return "TRAILING_DATA";
        case RestartCodecStatus::DUPLICATE_FIELD: return "DUPLICATE_FIELD";
        case RestartCodecStatus::UNKNOWN_FIELD: return "UNKNOWN_FIELD";
        case RestartCodecStatus::MISSING_FIELD: return "MISSING_FIELD";
        case RestartCodecStatus::NONCANONICAL_FIELD_ORDER:
            return "NONCANONICAL_FIELD_ORDER";
        case RestartCodecStatus::TYPE_MISMATCH: return "TYPE_MISMATCH";
        case RestartCodecStatus::RANGE_ERROR: return "RANGE_ERROR";
        case RestartCodecStatus::RESOURCE_LIMIT: return "RESOURCE_LIMIT";
        case RestartCodecStatus::MALFORMED_RECORD:
            return "MALFORMED_RECORD";
        case RestartCodecStatus::MANIFEST_REJECTED:
            return "MANIFEST_REJECTED";
    }
    return "INVALID_INPUT";
}

RestartCodecEncodeResult RestartManifestCodec::encode(
    const RestartCandidatePackage& candidate,
    const RestartManifestExpectation& expected) noexcept
{
    try
    {
        const RestartManifestValidationResult validation=
            RestartManifestValidator::validate(candidate,expected);
        if(!validation.success())
            return encodeFailure(RestartCodecStatus::MANIFEST_REJECTED,
                std::string("manifest validation rejected encode: ")
                    +restartManifestStatusName(validation.status)+": "
                    +validation.detail);

        std::vector<std::string> owners;
        owners.reserve(candidate.manifest.owners.size());
        for(const RestartOwnerRecord& owner:candidate.manifest.owners)
            owners.push_back(encodeOwner(owner));

        std::vector<std::string> seeds;
        seeds.reserve(candidate.manifest.recoverySeeds.size());
        for(const RestartRecoverySeed& seed:candidate.manifest.recoverySeeds)
            seeds.push_back(encodeSeed(seed));

        std::vector<std::string> payloads;
        payloads.reserve(candidate.payloads.size());
        for(const RestartPayloadCandidate& payload:candidate.payloads)
        {
            if(payload.bytes.size()>maximumPayloadBytes)
                return encodeFailure(RestartCodecStatus::RESOURCE_LIMIT,
                    "payload exceeds codec limit");
            payloads.push_back(encodePayload(payload));
        }

        const RestartManifest& manifest=candidate.manifest;
        std::string root=encodedRecord({
            {1U,WireType::STRING,manifest.schema},
            {2U,WireType::U32,scalarU32(manifest.version)},
            {3U,WireType::STRING,manifest.generationId},
            {4U,WireType::STRING,manifest.generationChecksumSha256},
            {5U,WireType::RECORD,encodeGenerationState(
                manifest.generationState)},
            {6U,WireType::U8,scalarU8(
                static_cast<std::uint8_t>(manifest.executionScope))},
            {7U,WireType::U32,scalarU32(manifest.rankCount)},
            {8U,WireType::U32,scalarU32(manifest.rankIndex)},
            {9U,WireType::BOOL,scalarBool(manifest.collectiveCommit)},
            {10U,WireType::STRING,manifest.decompositionId},
            {11U,WireType::U8,scalarU8(
                static_cast<std::uint8_t>(manifest.byteOrder))},
            {12U,WireType::STRING,manifest.floatingPointFormat},
            {13U,WireType::BOOL,scalarBool(manifest.trailingDataPermitted)},
            {14U,WireType::U64,scalarU64(manifest.cellCount)},
            {15U,WireType::RECORD,encodeIdentities(manifest.identities)},
            {16U,WireType::RECORD_LIST,encodedRecordList(owners)},
            {17U,WireType::RECORD_LIST,encodedRecordList(seeds)},
            {18U,WireType::RECORD_LIST,encodedRecordList(payloads)}
        });

        std::string bytes;
        if(root.size()>maximumDocumentBytes-diskHeaderBytes
           ||root.size()>std::numeric_limits<std::size_t>::max()
                -diskHeaderBytes)
            return encodeFailure(RestartCodecStatus::RESOURCE_LIMIT,
                "restart root exceeds codec limit");
        bytes.reserve(diskHeaderBytes+root.size());
        appendRaw(bytes,std::string_view(
            reinterpret_cast<const char*>(diskMagic.data()),diskMagic.size()));
        appendU32(bytes,diskFormatVersion);
        appendU64(bytes,root.size());
        appendRaw(bytes,root);

        RestartCodecEncodeResult result;
        result.status=RestartCodecStatus::SUCCESS;
        result.detail=
            "candidate-only restart package encoded; no file written";
        result.bytes=std::move(bytes);
        return result;
    }
    catch(const std::bad_alloc&)
    {
        return encodeFailure(RestartCodecStatus::RESOURCE_LIMIT,
            "restart codec allocation failed");
    }
    catch(const std::length_error& error)
    {
        return encodeFailure(RestartCodecStatus::RESOURCE_LIMIT,error.what());
    }
    catch(const std::exception& error)
    {
        return encodeFailure(RestartCodecStatus::INVALID_INPUT,error.what());
    }
    catch(...)
    {
        return encodeFailure(RestartCodecStatus::INVALID_INPUT,
            "restart codec encode failed with an unknown exception");
    }
}

RestartCodecDecodeResult RestartManifestCodec::decode(
    const std::string& bytes,
    const RestartManifestExpectation& expected) noexcept
{
    try
    {
        if(bytes.size()>maximumDocumentBytes)
            return decodeFailure(RestartCodecStatus::RESOURCE_LIMIT,
                "restart document exceeds codec limit");
        Reader reader(bytes);
        std::string_view magic;
        if(!reader.readView(diskMagic.size(),magic))
            return decodeFailure(RestartCodecStatus::TRUNCATED_INPUT,
                "restart document is shorter than the magic header");
        if(!std::equal(diskMagic.begin(),diskMagic.end(),
                reinterpret_cast<const unsigned char*>(magic.data())))
            return decodeFailure(RestartCodecStatus::INVALID_MAGIC,
                "restart document magic mismatch");
        std::uint32_t version=0U;
        if(!reader.readU32(version))
            return decodeFailure(RestartCodecStatus::TRUNCATED_INPUT,
                "restart document is missing the codec version");
        if(version!=diskFormatVersion)
            return decodeFailure(
                RestartCodecStatus::UNSUPPORTED_CODEC_VERSION,
                "restart disk codec version is unsupported");
        std::uint64_t rootLength=0U;
        if(!reader.readU64(rootLength))
            return decodeFailure(RestartCodecStatus::TRUNCATED_INPUT,
                "restart document is missing the root length");
        if(rootLength>static_cast<std::uint64_t>(reader.remaining()))
            return decodeFailure(RestartCodecStatus::TRUNCATED_INPUT,
                "restart root record is truncated");
        if(rootLength<static_cast<std::uint64_t>(reader.remaining()))
            return decodeFailure(RestartCodecStatus::TRAILING_DATA,
                "restart document has trailing bytes");
        std::string_view rootBytes;
        if(!reader.readView(rootLength,rootBytes))
            return decodeFailure(RestartCodecStatus::TRUNCATED_INPUT,
                "restart root record is truncated");

        DecodeError error;
        RecordView root;
        if(!parseRecord(rootBytes,root,error,"root")
           ||!exactFields(root,{
                {1U,WireType::STRING},{2U,WireType::U32},
                {3U,WireType::STRING},{4U,WireType::STRING},
                {5U,WireType::RECORD},{6U,WireType::U8},
                {7U,WireType::U32},{8U,WireType::U32},
                {9U,WireType::BOOL},{10U,WireType::STRING},
                {11U,WireType::U8},{12U,WireType::STRING},
                {13U,WireType::BOOL},{14U,WireType::U64},
                {15U,WireType::RECORD},{16U,WireType::RECORD_LIST},
                {17U,WireType::RECORD_LIST},
                {18U,WireType::RECORD_LIST}},error,"root"))
            return decodeFailure(error.status,error.detail);

        RestartCandidatePackage candidate;
        RestartManifest& manifest=candidate.manifest;
        std::uint8_t rawExecutionScope=0U;
        std::uint8_t rawByteOrder=0U;
        if(!decodeString(field(root,1U),manifest.schema,error,"root.schema")
           ||!decodeU32(field(root,2U),manifest.version,error,"root.version")
           ||!decodeString(field(root,3U),manifest.generationId,error,
                "root.generation")
           ||!decodeString(field(root,4U),
                manifest.generationChecksumSha256,error,
                "root.generation-checksum")
           ||!decodeGenerationState(field(root,5U),
                manifest.generationState,error)
           ||!decodeU8(field(root,6U),rawExecutionScope,error,
                "root.execution-scope")
           ||!decodeU32(field(root,7U),manifest.rankCount,error,
                "root.rank-count")
           ||!decodeU32(field(root,8U),manifest.rankIndex,error,
                "root.rank-index")
           ||!decodeBool(field(root,9U),manifest.collectiveCommit,error,
                "root.collective-commit")
           ||!decodeString(field(root,10U),manifest.decompositionId,error,
                "root.decomposition")
           ||!decodeU8(field(root,11U),rawByteOrder,error,
                "root.byte-order")
           ||!decodeString(field(root,12U),manifest.floatingPointFormat,error,
                "root.floating-point")
           ||!decodeBool(field(root,13U),manifest.trailingDataPermitted,error,
                "root.trailing-data-permitted")
           ||!decodeU64(field(root,14U),manifest.cellCount,error,
                "root.cell-count")
           ||!decodeIdentities(field(root,15U),manifest.identities,error))
            return decodeFailure(error.status,error.detail);
        if(rawExecutionScope>static_cast<std::uint8_t>(
                RestartExecutionScope::MPI_PARTITIONED)
           ||rawByteOrder>static_cast<std::uint8_t>(
                RestartByteOrder::BIG_ENDIAN_IEEE754))
            return decodeFailure(RestartCodecStatus::RANGE_ERROR,
                "root enum value is outside the declared domain");
        manifest.executionScope=
            static_cast<RestartExecutionScope>(rawExecutionScope);
        manifest.byteOrder=static_cast<RestartByteOrder>(rawByteOrder);

        std::vector<std::string_view> records;
        if(!parseRecordList(field(root,16U),records,error,"owners"))
            return decodeFailure(error.status,error.detail);
        manifest.owners.resize(records.size());
        for(std::size_t i=0;i<records.size();++i)
            if(!decodeOwner(records[i],manifest.owners[i],error,i))
                return decodeFailure(error.status,error.detail);

        if(!parseRecordList(field(root,17U),records,error,"recovery-seeds"))
            return decodeFailure(error.status,error.detail);
        manifest.recoverySeeds.resize(records.size());
        for(std::size_t i=0;i<records.size();++i)
            if(!decodeSeed(records[i],manifest.recoverySeeds[i],error,i))
                return decodeFailure(error.status,error.detail);

        if(!parseRecordList(field(root,18U),records,error,"payloads"))
            return decodeFailure(error.status,error.detail);
        candidate.payloads.resize(records.size());
        for(std::size_t i=0;i<records.size();++i)
            if(!decodePayload(records[i],candidate.payloads[i],error,i))
                return decodeFailure(error.status,error.detail);
        candidate.manifestTrailingBytes=0U;

        const RestartManifestValidationResult validation=
            RestartManifestValidator::validate(candidate,expected);
        if(!validation.success())
            return decodeFailure(RestartCodecStatus::MANIFEST_REJECTED,
                std::string("decoded manifest rejected: ")
                    +restartManifestStatusName(validation.status)+": "
                    +validation.detail);

        RestartCodecDecodeResult result;
        result.status=RestartCodecStatus::SUCCESS;
        result.detail=
            "candidate-only restart package decoded and validated; no state published";
        result.candidate=std::move(candidate);
        return result;
    }
    catch(const std::bad_alloc&)
    {
        return decodeFailure(RestartCodecStatus::RESOURCE_LIMIT,
            "restart codec allocation failed");
    }
    catch(const std::length_error& error)
    {
        return decodeFailure(RestartCodecStatus::RESOURCE_LIMIT,error.what());
    }
    catch(const std::exception& error)
    {
        return decodeFailure(RestartCodecStatus::MALFORMED_RECORD,error.what());
    }
    catch(...)
    {
        return decodeFailure(RestartCodecStatus::MALFORMED_RECORD,
            "restart codec decode failed with an unknown exception");
    }
}

} // End namespace Foam
