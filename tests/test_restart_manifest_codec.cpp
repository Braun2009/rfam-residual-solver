#include "RestartManifestCodec.H"
#include "ThermoDataProvenance.H"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Foam;

namespace
{

int checks=0;
int failures=0;

void check(bool condition,const std::string& label)
{
    ++checks;
    std::cout << "  [" << (condition?"PASS":"FAIL") << "] "
              << label << '\n';
    if(!condition) ++failures;
}

std::string hash(const std::string& value)
{
    return thermoProvenance::sha256Hex(value);
}

RestartCandidatePackage validCandidate()
{
    RestartCandidatePackage candidate;
    RestartManifest& manifest=candidate.manifest;
    manifest.generationId="generation-20260727-codec-0001";
    manifest.generationState.physicalTime=1.25e-5;
    manifest.generationState.acceptedIteration=42U;
    manifest.generationState.previousAcceptedDtPresent=true;
    manifest.generationState.previousAcceptedDt=1.0e-8;
    manifest.byteOrder=nativeRestartByteOrder();
    manifest.cellCount=2U;
    manifest.identities.configurationSha256=hash("configuration");
    manifest.identities.buildSha256=hash("build");
    manifest.identities.meshSha256=hash("mesh");
    manifest.identities.thermoModelId=
        "peng-robinson-rho-e:sha256="+hash("thermo");
    manifest.identities.kijManifestSha256=hash("kij");
    manifest.identities.mechanismSha256=hash("mechanism");
    manifest.identities.chemistryIntegratorSha256=
        hash("chemistry-integrator");
    manifest.identities.transportSha256=hash("transport");
    manifest.identities.topologyEnvelopeSha256=hash("topology-envelope");

    const auto& roster=restartOwnerRoster();
    manifest.owners.resize(roster.size());
    for(std::size_t i=0;i<roster.size();++i)
    {
        RestartOwnerRecord& owner=manifest.owners[i];
        owner.owner=roster[i];
        owner.publicationOrder=static_cast<std::uint32_t>(i);
        owner.generationId=manifest.generationId;
        owner.identitySha256=hash(
            std::string("owner-")+restartOwnerName(owner.owner));
        owner.payloadSha256=hash(std::string());
    }
    manifest.owners[0].disposition=RestartOwnerDisposition::STATE_PAYLOAD;
    manifest.owners[1].disposition=RestartOwnerDisposition::STATE_PAYLOAD;
    manifest.owners[2].disposition=RestartOwnerDisposition::DISABLED;
    manifest.owners[3].disposition=RestartOwnerDisposition::DISABLED;
    manifest.owners[4].disposition=RestartOwnerDisposition::IDENTITY_ONLY;
    manifest.owners[5].disposition=RestartOwnerDisposition::IDENTITY_ONLY;
    manifest.owners[6].disposition=RestartOwnerDisposition::IDENTITY_ONLY;

    candidate.payloads={
        {RestartOwner::FLUID,manifest.generationId,
            std::string("fluid\0state\xff",12),0U},
        {RestartOwner::TURBULENCE,manifest.generationId,
            "turbulence-state-v2",0U}
    };
    for(RestartPayloadCandidate& payload:candidate.payloads)
    {
        payload.consumedBytes=payload.bytes.size();
        RestartOwnerRecord& owner=
            manifest.owners[static_cast<std::size_t>(payload.owner)];
        owner.payloadBytes=payload.bytes.size();
        owner.payloadSha256=hash(payload.bytes);
    }

    manifest.recoverySeeds={
        {0U,5.0e6,155.0,RestartPhaseAssertion::SINGLE_LIQUID,
            hash("phase-cell-0")},
        {1U,5.0e6,180.0,RestartPhaseAssertion::TWO_PHASE,
            hash("phase-cell-1-quality-x-y")}
    };
    manifest.generationChecksumSha256=restartManifestChecksum(manifest);
    return candidate;
}

RestartManifestExpectation expectation(
    const RestartCandidatePackage& candidate,
    bool requirePreviousDt=true)
{
    RestartManifestExpectation expected;
    expected.identities=candidate.manifest.identities;
    expected.byteOrder=nativeRestartByteOrder();
    expected.cellCount=candidate.manifest.cellCount;
    expected.requirePreviousAcceptedDt=requirePreviousDt;
    const auto& roster=restartOwnerRoster();
    for(std::size_t i=0;i<roster.size();++i)
    {
        expected.owners[i].owner=roster[i];
        expected.owners[i].disposition=
            candidate.manifest.owners[i].disposition;
        expected.owners[i].identitySha256=
            candidate.manifest.owners[i].identitySha256;
    }
    return expected;
}

void refreshChecksum(RestartCandidatePackage& candidate)
{
    candidate.manifest.generationChecksumSha256=
        restartManifestChecksum(candidate.manifest);
}

std::uint16_t getU16(const std::string& bytes,std::size_t offset)
{
    if(offset+2U>bytes.size()) throw std::out_of_range("getU16");
    std::uint16_t value=0U;
    for(unsigned int shift=0U;shift<16U;shift+=8U)
        value|=static_cast<std::uint16_t>(
            static_cast<unsigned char>(bytes[offset++]))<<shift;
    return value;
}

std::uint32_t getU32(const std::string& bytes,std::size_t offset)
{
    if(offset+4U>bytes.size()) throw std::out_of_range("getU32");
    std::uint32_t value=0U;
    for(unsigned int shift=0U;shift<32U;shift+=8U)
        value|=static_cast<std::uint32_t>(
            static_cast<unsigned char>(bytes[offset++]))<<shift;
    return value;
}

std::uint64_t getU64(const std::string& bytes,std::size_t offset)
{
    if(offset+8U>bytes.size()) throw std::out_of_range("getU64");
    std::uint64_t value=0U;
    for(unsigned int shift=0U;shift<64U;shift+=8U)
        value|=static_cast<std::uint64_t>(
            static_cast<unsigned char>(bytes[offset++]))<<shift;
    return value;
}

void putU16(std::string& bytes,std::size_t offset,std::uint16_t value)
{
    if(offset+2U>bytes.size()) throw std::out_of_range("putU16");
    for(unsigned int shift=0U;shift<16U;shift+=8U)
        bytes[offset++]=static_cast<char>((value>>shift)&0xffU);
}

void putU32(std::string& bytes,std::size_t offset,std::uint32_t value)
{
    if(offset+4U>bytes.size()) throw std::out_of_range("putU32");
    for(unsigned int shift=0U;shift<32U;shift+=8U)
        bytes[offset++]=static_cast<char>((value>>shift)&0xffU);
}

void putU64(std::string& bytes,std::size_t offset,std::uint64_t value)
{
    if(offset+8U>bytes.size()) throw std::out_of_range("putU64");
    for(unsigned int shift=0U;shift<64U;shift+=8U)
        bytes[offset++]=static_cast<char>((value>>shift)&0xffU);
}

std::uint64_t doubleBits(double value)
{
    std::uint64_t bits=0U;
    std::memcpy(&bits,&value,sizeof(bits));
    return bits;
}

struct FieldLocation
{
    std::size_t header=0U;
    std::size_t type=0U;
    std::size_t reserved=0U;
    std::size_t length=0U;
    std::size_t data=0U;
    std::size_t end=0U;
};

FieldLocation locateField(const std::string& record,std::uint16_t wanted)
{
    const std::uint16_t count=getU16(record,0U);
    std::size_t offset=2U;
    for(std::uint16_t i=0U;i<count;++i)
    {
        const std::size_t header=offset;
        const std::uint16_t id=getU16(record,offset);
        const std::uint64_t length=getU64(record,offset+4U);
        const std::size_t data=offset+12U;
        if(length>record.size()-data)
            throw std::out_of_range("field length");
        const std::size_t end=data+static_cast<std::size_t>(length);
        if(id==wanted)
            return {header,offset+2U,offset+3U,offset+4U,data,end};
        offset=end;
    }
    throw std::out_of_range("field id");
}

std::string rootRecord(const std::string& document)
{
    if(document.size()<20U) throw std::out_of_range("root header");
    const std::uint64_t length=getU64(document,12U);
    if(length>document.size()-20U) throw std::out_of_range("root length");
    return document.substr(20U,static_cast<std::size_t>(length));
}

std::string withRootRecord(std::string document,const std::string& root)
{
    document.resize(20U);
    putU64(document,12U,root.size());
    document+=root;
    return document;
}

void replaceFieldData(
    std::string& record,
    std::uint16_t id,
    const std::string& data)
{
    const FieldLocation location=locateField(record,id);
    record.replace(location.data,location.end-location.data,data);
    putU64(record,location.length,data.size());
}

void appendField(
    std::string& record,
    std::uint16_t id,
    std::uint8_t type,
    const std::string& data)
{
    const std::uint16_t count=getU16(record,0U);
    putU16(record,0U,static_cast<std::uint16_t>(count+1U));
    record.push_back(static_cast<char>(id&0xffU));
    record.push_back(static_cast<char>((id>>8U)&0xffU));
    record.push_back(static_cast<char>(type));
    record.push_back('\0');
    const std::size_t lengthOffset=record.size();
    record.resize(record.size()+8U,'\0');
    putU64(record,lengthOffset,data.size());
    record+=data;
}

void appendDuplicateField(std::string& record,std::uint16_t id)
{
    const FieldLocation location=locateField(record,id);
    const std::string copy=record.substr(
        location.header,location.end-location.header);
    const std::uint16_t count=getU16(record,0U);
    putU16(record,0U,static_cast<std::uint16_t>(count+1U));
    record+=copy;
}

void removeField(std::string& record,std::uint16_t id)
{
    const FieldLocation location=locateField(record,id);
    record.erase(location.header,location.end-location.header);
    putU16(record,0U,static_cast<std::uint16_t>(getU16(record,0U)-1U));
}

void swapFields(std::string& record,std::uint16_t first,std::uint16_t second)
{
    const std::uint16_t count=getU16(record,0U);
    std::vector<std::pair<std::uint16_t,std::string>> fields;
    std::size_t offset=2U;
    for(std::uint16_t i=0U;i<count;++i)
    {
        const std::uint16_t id=getU16(record,offset);
        const std::uint64_t length=getU64(record,offset+4U);
        const std::size_t end=offset+12U+static_cast<std::size_t>(length);
        fields.push_back({id,record.substr(offset,end-offset)});
        offset=end;
    }
    const auto a=std::find_if(fields.begin(),fields.end(),[&](const auto& item)
        { return item.first==first; });
    const auto b=std::find_if(fields.begin(),fields.end(),[&](const auto& item)
        { return item.first==second; });
    if(a==fields.end()||b==fields.end()) throw std::out_of_range("swap field");
    std::iter_swap(a,b);
    record.resize(2U);
    for(const auto& item:fields) record+=item.second;
}

std::string mutateRoot(
    const std::string& document,
    const std::function<void(std::string&)>& mutation)
{
    std::string root=rootRecord(document);
    mutation(root);
    return withRootRecord(document,root);
}

std::string mutateNestedRecord(
    const std::string& document,
    std::uint16_t rootField,
    const std::function<void(std::string&)>& mutation)
{
    return mutateRoot(document,[&](std::string& root)
    {
        const FieldLocation location=locateField(root,rootField);
        std::string nested=root.substr(
            location.data,location.end-location.data);
        mutation(nested);
        replaceFieldData(root,rootField,nested);
    });
}

std::string mutateListRecord(
    const std::string& document,
    std::uint16_t rootField,
    std::size_t wanted,
    const std::function<void(std::string&)>& mutation)
{
    return mutateRoot(document,[&](std::string& root)
    {
        const FieldLocation location=locateField(root,rootField);
        std::string list=root.substr(location.data,location.end-location.data);
        const std::uint32_t count=getU32(list,0U);
        if(wanted>=count) throw std::out_of_range("list index");
        std::size_t offset=4U;
        for(std::size_t i=0U;i<count;++i)
        {
            const std::size_t lengthOffset=offset;
            const std::uint64_t length=getU64(list,offset);
            const std::size_t data=offset+8U;
            const std::size_t end=data+static_cast<std::size_t>(length);
            if(i==wanted)
            {
                std::string nested=list.substr(data,end-data);
                mutation(nested);
                list.replace(data,end-data,nested);
                putU64(list,lengthOffset,nested.size());
                break;
            }
            offset=end;
        }
        replaceFieldData(root,rootField,list);
    });
}

void expectDecodeStatus(
    const std::string& bytes,
    const RestartManifestExpectation& expected,
    RestartCodecStatus status,
    const std::string& label)
{
    const RestartCodecDecodeResult result=
        RestartManifestCodec::decode(bytes,expected);
    check(result.status==status,
        label+" expected="+restartCodecStatusName(status)
            +" actual="+restartCodecStatusName(result.status));
    check(!result.candidate.has_value(),label+" exposes no decoded candidate");
}

} // namespace

int main()
{
    try
    {
        const RestartCandidatePackage baseline=validCandidate();
        const RestartManifestExpectation expected=expectation(baseline);
        check(RestartManifestValidator::validate(baseline,expected).success(),
            "baseline manifest validates before encoding");

        const RestartCodecEncodeResult encoded=
            RestartManifestCodec::encode(baseline,expected);
        const std::string canonicalDocumentSha256=hash(encoded.bytes);
        check(encoded.success(),"valid candidate encodes");
        check(canonicalDocumentSha256
                =="91c32b4e53ac7a8cf962c23147ccc69d376c290e4c4ba662a7dd5608eb8a0d16",
            "canonical fixture has the frozen disk-byte identity");
        check(encoded.detail.find("no file written")!=std::string::npos,
            "encode result remains explicitly non-I/O");
        check(encoded.bytes.size()>20U,"encoded document has header and root");
        check(encoded.bytes.substr(0U,7U)=="RFAMRM2",
            "encoded document has the v2 codec magic");
        check(getU32(encoded.bytes,8U)==RestartManifestCodec::diskFormatVersion,
            "encoded document binds the codec version");
        check(getU64(encoded.bytes,12U)==encoded.bytes.size()-20U,
            "encoded root length covers the complete document");

        const RestartCodecEncodeResult repeated=
            RestartManifestCodec::encode(baseline,expected);
        check(repeated.success()&&repeated.bytes==encoded.bytes,
            "repeated encoding is byte deterministic");

        const RestartCodecDecodeResult decoded=
            RestartManifestCodec::decode(encoded.bytes,expected);
        check(decoded.success(),"encoded candidate decodes and validates");
        check(decoded.detail.find("no state published")!=std::string::npos,
            "decode result remains explicitly nonpublication");
        check(decoded.candidate->manifestTrailingBytes==0U,
            "decoder derives zero manifest trailing bytes");
        check(decoded.candidate->payloads.size()==baseline.payloads.size()
                &&decoded.candidate->payloads[0].bytes
                    ==baseline.payloads[0].bytes,
            "binary payload bytes including NUL and 0xff round-trip exactly");
        check(decoded.candidate->payloads[0].consumedBytes
                ==decoded.candidate->payloads[0].bytes.size(),
            "decoder derives payload consumed bytes from the wire record");
        check(decoded.candidate->manifest.generationState.physicalTime
                ==baseline.manifest.generationState.physicalTime
                &&decoded.candidate->manifest.generationState.acceptedIteration
                    ==baseline.manifest.generationState.acceptedIteration
                &&decoded.candidate->manifest.generationState
                    .previousAcceptedDtPresent
                &&decoded.candidate->manifest.generationState.previousAcceptedDt
                    ==baseline.manifest.generationState.previousAcceptedDt,
            "generation-state values round-trip exactly");
        check(RestartManifestValidator::validate(
                *decoded.candidate,expected).success(),
            "decoded candidate remains validator-complete");
        const RestartCodecEncodeResult reencoded=
            RestartManifestCodec::encode(*decoded.candidate,expected);
        check(reencoded.success()&&reencoded.bytes==encoded.bytes,
            "encode-decode-encode round-trip is byte identical");

        RestartCandidatePackage absentDt=baseline;
        absentDt.manifest.generationState.previousAcceptedDtPresent=false;
        absentDt.manifest.generationState.previousAcceptedDt=0.0;
        refreshChecksum(absentDt);
        const RestartManifestExpectation absentExpected=
            expectation(absentDt,false);
        const RestartCodecEncodeResult absentEncoded=
            RestartManifestCodec::encode(absentDt,absentExpected);
        const RestartCodecDecodeResult absentDecoded=
            RestartManifestCodec::decode(absentEncoded.bytes,absentExpected);
        check(absentEncoded.success()&&absentDecoded.success()
                &&!absentDecoded.candidate->manifest.generationState
                    .previousAcceptedDtPresent
                &&absentDecoded.candidate->manifest.generationState
                    .previousAcceptedDt==0.0,
            "canonical absent previous-dt state round-trips");
        check(RestartManifestCodec::encode(
                *absentDecoded.candidate,absentExpected).bytes
                    ==absentEncoded.bytes,
            "absent previous-dt encoding is deterministic");

        RestartCandidatePackage invalid=baseline;
        invalid.payloads[0].consumedBytes-=1U;
        const RestartCodecEncodeResult invalidEncode=
            RestartManifestCodec::encode(invalid,expected);
        check(invalidEncode.status==RestartCodecStatus::MANIFEST_REJECTED
                &&invalidEncode.bytes.empty(),
            "encoder refuses a manifest-invalid candidate");

        invalid=baseline;
        invalid.manifest.generationState.previousAcceptedDt=2.0e-8;
        const RestartCodecEncodeResult staleGeneration=
            RestartManifestCodec::encode(invalid,expected);
        check(staleGeneration.status==RestartCodecStatus::MANIFEST_REJECTED,
            "encoder rejects stale generation checksum after dt change");

        RestartCandidatePackage oversizedString=baseline;
        oversizedString.manifest.generationId.assign(
            4U*1024U*1024U+1U,'g');
        for(RestartOwnerRecord& owner:oversizedString.manifest.owners)
            owner.generationId=oversizedString.manifest.generationId;
        for(RestartPayloadCandidate& payload:oversizedString.payloads)
            payload.generationId=oversizedString.manifest.generationId;
        refreshChecksum(oversizedString);
        const RestartManifestExpectation oversizedExpected=
            expectation(oversizedString);
        const RestartCodecEncodeResult oversizedStringEncode=
            RestartManifestCodec::encode(
                oversizedString,oversizedExpected);
        check(oversizedStringEncode.status
                    ==RestartCodecStatus::MANIFEST_REJECTED
                &&oversizedStringEncode.bytes.empty(),
            "oversized semantic string cannot emit a disk document");

        std::string mutation=encoded.bytes;
        mutation[0]='X';
        expectDecodeStatus(mutation,expected,RestartCodecStatus::INVALID_MAGIC,
            "wrong disk magic");

        mutation=encoded.bytes;
        putU32(mutation,8U,RestartManifestCodec::diskFormatVersion+1U);
        expectDecodeStatus(mutation,expected,
            RestartCodecStatus::UNSUPPORTED_CODEC_VERSION,
            "future codec version");

        mutation=encoded.bytes+std::string("x",1U);
        expectDecodeStatus(mutation,expected,RestartCodecStatus::TRAILING_DATA,
            "document trailing byte");

        bool everyPrefixRejected=true;
        for(std::size_t length=0U;length<encoded.bytes.size();++length)
        {
            const RestartCodecDecodeResult result=RestartManifestCodec::decode(
                encoded.bytes.substr(0U,length),expected);
            if(result.success()||result.candidate.has_value())
            {
                everyPrefixRejected=false;
                break;
            }
        }
        check(everyPrefixRejected,
            "every strict prefix is rejected without a candidate");

        bool everySingleBitFlipRejected=true;
        for(std::size_t offset=0U;offset<encoded.bytes.size()
            &&everySingleBitFlipRejected;++offset)
        {
            for(unsigned int bit=0U;bit<8U;++bit)
            {
                std::string flipped=encoded.bytes;
                flipped[offset]=static_cast<char>(
                    static_cast<unsigned char>(flipped[offset])^(1U<<bit));
                const RestartCodecDecodeResult result=
                    RestartManifestCodec::decode(flipped,expected);
                if(result.success()||result.candidate.has_value())
                {
                    everySingleBitFlipRejected=false;
                    break;
                }
            }
        }
        check(everySingleBitFlipRejected,
            "every single-bit disk mutation fails closed");

        mutation=mutateRoot(encoded.bytes,[](std::string& root)
        { appendDuplicateField(root,1U); });
        expectDecodeStatus(mutation,expected,
            RestartCodecStatus::DUPLICATE_FIELD,"duplicate root field");

        mutation=mutateRoot(encoded.bytes,[](std::string& root)
        { appendField(root,99U,6U,"unknown"); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::UNKNOWN_FIELD,
            "unknown root field");

        mutation=mutateRoot(encoded.bytes,[](std::string& root)
        { removeField(root,4U); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::MISSING_FIELD,
            "missing root field");

        mutation=mutateRoot(encoded.bytes,[](std::string& root)
        { swapFields(root,1U,2U); });
        expectDecodeStatus(mutation,expected,
            RestartCodecStatus::NONCANONICAL_FIELD_ORDER,
            "reordered root fields");

        mutation=mutateRoot(encoded.bytes,[](std::string& root)
        { root[locateField(root,1U).type]=static_cast<char>(3U); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::TYPE_MISMATCH,
            "wrong root field type");

        mutation=mutateRoot(encoded.bytes,[](std::string& root)
        { root[locateField(root,1U).type]=static_cast<char>(255U); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::TYPE_MISMATCH,
            "unknown root wire type");

        mutation=mutateRoot(encoded.bytes,[](std::string& root)
        { root[locateField(root,1U).reserved]=1; });
        expectDecodeStatus(mutation,expected,
            RestartCodecStatus::MALFORMED_RECORD,
            "nonzero reserved field byte");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        { appendDuplicateField(generation,1U); });
        expectDecodeStatus(mutation,expected,
            RestartCodecStatus::DUPLICATE_FIELD,
            "duplicate generation-state field");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        { appendField(generation,99U,3U,std::string(8U,'\0')); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::UNKNOWN_FIELD,
            "unknown generation-state field");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        { removeField(generation,2U); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::MISSING_FIELD,
            "missing generation-state field");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        { generation[locateField(generation,1U).type]=static_cast<char>(3U); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::TYPE_MISMATCH,
            "generation-state type mismatch");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        { swapFields(generation,1U,2U); });
        expectDecodeStatus(mutation,expected,
            RestartCodecStatus::NONCANONICAL_FIELD_ORDER,
            "generation-state field reorder");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        { generation.push_back('x'); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::TRAILING_DATA,
            "generation-state nested trailing byte");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        {
            putU64(generation,locateField(generation,1U).data,
                doubleBits(-1.0));
        });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "negative physical time on disk");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        {
            putU64(generation,locateField(generation,1U).data,
                doubleBits(-0.0));
        });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "negative-zero physical time on disk");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        {
            putU64(generation,locateField(generation,1U).data,
                doubleBits(std::numeric_limits<double>::quiet_NaN()));
        });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "NaN physical time on disk");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        {
            putU64(generation,locateField(generation,2U).data,
                static_cast<std::uint64_t>(std::numeric_limits<int>::max())+1U);
        });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "accepted iteration outside current solver range");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        { generation[locateField(generation,3U).data]=2; });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "previous-dt presence outside boolean range");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        {
            putU64(generation,locateField(generation,4U).data,doubleBits(0.0));
        });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "present previous dt is nonpositive");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        {
            generation[locateField(generation,3U).data]=0;
            putU64(generation,locateField(generation,4U).data,
                doubleBits(-0.0));
        });
        expectDecodeStatus(mutation,absentExpected,
            RestartCodecStatus::RANGE_ERROR,
            "absent previous dt rejects negative zero");

        mutation=mutateNestedRecord(encoded.bytes,5U,[](std::string& generation)
        {
            putU64(generation,locateField(generation,4U).data,
                doubleBits(2.0e-8));
        });
        expectDecodeStatus(mutation,expected,
            RestartCodecStatus::MANIFEST_REJECTED,
            "generation dt mutation is checksum-bound");

        mutation=mutateRoot(encoded.bytes,[](std::string& root)
        { root[locateField(root,6U).data]=static_cast<char>(255U); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "execution-scope enum range");

        mutation=mutateRoot(encoded.bytes,[](std::string& root)
        { root[locateField(root,11U).data]=static_cast<char>(255U); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "byte-order enum range");

        mutation=mutateRoot(encoded.bytes,[](std::string& root)
        {
            const FieldLocation location=locateField(root,17U);
            putU32(root,location.data,std::numeric_limits<std::uint32_t>::max());
        });
        expectDecodeStatus(mutation,expected,
            RestartCodecStatus::RESOURCE_LIMIT,"record-list count limit");

        mutation=mutateListRecord(encoded.bytes,16U,0U,[](std::string& owner)
        { appendField(owner,99U,6U,"unknown"); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::UNKNOWN_FIELD,
            "unknown owner field");

        mutation=mutateListRecord(encoded.bytes,16U,0U,[](std::string& owner)
        { owner[locateField(owner,2U).data]=static_cast<char>(255U); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "owner disposition enum range");

        mutation=mutateListRecord(encoded.bytes,17U,0U,[](std::string& seed)
        { appendDuplicateField(seed,5U); });
        expectDecodeStatus(mutation,expected,
            RestartCodecStatus::DUPLICATE_FIELD,
            "duplicate recovery-seed field");

        mutation=mutateListRecord(encoded.bytes,17U,0U,[](std::string& seed)
        { seed[locateField(seed,4U).data]=static_cast<char>(255U); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "recovery phase enum range");

        mutation=mutateListRecord(encoded.bytes,18U,0U,[](std::string& payload)
        { appendField(payload,99U,7U,"unknown"); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::UNKNOWN_FIELD,
            "unknown payload field");

        mutation=mutateListRecord(encoded.bytes,18U,0U,[](std::string& payload)
        { payload[locateField(payload,1U).data]=static_cast<char>(255U); });
        expectDecodeStatus(mutation,expected,RestartCodecStatus::RANGE_ERROR,
            "payload owner enum range");

        mutation=mutateListRecord(encoded.bytes,18U,0U,[](std::string& payload)
        { payload[locateField(payload,3U).data]^=1; });
        expectDecodeStatus(mutation,expected,
            RestartCodecStatus::MANIFEST_REJECTED,
            "payload corruption reaches checksum fail-close");

        RestartManifestExpectation wrongExpected=expected;
        wrongExpected.identities.configurationSha256=hash("wrong-config");
        expectDecodeStatus(encoded.bytes,wrongExpected,
            RestartCodecStatus::MANIFEST_REJECTED,
            "independent expectation identity mismatch");

        if(failures==0)
        {
            std::cout << "RESTART_MANIFEST_CODEC_CANONICAL_SHA256 "
                      << canonicalDocumentSha256 << '\n';
            std::cout << "RESTART_MANIFEST_CODEC_TEST_PASS checks="
                      << checks << '\n';
            return 0;
        }
        std::cout << "RESTART_MANIFEST_CODEC_TEST_FAIL checks="
                  << checks << " failures=" << failures << '\n';
        return 1;
    }
    catch(const std::exception& error)
    {
        std::cout << "RESTART_MANIFEST_CODEC_TEST_FAIL exception="
                  << error.what() << '\n';
        return 1;
    }
}
