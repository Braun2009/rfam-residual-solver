#include "RestartGenerationStore.H"

#include "ThermoDataProvenance.H"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <limits>
#include <mutex>
#include <set>
#include <sstream>
#include <system_error>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#include <sys/stat.h>
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/stat.h>
#ifdef __linux__
#include <sys/syscall.h>
#endif
#include <unistd.h>
#endif

namespace Foam
{
namespace
{

namespace fs = std::filesystem;

constexpr std::uint64_t maximumPayloadBytes =
    16ULL*1024ULL*1024ULL*1024ULL;
constexpr std::uint64_t maximumManifestBytes =
    32ULL*1024ULL*1024ULL*1024ULL;
constexpr std::uint64_t maximumOwnerFileBytes = maximumPayloadBytes+256ULL;
constexpr std::size_t ownerFileFixedBytes = 92U;
constexpr char ownerFileMagic[8]={'R','F','A','M','O','G','2','\0'};
constexpr std::uint32_t ownerFileVersion = 1U;

std::mutex generationStoreMutex;

struct FileReadResult
{
    bool ok = false;
    bool notFound = false;
    std::string detail;
    std::string bytes;
};

struct StageResult
{
    RestartGenerationStoreStatus status =
        RestartGenerationStoreStatus::INVALID_INPUT;
    std::string detail;

    bool success() const
    {
        return status == RestartGenerationStoreStatus::SUCCESS;
    }
};

struct OwnerFilesResult
{
    RestartGenerationStoreStatus status =
        RestartGenerationStoreStatus::INVALID_INPUT;
    std::string detail;
    std::optional<RestartCandidatePackage> candidate;

    bool success() const
    {
        return status == RestartGenerationStoreStatus::SUCCESS
            && candidate.has_value();
    }
};

struct OwnerFileDecodeResult
{
    RestartGenerationStoreStatus status =
        RestartGenerationStoreStatus::OWNER_CONTENT_MISMATCH;
    std::string detail;
    std::string generationId;
    std::string payload;

    bool success() const
    {
        return status == RestartGenerationStoreStatus::SUCCESS;
    }
};

RestartGenerationStorePublishResult publishFailure(
    RestartGenerationStoreStatus status,
    const std::string& detail)
{
    RestartGenerationStorePublishResult result;
    result.status=status;
    result.detail=detail;
    return result;
}

RestartGenerationStoreLoadResult loadFailure(
    RestartGenerationStoreStatus status,
    const std::string& detail)
{
    RestartGenerationStoreLoadResult result;
    result.status=status;
    result.detail=detail;
    return result;
}

StageResult stageFailure(
    RestartGenerationStoreStatus status,
    const std::string& detail)
{
    StageResult result;
    result.status=status;
    result.detail=detail;
    return result;
}

OwnerFilesResult ownerFilesFailure(
    RestartGenerationStoreStatus status,
    const std::string& detail)
{
    OwnerFilesResult result;
    result.status=status;
    result.detail=detail;
    return result;
}

bool safeGenerationId(const std::string& generationId)
{
    if(generationId.empty()||generationId.size()>128U) return false;
    const auto alphaNumeric=[](unsigned char c)
    {
        return (c>='a'&&c<='z')||(c>='A'&&c<='Z')
            ||(c>='0'&&c<='9');
    };
    if(!alphaNumeric(static_cast<unsigned char>(generationId.front()))
       ||!alphaNumeric(static_cast<unsigned char>(generationId.back())))
        return false;
    for(unsigned char c:generationId)
    {
        if(!alphaNumeric(c)&&c!='-'&&c!='_'&&c!='.') return false;
    }
    return true;
}

struct DirectoryIdentity
{
    std::uint64_t namespaceId = 0U;
    std::uint64_t objectId = 0U;

    bool operator==(const DirectoryIdentity& other) const
    {
        return namespaceId==other.namespaceId&&objectId==other.objectId;
    }
};

#ifdef __linux__
bool openDirectoryComponentsNoFollow(
    const fs::path& path,
    int& descriptor,
    std::string& detail)
{
    descriptor=-1;
    if(path.empty())
    {
        detail="directory path is empty";
        return false;
    }

    const fs::path relative=path.relative_path();
    for(const fs::path& component:relative)
    {
        if(component=="."||component=="..")
        {
            detail="directory path contains an unsupported alias component: "
                +path.string();
            return false;
        }
    }

    int flags=O_RDONLY;
#ifdef O_DIRECTORY
    flags|=O_DIRECTORY;
#endif
#ifdef O_CLOEXEC
    flags|=O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
    flags|=O_NOFOLLOW;
#endif

    int current=::open(path.is_absolute()?"/":".",flags);
    if(current<0)
    {
        detail="cannot open directory traversal anchor for "+path.string()
            +": "+std::strerror(errno);
        return false;
    }

    for(const fs::path& component:relative)
    {
        if(component.empty()) continue;
        const int next=::openat(current,component.c_str(),flags);
        const int openError=errno;
        const int closeStatus=::close(current);
        if(next<0)
        {
            detail="cannot open no-follow directory component "
                +component.string()+" in "+path.string()+": "
                +std::strerror(openError);
            return false;
        }
        if(closeStatus!=0)
        {
            ::close(next);
            detail="cannot close directory traversal component for "
                +path.string()+": "+std::strerror(errno);
            return false;
        }
        current=next;
    }

    descriptor=current;
    detail="component-wise no-follow directory handle retained";
    return true;
}
#endif

bool readDirectoryIdentity(
    const fs::path& path,
    DirectoryIdentity& identity,
    std::string& detail)
{
#ifdef _WIN32
    HANDLE directory=CreateFileW(
        path.c_str(),
        GENERIC_READ,
        FILE_SHARE_READ|FILE_SHARE_WRITE|FILE_SHARE_DELETE,
        nullptr,
        OPEN_EXISTING,
        FILE_FLAG_BACKUP_SEMANTICS|FILE_FLAG_OPEN_REPARSE_POINT,
        nullptr);
    if(directory==INVALID_HANDLE_VALUE)
    {
        detail="cannot open stable directory handle: "+path.string();
        return false;
    }
    BY_HANDLE_FILE_INFORMATION information{};
    const bool valid=GetFileInformationByHandle(directory,&information)
        &&(information.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY)!=0U
        &&(information.dwFileAttributes&FILE_ATTRIBUTE_REPARSE_POINT)==0U;
    const BOOL closeStatus=CloseHandle(directory);
    if(!valid||!closeStatus)
    {
        detail="path is not a stable non-reparse directory: "+path.string();
        return false;
    }
    identity.namespaceId=information.dwVolumeSerialNumber;
    identity.objectId=
        (static_cast<std::uint64_t>(information.nFileIndexHigh)<<32U)
        |information.nFileIndexLow;
#elif defined(__linux__)
    int descriptor=-1;
    if(!openDirectoryComponentsNoFollow(path,descriptor,detail)) return false;
    struct stat information{};
    const bool valid=::fstat(descriptor,&information)==0
        &&S_ISDIR(information.st_mode);
    const int closeStatus=::close(descriptor);
    if(!valid||closeStatus!=0)
    {
        detail="path is not a stable component-wise no-follow directory: "
            +path.string();
        return false;
    }
    identity.namespaceId=static_cast<std::uint64_t>(information.st_dev);
    identity.objectId=static_cast<std::uint64_t>(information.st_ino);
#else
    int flags=O_RDONLY;
#ifdef O_DIRECTORY
    flags|=O_DIRECTORY;
#endif
#ifdef O_CLOEXEC
    flags|=O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
    flags|=O_NOFOLLOW;
#endif
    const int descriptor=::open(path.c_str(),flags);
    if(descriptor<0)
    {
        detail="cannot open stable directory handle: "+path.string()+": "
            +std::strerror(errno);
        return false;
    }
    struct stat information{};
    const bool valid=::fstat(descriptor,&information)==0
        &&S_ISDIR(information.st_mode);
    const int closeStatus=::close(descriptor);
    if(!valid||closeStatus!=0)
    {
        detail="path is not a stable no-follow directory: "+path.string();
        return false;
    }
    identity.namespaceId=static_cast<std::uint64_t>(information.st_dev);
    identity.objectId=static_cast<std::uint64_t>(information.st_ino);
#endif
    detail="stable directory identity read";
    return true;
}

class StableDirectoryHandle
{
#ifdef _WIN32
    HANDLE handle_ = INVALID_HANDLE_VALUE;
#else
    int descriptor_ = -1;
#endif
    DirectoryIdentity identity_;

public:
    StableDirectoryHandle() = default;
    StableDirectoryHandle(const StableDirectoryHandle&) = delete;
    StableDirectoryHandle& operator=(const StableDirectoryHandle&) = delete;

    ~StableDirectoryHandle()
    {
#ifdef _WIN32
        if(handle_!=INVALID_HANDLE_VALUE) CloseHandle(handle_);
#else
        if(descriptor_>=0) ::close(descriptor_);
#endif
    }

    bool open(const fs::path& path,std::string& detail)
    {
#ifdef _WIN32
        handle_=CreateFileW(
            path.c_str(),
            GENERIC_READ,
            FILE_SHARE_READ|FILE_SHARE_WRITE|FILE_SHARE_DELETE,
            nullptr,
            OPEN_EXISTING,
            FILE_FLAG_BACKUP_SEMANTICS|FILE_FLAG_OPEN_REPARSE_POINT,
            nullptr);
        if(handle_==INVALID_HANDLE_VALUE)
        {
            detail="cannot open stable directory handle: "+path.string();
            return false;
        }
        BY_HANDLE_FILE_INFORMATION information{};
        if(!GetFileInformationByHandle(handle_,&information)
           ||(information.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY)==0U
           ||(information.dwFileAttributes&FILE_ATTRIBUTE_REPARSE_POINT)!=0U)
        {
            detail="path is not a stable non-reparse directory: "
                +path.string();
            CloseHandle(handle_);
            handle_=INVALID_HANDLE_VALUE;
            return false;
        }
        identity_.namespaceId=information.dwVolumeSerialNumber;
        identity_.objectId=
            (static_cast<std::uint64_t>(information.nFileIndexHigh)<<32U)
            |information.nFileIndexLow;
#elif defined(__linux__)
        if(!openDirectoryComponentsNoFollow(path,descriptor_,detail))
            return false;
        struct stat information{};
        if(::fstat(descriptor_,&information)!=0
           ||!S_ISDIR(information.st_mode))
        {
            detail="path is not a stable component-wise no-follow directory: "
                +path.string();
            ::close(descriptor_);
            descriptor_=-1;
            return false;
        }
        identity_.namespaceId=static_cast<std::uint64_t>(information.st_dev);
        identity_.objectId=static_cast<std::uint64_t>(information.st_ino);
#else
        int flags=O_RDONLY;
#ifdef O_DIRECTORY
        flags|=O_DIRECTORY;
#endif
#ifdef O_CLOEXEC
        flags|=O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
        flags|=O_NOFOLLOW;
#endif
        descriptor_=::open(path.c_str(),flags);
        if(descriptor_<0)
        {
            detail="cannot open stable directory handle: "+path.string()+": "
                +std::strerror(errno);
            return false;
        }
        struct stat information{};
        if(::fstat(descriptor_,&information)!=0
           ||!S_ISDIR(information.st_mode))
        {
            detail="path is not a stable no-follow directory: "+path.string();
            ::close(descriptor_);
            descriptor_=-1;
            return false;
        }
        identity_.namespaceId=static_cast<std::uint64_t>(information.st_dev);
        identity_.objectId=static_cast<std::uint64_t>(information.st_ino);
#endif
        detail="stable directory handle retained";
        return true;
    }

    bool openChild(
        const StableDirectoryHandle& parent,
        const fs::path& admittedParent,
        const fs::path& child,
        std::string& detail)
    {
        if(child.empty()||child.has_root_path()
           ||!child.parent_path().empty()||child=="."||child=="..")
        {
            detail="child directory name is not one safe path component";
            return false;
        }
#ifdef __linux__
        (void)admittedParent;
        int flags=O_RDONLY;
#ifdef O_DIRECTORY
        flags|=O_DIRECTORY;
#endif
#ifdef O_CLOEXEC
        flags|=O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
        flags|=O_NOFOLLOW;
#endif
        descriptor_=::openat(parent.descriptor_,child.c_str(),flags);
        if(descriptor_<0)
        {
            detail="cannot open retained-parent child directory "
                +child.string()+": "+std::strerror(errno);
            return false;
        }
        struct stat information{};
        if(::fstat(descriptor_,&information)!=0
           ||!S_ISDIR(information.st_mode))
        {
            detail="retained-parent child is not a stable no-follow "
                "directory: "+child.string();
            ::close(descriptor_);
            descriptor_=-1;
            return false;
        }
        identity_.namespaceId=static_cast<std::uint64_t>(information.st_dev);
        identity_.objectId=static_cast<std::uint64_t>(information.st_ino);
        detail="retained-parent child directory handle retained";
        return true;
#else
        return open(parent.anchoredPath(admittedParent)/child,detail);
#endif
    }

    fs::path anchoredPath(const fs::path& admittedPath) const
    {
#ifdef __linux__
        (void)admittedPath;
        return fs::path("/proc/self/fd")/std::to_string(descriptor_);
#else
        return admittedPath;
#endif
    }

    bool stillNamedBy(const fs::path& path,std::string& detail) const
    {
        DirectoryIdentity current;
        if(!readDirectoryIdentity(path,current,detail)) return false;
        if(!(current==identity_))
        {
            detail="directory path no longer names the retained identity: "
                +path.string();
            return false;
        }
        detail="directory path still names the retained identity";
        return true;
    }

    bool stillNamedByChild(
        const StableDirectoryHandle& parent,
        const fs::path& admittedParent,
        const fs::path& child,
        std::string& detail) const
    {
#ifdef __linux__
        (void)admittedParent;
        if(child.empty()||child.has_root_path()
           ||!child.parent_path().empty()||child=="."||child=="..")
        {
            detail="child directory name is not one safe path component";
            return false;
        }
        int flags=O_RDONLY;
#ifdef O_DIRECTORY
        flags|=O_DIRECTORY;
#endif
#ifdef O_CLOEXEC
        flags|=O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
        flags|=O_NOFOLLOW;
#endif
        const int descriptor=::openat(parent.descriptor_,child.c_str(),flags);
        if(descriptor<0)
        {
            detail="cannot reopen retained-parent child directory "
                +child.string()+": "+std::strerror(errno);
            return false;
        }
        struct stat information{};
        const bool valid=::fstat(descriptor,&information)==0
            &&S_ISDIR(information.st_mode);
        const int closeStatus=::close(descriptor);
        if(!valid||closeStatus!=0)
        {
            detail="retained-parent child identity is unavailable: "
                +child.string();
            return false;
        }
        const DirectoryIdentity current{
            static_cast<std::uint64_t>(information.st_dev),
            static_cast<std::uint64_t>(information.st_ino)};
        if(!(current==identity_))
        {
            detail="retained parent no longer names the retained child: "
                +child.string();
            return false;
        }
        detail="retained parent still names the retained child identity";
        return true;
#else
        return stillNamedBy(
            parent.anchoredPath(admittedParent)/child,detail);
#endif
    }

    bool sync(std::string& detail) const
    {
#ifdef _WIN32
        const BOOL flushStatus=FlushFileBuffers(handle_);
        const DWORD flushError=GetLastError();
        if(!flushStatus)
        {
            detail="FlushFileBuffers failed for retained directory: win32="
                +std::to_string(flushError);
            return false;
        }
#else
        if(::fsync(descriptor_)!=0)
        {
            detail="fsync failed for retained directory: "
                +std::string(std::strerror(errno));
            return false;
        }
#endif
        detail="retained directory synchronized";
        return true;
    }
};

void appendU32(std::string& bytes,std::uint32_t value)
{
    for(unsigned int shift=0U;shift<32U;shift+=8U)
        bytes.push_back(static_cast<char>((value>>shift)&0xffU));
}

void appendU64(std::string& bytes,std::uint64_t value)
{
    for(unsigned int shift=0U;shift<64U;shift+=8U)
        bytes.push_back(static_cast<char>((value>>shift)&0xffU));
}

bool readU32(
    const std::string& bytes,
    std::size_t offset,
    std::uint32_t& value)
{
    if(offset>bytes.size()||bytes.size()-offset<4U) return false;
    value=0U;
    for(unsigned int shift=0U;shift<32U;shift+=8U)
        value|=static_cast<std::uint32_t>(
            static_cast<unsigned char>(bytes[offset++]))<<shift;
    return true;
}

bool readU64(
    const std::string& bytes,
    std::size_t offset,
    std::uint64_t& value)
{
    if(offset>bytes.size()||bytes.size()-offset<8U) return false;
    value=0U;
    for(unsigned int shift=0U;shift<64U;shift+=8U)
        value|=static_cast<std::uint64_t>(
            static_cast<unsigned char>(bytes[offset++]))<<shift;
    return true;
}

std::string encodeOwnerFile(
    const RestartOwnerRecord& owner,
    const RestartPayloadCandidate& payload)
{
    std::string bytes;
    bytes.reserve(ownerFileFixedBytes+payload.generationId.size()
        +payload.bytes.size());
    bytes.append(ownerFileMagic,sizeof(ownerFileMagic));
    appendU32(bytes,ownerFileVersion);
    bytes.push_back(static_cast<char>(owner.owner));
    bytes.append(3U,'\0');
    appendU32(bytes,static_cast<std::uint32_t>(payload.generationId.size()));
    appendU64(bytes,static_cast<std::uint64_t>(payload.bytes.size()));
    bytes.append(owner.payloadSha256);
    bytes.append(payload.generationId);
    bytes.append(payload.bytes);
    return bytes;
}

OwnerFileDecodeResult ownerFileFailure(
    RestartGenerationStoreStatus status,
    const std::string& detail)
{
    OwnerFileDecodeResult result;
    result.status=status;
    result.detail=detail;
    return result;
}

OwnerFileDecodeResult decodeOwnerFile(
    const std::string& bytes,
    const RestartOwnerRecord& expected)
{
    if(bytes.size()<ownerFileFixedBytes)
        return ownerFileFailure(
            RestartGenerationStoreStatus::OWNER_SIZE_MISMATCH,
            "owner envelope is truncated");
    if(!std::equal(
        std::begin(ownerFileMagic),std::end(ownerFileMagic),bytes.begin()))
        return ownerFileFailure(
            RestartGenerationStoreStatus::OWNER_CONTENT_MISMATCH,
            "owner envelope magic mismatch");
    std::uint32_t version=0U;
    std::uint32_t generationBytes=0U;
    std::uint64_t payloadBytes=0U;
    if(!readU32(bytes,8U,version)||version!=ownerFileVersion
       ||static_cast<unsigned char>(bytes[12])
            !=static_cast<std::uint8_t>(expected.owner)
       ||bytes[13]!='\0'||bytes[14]!='\0'||bytes[15]!='\0')
        return ownerFileFailure(
            RestartGenerationStoreStatus::OWNER_CONTENT_MISMATCH,
            "owner envelope version, owner, or reserved bytes mismatch");
    if(!readU32(bytes,16U,generationBytes)
       ||!readU64(bytes,20U,payloadBytes)
       ||generationBytes==0U||generationBytes>128U
       ||payloadBytes>maximumPayloadBytes)
        return ownerFileFailure(
            RestartGenerationStoreStatus::OWNER_SIZE_MISMATCH,
            "owner envelope declares invalid lengths");
    const std::uint64_t variableBytes=
        static_cast<std::uint64_t>(generationBytes)+payloadBytes;
    if(variableBytes>std::numeric_limits<std::size_t>::max()
       ||ownerFileFixedBytes+static_cast<std::size_t>(variableBytes)
            !=bytes.size())
        return ownerFileFailure(
            RestartGenerationStoreStatus::OWNER_SIZE_MISMATCH,
            "owner envelope length does not consume the file exactly");

    const std::string declaredSha=bytes.substr(28U,64U);
    if(!validRestartSha256(declaredSha)
       ||declaredSha!=expected.payloadSha256)
        return ownerFileFailure(
            RestartGenerationStoreStatus::OWNER_CHECKSUM_MISMATCH,
            "owner envelope SHA256 does not match the manifest");
    const std::string generationId=bytes.substr(92U,generationBytes);
    if(!safeGenerationId(generationId)
       ||generationId!=expected.generationId)
        return ownerFileFailure(
            RestartGenerationStoreStatus::OWNER_GENERATION_MISMATCH,
            "owner envelope generation does not match the manifest");
    if(payloadBytes!=expected.payloadBytes)
        return ownerFileFailure(
            RestartGenerationStoreStatus::OWNER_SIZE_MISMATCH,
            "owner envelope payload length does not match the manifest");

    const std::string payload=bytes.substr(
        ownerFileFixedBytes+generationBytes,
        static_cast<std::size_t>(payloadBytes));
    if(thermoProvenance::sha256Hex(payload)!=declaredSha)
        return ownerFileFailure(
            RestartGenerationStoreStatus::OWNER_CHECKSUM_MISMATCH,
            "owner envelope payload SHA256 mismatch");

    OwnerFileDecodeResult result;
    result.status=RestartGenerationStoreStatus::SUCCESS;
    result.detail="owner envelope validated";
    result.generationId=generationId;
    result.payload=payload;
    return result;
}

bool pathExists(const fs::path& path,std::error_code& error)
{
    const fs::file_status status=fs::symlink_status(path,error);
    if(error==std::errc::no_such_file_or_directory)
    {
        error.clear();
        return false;
    }
    return !error&&status.type()!=fs::file_type::not_found;
}

FileReadResult readRegularFile(
    const fs::path& path,
    std::uint64_t maximumBytes)
{
    FileReadResult result;
#ifdef _WIN32
    HANDLE file=CreateFileW(
        path.c_str(),
        GENERIC_READ,
        FILE_SHARE_READ,
        nullptr,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL|FILE_FLAG_OPEN_REPARSE_POINT
            |FILE_FLAG_SEQUENTIAL_SCAN,
        nullptr);
    if(file==INVALID_HANDLE_VALUE)
    {
        const DWORD error=GetLastError();
        result.notFound=error==ERROR_FILE_NOT_FOUND
            ||error==ERROR_PATH_NOT_FOUND;
        result.detail="cannot open stable file handle: "+path.string()
            +": win32="+std::to_string(error);
        return result;
    }
    BY_HANDLE_FILE_INFORMATION before{};
    if(!GetFileInformationByHandle(file,&before)
       ||(before.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY)!=0U
       ||(before.dwFileAttributes&FILE_ATTRIBUTE_REPARSE_POINT)!=0U
       ||before.nNumberOfLinks!=1U)
    {
        result.detail="file is not a regular, single-link, non-reparse file: "
            +path.string();
        CloseHandle(file);
        return result;
    }
    const std::uint64_t fileBytes=
        (static_cast<std::uint64_t>(before.nFileSizeHigh)<<32U)
        |before.nFileSizeLow;
    if(fileBytes>maximumBytes
       ||fileBytes>std::numeric_limits<std::size_t>::max())
    {
        result.detail="file exceeds the store resource limit: "+path.string();
        CloseHandle(file);
        return result;
    }
    result.bytes.resize(static_cast<std::size_t>(fileBytes));
    std::size_t offset=0U;
    while(offset<result.bytes.size())
    {
        const DWORD request=static_cast<DWORD>(std::min<std::size_t>(
            result.bytes.size()-offset,
            static_cast<std::size_t>(std::numeric_limits<DWORD>::max())));
        DWORD received=0U;
        if(!ReadFile(file,&result.bytes[offset],request,&received,nullptr)
           ||received==0U)
        {
            result.detail="short stable-handle read: "+path.string();
            CloseHandle(file);
            result.bytes.clear();
            return result;
        }
        offset+=received;
    }
    char trailing=0;
    DWORD trailingBytes=0U;
    if(!ReadFile(file,&trailing,1U,&trailingBytes,nullptr)
       ||trailingBytes!=0U)
    {
        result.detail="file changed length during stable-handle read: "
            +path.string();
        CloseHandle(file);
        result.bytes.clear();
        return result;
    }
    BY_HANDLE_FILE_INFORMATION after{};
    const bool unchanged=GetFileInformationByHandle(file,&after)
        &&before.dwVolumeSerialNumber==after.dwVolumeSerialNumber
        &&before.nFileIndexHigh==after.nFileIndexHigh
        &&before.nFileIndexLow==after.nFileIndexLow
        &&before.nFileSizeHigh==after.nFileSizeHigh
        &&before.nFileSizeLow==after.nFileSizeLow
        &&before.ftLastWriteTime.dwHighDateTime
            ==after.ftLastWriteTime.dwHighDateTime
        &&before.ftLastWriteTime.dwLowDateTime
            ==after.ftLastWriteTime.dwLowDateTime;
    const BOOL closeStatus=CloseHandle(file);
    if(!unchanged||!closeStatus)
    {
        result.detail="file identity changed during stable-handle read: "
            +path.string();
        result.bytes.clear();
        return result;
    }
#else
    int flags=O_RDONLY;
#ifdef O_CLOEXEC
    flags|=O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
    flags|=O_NOFOLLOW;
#endif
    const int descriptor=::open(path.c_str(),flags);
    if(descriptor<0)
    {
        result.notFound=errno==ENOENT;
        result.detail="cannot open stable file handle: "+path.string()+": "
            +std::strerror(errno);
        return result;
    }
    struct stat before{};
    if(::fstat(descriptor,&before)!=0||!S_ISREG(before.st_mode)
       ||before.st_nlink!=1||before.st_size<0)
    {
        result.detail="file is not a regular, single-link, no-follow file: "
            +path.string();
        ::close(descriptor);
        return result;
    }
    const std::uint64_t fileBytes=static_cast<std::uint64_t>(before.st_size);
    if(fileBytes>maximumBytes
       ||fileBytes>std::numeric_limits<std::size_t>::max())
    {
        result.detail="file exceeds the store resource limit: "+path.string();
        ::close(descriptor);
        return result;
    }
    result.bytes.resize(static_cast<std::size_t>(fileBytes));
    std::size_t offset=0U;
    while(offset<result.bytes.size())
    {
        const ssize_t received=::read(
            descriptor,&result.bytes[offset],result.bytes.size()-offset);
        if(received<=0)
        {
            result.detail="short stable-handle read: "+path.string();
            ::close(descriptor);
            result.bytes.clear();
            return result;
        }
        offset+=static_cast<std::size_t>(received);
    }
    char trailing=0;
    const ssize_t trailingBytes=::read(descriptor,&trailing,1U);
    struct stat after{};
    const bool unchanged=trailingBytes==0
        &&::fstat(descriptor,&after)==0
        &&before.st_dev==after.st_dev
        &&before.st_ino==after.st_ino
        &&before.st_size==after.st_size
        &&before.st_mtime==after.st_mtime
        &&before.st_ctime==after.st_ctime;
    const int closeStatus=::close(descriptor);
    if(!unchanged||closeStatus!=0)
    {
        result.detail="file identity changed during stable-handle read: "
            +path.string();
        result.bytes.clear();
        return result;
    }
#endif
    result.ok=true;
    result.detail="exact stable-handle file read completed";
    return result;
}

bool syncFile(const fs::path& path,std::string& detail)
{
#ifdef _WIN32
    HANDLE file=CreateFileW(
        path.c_str(),
        GENERIC_READ|GENERIC_WRITE,
        FILE_SHARE_READ,
        nullptr,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL|FILE_FLAG_OPEN_REPARSE_POINT,
        nullptr);
    if(file==INVALID_HANDLE_VALUE)
    {
        detail="cannot open staged file for FlushFileBuffers: "+path.string();
        return false;
    }
    BY_HANDLE_FILE_INFORMATION information{};
    if(!GetFileInformationByHandle(file,&information)
       ||(information.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY)!=0U
       ||(information.dwFileAttributes&FILE_ATTRIBUTE_REPARSE_POINT)!=0U
       ||information.nNumberOfLinks!=1U)
    {
        detail="staged file sync target is not a regular single-link file";
        CloseHandle(file);
        return false;
    }
    const BOOL flushStatus=FlushFileBuffers(file);
    const DWORD flushError=GetLastError();
    const BOOL closeStatus=CloseHandle(file);
    if(!flushStatus)
    {
        detail="FlushFileBuffers failed for staged file: "+path.string()
            +": win32="+std::to_string(flushError);
        return false;
    }
    if(!closeStatus)
    {
        detail="file handle close failed after FlushFileBuffers: "
            +path.string();
        return false;
    }
#else
    int flags=O_RDWR;
#ifdef O_CLOEXEC
    flags|=O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
    flags|=O_NOFOLLOW;
#endif
    const int descriptor=::open(path.c_str(),flags);
    if(descriptor<0)
    {
        detail="cannot open staged file for fsync: "+path.string()+": "
            +std::strerror(errno);
        return false;
    }
    struct stat information{};
    if(::fstat(descriptor,&information)!=0
       ||!S_ISREG(information.st_mode)||information.st_nlink!=1)
    {
        detail="staged file sync target is not a regular single-link file";
        ::close(descriptor);
        return false;
    }
    const int syncStatus=::fsync(descriptor);
    const int syncError=errno;
    const int closeStatus=::close(descriptor);
    if(syncStatus!=0)
    {
        detail="fsync failed for staged file: "+path.string()+": "
            +std::strerror(syncError);
        return false;
    }
    if(closeStatus!=0)
    {
        detail="descriptor close failed after fsync: "+path.string();
        return false;
    }
#endif
    detail="staged file synchronized";
    return true;
}

fs::path temporaryPath(const fs::path& finalPath)
{
    fs::path temporary=finalPath;
    temporary += ".tmp";
    return temporary;
}

bool payloadFault(
    const RestartGenerationStoreFault& fault,
    RestartGenerationStoreFaultPoint point,
    RestartOwner owner)
{
    return fault.point==point&&fault.owner==owner;
}

bool manifestFault(
    const RestartGenerationStoreFault& fault,
    RestartGenerationStoreFaultPoint point)
{
    return fault.point==point;
}

StageResult stageFile(
    const fs::path& path,
    const std::string& bytes,
    bool manifest,
    RestartOwner owner,
    const RestartGenerationStoreFault& fault)
{
    const bool openFailure=manifest
        ?manifestFault(fault,RestartGenerationStoreFaultPoint::MANIFEST_OPEN)
        :payloadFault(
            fault,RestartGenerationStoreFaultPoint::PAYLOAD_OPEN,owner);
    if(openFailure)
        return stageFailure(RestartGenerationStoreStatus::STAGE_OPEN_FAILED,
            "injected staged-file open failure");

    int descriptor=-1;
#ifdef _WIN32
    descriptor=_wopen(
        path.c_str(),
        _O_CREAT|_O_EXCL|_O_WRONLY|_O_BINARY,
        _S_IREAD|_S_IWRITE);
#else
    int flags=O_CREAT|O_EXCL|O_WRONLY;
#ifdef O_CLOEXEC
    flags|=O_CLOEXEC;
#endif
#ifdef O_NOFOLLOW
    flags|=O_NOFOLLOW;
#endif
    descriptor=::open(path.c_str(),flags,0600);
#endif
    if(descriptor<0)
        return stageFailure(RestartGenerationStoreStatus::STAGE_OPEN_FAILED,
            "cannot open staged file: "+path.string());
    std::FILE* output=
#ifdef _WIN32
        _fdopen(descriptor,"wb");
#else
        ::fdopen(descriptor,"wb");
#endif
    if(output==nullptr)
    {
#ifdef _WIN32
        _close(descriptor);
#else
        ::close(descriptor);
#endif
        return stageFailure(RestartGenerationStoreStatus::STAGE_OPEN_FAILED,
            "cannot attach a buffered stream to staged file: "+path.string());
    }

    RestartGenerationStoreStatus failureStatus=
        RestartGenerationStoreStatus::SUCCESS;
    std::string failureDetail;
    const bool writeFailure=manifest
        ?manifestFault(fault,RestartGenerationStoreFaultPoint::MANIFEST_WRITE)
        :payloadFault(
            fault,RestartGenerationStoreFaultPoint::PAYLOAD_WRITE,owner);
    if(writeFailure)
    {
        const std::size_t partial=bytes.size()/2U;
        if(partial!=0U&&std::fwrite(bytes.data(),1U,partial,output)!=partial)
            failureDetail="partial injected write also failed at the OS";
        failureStatus=RestartGenerationStoreStatus::STAGE_WRITE_FAILED;
        if(failureDetail.empty())
            failureDetail="injected staged-file write failure";
    }
    else
    {
        if(std::fwrite(bytes.data(),1U,bytes.size(),output)!=bytes.size())
        {
            failureStatus=RestartGenerationStoreStatus::STAGE_WRITE_FAILED;
            failureDetail="staged-file write failed: "+path.string();
        }
    }

    const int flushStatus=std::fflush(output);
    const bool flushFailure=manifest
        ?manifestFault(fault,RestartGenerationStoreFaultPoint::MANIFEST_FLUSH)
        :payloadFault(
            fault,RestartGenerationStoreFaultPoint::PAYLOAD_FLUSH,owner);
    if(failureStatus==RestartGenerationStoreStatus::SUCCESS
       &&(flushFailure||flushStatus!=0))
    {
        failureStatus=RestartGenerationStoreStatus::STAGE_FLUSH_FAILED;
        failureDetail=flushFailure
            ?"injected staged-file flush failure"
            :"staged-file flush failed: "+path.string();
    }

    const int closeStatus=std::fclose(output);
    const bool closeFailure=manifest
        ?manifestFault(fault,RestartGenerationStoreFaultPoint::MANIFEST_CLOSE)
        :payloadFault(
            fault,RestartGenerationStoreFaultPoint::PAYLOAD_CLOSE,owner);
    if(failureStatus==RestartGenerationStoreStatus::SUCCESS
       &&(closeFailure||closeStatus!=0))
    {
        failureStatus=RestartGenerationStoreStatus::STAGE_CLOSE_FAILED;
        failureDetail=closeFailure
            ?"injected staged-file close failure"
            :"staged-file close failed: "+path.string();
    }
    if(failureStatus!=RestartGenerationStoreStatus::SUCCESS)
        return stageFailure(failureStatus,failureDetail);

    const bool fileSyncFailure=manifest
        ?manifestFault(
            fault,RestartGenerationStoreFaultPoint::MANIFEST_FILE_SYNC)
        :payloadFault(
            fault,RestartGenerationStoreFaultPoint::PAYLOAD_FILE_SYNC,owner);
    if(fileSyncFailure)
        return stageFailure(RestartGenerationStoreStatus::FILE_SYNC_FAILED,
            "injected staged-file synchronization failure");
    std::string detail;
    if(!syncFile(path,detail))
        return stageFailure(
            RestartGenerationStoreStatus::FILE_SYNC_FAILED,detail);

    const bool readbackFailure=manifest
        ?manifestFault(
            fault,RestartGenerationStoreFaultPoint::MANIFEST_READBACK)
        :payloadFault(
            fault,RestartGenerationStoreFaultPoint::PAYLOAD_READBACK,owner);
    if(readbackFailure)
        return stageFailure(RestartGenerationStoreStatus::STAGE_READBACK_FAILED,
            "injected staged-file readback failure");
    const FileReadResult readback=readRegularFile(
        path,manifest?maximumManifestBytes:maximumOwnerFileBytes);
    if(!readback.ok||readback.bytes!=bytes)
        return stageFailure(RestartGenerationStoreStatus::STAGE_READBACK_FAILED,
            readback.ok
                ?"staged-file readback differs from written bytes"
                :readback.detail);

    StageResult result;
    result.status=RestartGenerationStoreStatus::SUCCESS;
    result.detail="staged file flushed, closed, synchronized, and verified";
    return result;
}

StageResult renameFile(
    const fs::path& temporary,
    const fs::path& finalPath,
    bool manifest,
    RestartOwner owner,
    const RestartGenerationStoreFault& fault)
{
    if(temporary.parent_path()!=finalPath.parent_path())
        return stageFailure(RestartGenerationStoreStatus::ATOMIC_RENAME_FAILED,
            "temporary and final paths are not in the same directory");
    const bool renameFailure=manifest
        ?manifestFault(fault,RestartGenerationStoreFaultPoint::MANIFEST_RENAME)
        :payloadFault(
            fault,RestartGenerationStoreFaultPoint::PAYLOAD_RENAME,owner);
    if(renameFailure)
        return stageFailure(RestartGenerationStoreStatus::ATOMIC_RENAME_FAILED,
            "injected atomic rename failure");

    std::error_code error;
    if(pathExists(finalPath,error))
        return stageFailure(RestartGenerationStoreStatus::ATOMIC_RENAME_FAILED,
            "final path already exists: "+finalPath.string());
    if(error)
        return stageFailure(RestartGenerationStoreStatus::ATOMIC_RENAME_FAILED,
            "cannot inspect final path: "+error.message());

#ifdef _WIN32
    if(!MoveFileExW(
            temporary.c_str(),finalPath.c_str(),MOVEFILE_WRITE_THROUGH))
        return stageFailure(RestartGenerationStoreStatus::ATOMIC_RENAME_FAILED,
            "same-directory no-replace rename failed: win32="
                +std::to_string(GetLastError()));
#elif defined(__linux__) && defined(SYS_renameat2)
#ifndef RENAME_NOREPLACE
#define RENAME_NOREPLACE (1U<<0U)
#endif
    if(::syscall(
            SYS_renameat2,
            AT_FDCWD,temporary.c_str(),
            AT_FDCWD,finalPath.c_str(),
            RENAME_NOREPLACE)!=0)
        return stageFailure(RestartGenerationStoreStatus::ATOMIC_RENAME_FAILED,
            "same-directory renameat2(RENAME_NOREPLACE) failed: "
                +std::string(std::strerror(errno)));
#else
    return stageFailure(RestartGenerationStoreStatus::ATOMIC_RENAME_FAILED,
        "atomic no-replace rename is unsupported on this platform");
#endif

    StageResult result;
    result.status=RestartGenerationStoreStatus::SUCCESS;
    result.detail="same-directory atomic rename completed";
    return result;
}

const RestartPayloadCandidate* findPayload(
    const RestartCandidatePackage& candidate,
    RestartOwner owner)
{
    for(const RestartPayloadCandidate& payload:candidate.payloads)
        if(payload.owner==owner) return &payload;
    return nullptr;
}

OwnerFilesResult validateOwnerFiles(
    const fs::path& generationDirectory,
    const RestartCandidatePackage& decoded,
    const RestartManifestExpectation& expected)
{
    RestartCandidatePackage external;
    external.manifest=decoded.manifest;
    external.manifestTrailingBytes=0U;

    for(const RestartOwnerRecord& owner:decoded.manifest.owners)
    {
        const fs::path payloadPath=generationDirectory
            /RestartGenerationStore::ownerPayloadFileName(owner.owner);
        if(owner.disposition!=RestartOwnerDisposition::STATE_PAYLOAD)
        {
            const FileReadResult unexpected=readRegularFile(
                payloadPath,maximumOwnerFileBytes);
            if(!unexpected.notFound)
                return ownerFilesFailure(
                    RestartGenerationStoreStatus::OWNER_FILE_UNEXPECTED,
                    std::string("non-state owner file exists: ")
                        +restartOwnerName(owner.owner));
            continue;
        }

        const FileReadResult file=readRegularFile(
            payloadPath,maximumOwnerFileBytes);
        if(!file.ok)
            return ownerFilesFailure(
                RestartGenerationStoreStatus::OWNER_FILE_MISSING,
                std::string("owner file is unavailable: ")
                    +restartOwnerName(owner.owner)+": "+file.detail);
        const OwnerFileDecodeResult ownerFile=decodeOwnerFile(file.bytes,owner);
        if(!ownerFile.success())
            return ownerFilesFailure(ownerFile.status,
                ownerFile.detail+": "+restartOwnerName(owner.owner));

        const RestartPayloadCandidate* embedded=findPayload(decoded,owner.owner);
        if(embedded==nullptr
           ||embedded->generationId!=ownerFile.generationId
           ||embedded->bytes!=ownerFile.payload)
            return ownerFilesFailure(
                RestartGenerationStoreStatus::OWNER_CONTENT_MISMATCH,
                std::string("owner file differs from codec payload: ")
                    +restartOwnerName(owner.owner));

        RestartPayloadCandidate payload;
        payload.owner=owner.owner;
        payload.generationId=ownerFile.generationId;
        payload.bytes=ownerFile.payload;
        payload.consumedBytes=payload.bytes.size();
        external.payloads.push_back(std::move(payload));
    }

    const RestartManifestValidationResult validation=
        RestartManifestValidator::validate(external,expected);
    if(!validation.success())
        return ownerFilesFailure(
            RestartGenerationStoreStatus::PACKAGE_REJECTED,
            std::string("external owner package rejected: ")
                +restartManifestStatusName(validation.status)+": "
                +validation.detail);

    OwnerFilesResult result;
    result.status=RestartGenerationStoreStatus::SUCCESS;
    result.detail="all seven owner records and external payload files validated";
    result.candidate=std::move(external);
    return result;
}

StageResult validateCommittedLayout(
    const fs::path& generationDirectory,
    const RestartCandidatePackage& decoded)
{
    std::set<std::string> finalNames;
    std::set<std::string> temporaryNames;
    finalNames.insert(RestartGenerationStore::manifestFileName);
    temporaryNames.insert(
        std::string(RestartGenerationStore::manifestFileName)+".tmp");
    for(RestartOwner owner:restartOwnerRoster())
    {
        const std::string name=
            RestartGenerationStore::ownerPayloadFileName(owner);
        temporaryNames.insert(name+".tmp");
    }
    for(const RestartOwnerRecord& owner:decoded.manifest.owners)
    {
        if(owner.disposition==RestartOwnerDisposition::STATE_PAYLOAD)
            finalNames.insert(
                RestartGenerationStore::ownerPayloadFileName(owner.owner));
    }

    std::error_code error;
    fs::directory_iterator iterator(generationDirectory,error);
    if(error)
        return stageFailure(RestartGenerationStoreStatus::FILE_LAYOUT_REJECTED,
            "cannot enumerate generation directory: "+error.message());
    const fs::directory_iterator end;
    for(;iterator!=end;iterator.increment(error))
    {
        if(error)
            return stageFailure(
                RestartGenerationStoreStatus::FILE_LAYOUT_REJECTED,
                "generation directory enumeration failed: "+error.message());
        const std::string name=iterator->path().filename().string();
        std::error_code statusError;
        const fs::file_status status=iterator->symlink_status(statusError);
        if(statusError||status.type()!=fs::file_type::regular)
            return stageFailure(
                RestartGenerationStoreStatus::FILE_LAYOUT_REJECTED,
                "generation entry is not a regular non-symlink file: "+name);
        if(finalNames.count(name)==0U&&temporaryNames.count(name)==0U)
            return stageFailure(
                RestartGenerationStoreStatus::FILE_LAYOUT_REJECTED,
                "unexpected generation file: "+name);
    }

    StageResult result;
    result.status=RestartGenerationStoreStatus::SUCCESS;
    result.detail="committed generation file layout is exact";
    return result;
}

StageResult prepareGenerationDirectory(
    const StableDirectoryHandle& storeRoot,
    const fs::path& generationDirectory)
{
    std::error_code error;
    const fs::file_status generationStatus=
        fs::symlink_status(generationDirectory,error);
    const bool generationMissing=
        error==std::errc::no_such_file_or_directory
        ||(!error&&generationStatus.type()==fs::file_type::not_found);
    if(generationMissing) error.clear();
    else if(error)
        return stageFailure(
            RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
            "cannot inspect generation directory: "+error.message());
    if(!generationMissing)
    {
        if(generationStatus.type()!=fs::file_type::directory)
            return stageFailure(
                RestartGenerationStoreStatus::FILE_LAYOUT_REJECTED,
                "generation path is not a non-symlink directory");
        const fs::path manifestPath=
            generationDirectory/RestartGenerationStore::manifestFileName;
        if(pathExists(manifestPath,error))
            return stageFailure(
                RestartGenerationStoreStatus::GENERATION_CONFLICT,
                "generation commit marker already exists");
        if(error)
            return stageFailure(
                RestartGenerationStoreStatus::FILE_LAYOUT_REJECTED,
                "cannot inspect generation commit marker: "+error.message());
        return stageFailure(
            RestartGenerationStoreStatus::INCOMPLETE_GENERATION,
            "an unpublished generation directory is never reused");
    }

    const bool created=fs::create_directory(generationDirectory,error);
    if(error||!created)
        return stageFailure(
            RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
            "cannot create fresh generation directory: "+error.message());
    std::string detail;
    if(!storeRoot.sync(detail))
        return stageFailure(
            RestartGenerationStoreStatus::DIRECTORY_SYNC_FAILED,detail);

    StageResult result;
    result.status=RestartGenerationStoreStatus::SUCCESS;
    result.detail="fresh immutable generation directory created";
    return result;
}

} // End anonymous namespace

const char* restartGenerationStoreStatusName(
    RestartGenerationStoreStatus status)
{
    switch(status)
    {
        case RestartGenerationStoreStatus::SUCCESS: return "SUCCESS";
        case RestartGenerationStoreStatus::INVALID_INPUT:
            return "INVALID_INPUT";
        case RestartGenerationStoreStatus::UNSAFE_GENERATION_ID:
            return "UNSAFE_GENERATION_ID";
        case RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE:
            return "STORE_ROOT_UNAVAILABLE";
        case RestartGenerationStoreStatus::GENERATION_CONFLICT:
            return "GENERATION_CONFLICT";
        case RestartGenerationStoreStatus::INCOMPLETE_GENERATION:
            return "INCOMPLETE_GENERATION";
        case RestartGenerationStoreStatus::CODEC_REJECTED:
            return "CODEC_REJECTED";
        case RestartGenerationStoreStatus::STAGE_OPEN_FAILED:
            return "STAGE_OPEN_FAILED";
        case RestartGenerationStoreStatus::STAGE_WRITE_FAILED:
            return "STAGE_WRITE_FAILED";
        case RestartGenerationStoreStatus::STAGE_FLUSH_FAILED:
            return "STAGE_FLUSH_FAILED";
        case RestartGenerationStoreStatus::STAGE_CLOSE_FAILED:
            return "STAGE_CLOSE_FAILED";
        case RestartGenerationStoreStatus::FILE_SYNC_FAILED:
            return "FILE_SYNC_FAILED";
        case RestartGenerationStoreStatus::STAGE_READBACK_FAILED:
            return "STAGE_READBACK_FAILED";
        case RestartGenerationStoreStatus::ATOMIC_RENAME_FAILED:
            return "ATOMIC_RENAME_FAILED";
        case RestartGenerationStoreStatus::DIRECTORY_SYNC_FAILED:
            return "DIRECTORY_SYNC_FAILED";
        case RestartGenerationStoreStatus::NOT_COMMITTED:
            return "NOT_COMMITTED";
        case RestartGenerationStoreStatus::MANIFEST_READ_FAILED:
            return "MANIFEST_READ_FAILED";
        case RestartGenerationStoreStatus::MANIFEST_REJECTED:
            return "MANIFEST_REJECTED";
        case RestartGenerationStoreStatus::GENERATION_MISMATCH:
            return "GENERATION_MISMATCH";
        case RestartGenerationStoreStatus::OWNER_FILE_MISSING:
            return "OWNER_FILE_MISSING";
        case RestartGenerationStoreStatus::OWNER_FILE_UNEXPECTED:
            return "OWNER_FILE_UNEXPECTED";
        case RestartGenerationStoreStatus::OWNER_GENERATION_MISMATCH:
            return "OWNER_GENERATION_MISMATCH";
        case RestartGenerationStoreStatus::OWNER_SIZE_MISMATCH:
            return "OWNER_SIZE_MISMATCH";
        case RestartGenerationStoreStatus::OWNER_CHECKSUM_MISMATCH:
            return "OWNER_CHECKSUM_MISMATCH";
        case RestartGenerationStoreStatus::OWNER_CONTENT_MISMATCH:
            return "OWNER_CONTENT_MISMATCH";
        case RestartGenerationStoreStatus::FILE_LAYOUT_REJECTED:
            return "FILE_LAYOUT_REJECTED";
        case RestartGenerationStoreStatus::PACKAGE_REJECTED:
            return "PACKAGE_REJECTED";
        case RestartGenerationStoreStatus::STORE_ROOT_CHANGED:
            return "STORE_ROOT_CHANGED";
        case RestartGenerationStoreStatus::COMMIT_STATE_UNCERTAIN:
            return "COMMIT_STATE_UNCERTAIN";
        case RestartGenerationStoreStatus::INTERNAL_ERROR:
            return "INTERNAL_ERROR";
    }
    return "UNKNOWN";
}

fs::path RestartGenerationStore::generationDirectory(
    const fs::path& storeRoot,
    const std::string& generationId)
{
    return storeRoot/(std::string("rfam-generation-")+generationId);
}

bool RestartGenerationStore::validGenerationId(
    const std::string& generationId) noexcept
{
    return safeGenerationId(generationId);
}

std::string RestartGenerationStore::ownerPayloadFileName(RestartOwner owner)
{
    const std::size_t index=static_cast<std::size_t>(owner);
    if(index>=restartOwnerCount) return std::string();
    std::ostringstream name;
    name << std::setw(2) << std::setfill('0') << index << '-'
         << restartOwnerName(owner) << ".payload";
    return name.str();
}

RestartGenerationStorePublishResult RestartGenerationStore::publish(
    const fs::path& storeRoot,
    const RestartCandidatePackage& candidate,
    const RestartManifestExpectation& expected,
    const RestartGenerationStoreFault& fault) noexcept
{
    try
    {
        const std::lock_guard<std::mutex> lock(generationStoreMutex);
        if(storeRoot.empty())
            return publishFailure(
                RestartGenerationStoreStatus::INVALID_INPUT,
                "store root is empty");
        const std::string& generationId=candidate.manifest.generationId;
        if(!safeGenerationId(generationId))
            return publishFailure(
                RestartGenerationStoreStatus::UNSAFE_GENERATION_ID,
                "generation ID is not a portable path token");

        StableDirectoryHandle storeRootHandle;
        std::string rootDetail;
        if(!storeRootHandle.open(storeRoot,rootDetail))
            return publishFailure(
                RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
                "store root must be an existing non-symlink directory: "
                    +rootDetail);

        const RestartCodecEncodeResult encoded=
            RestartManifestCodec::encode(candidate,expected);
        if(!encoded.success())
            return publishFailure(
                RestartGenerationStoreStatus::CODEC_REJECTED,
                std::string("codec rejected candidate: ")
                    +restartCodecStatusName(encoded.status)+": "+encoded.detail);

        const fs::path visibleGenerationPath=
            generationDirectory(storeRoot,generationId);
        const fs::path generationName=visibleGenerationPath.filename();
        const fs::path generationPath=generationDirectory(
            storeRootHandle.anchoredPath(storeRoot),generationId);
        const StageResult prepared=
            prepareGenerationDirectory(storeRootHandle,generationPath);
        if(!prepared.success())
            return publishFailure(prepared.status,prepared.detail);

        StableDirectoryHandle generationHandle;
        std::string generationDetail;
        if(!generationHandle.openChild(
            storeRootHandle,storeRoot,generationName,generationDetail))
            return publishFailure(
                RestartGenerationStoreStatus::FILE_LAYOUT_REJECTED,
                "fresh generation directory identity is unavailable: "
                    +generationDetail);
        const fs::path retainedGenerationPath=
            generationHandle.anchoredPath(generationPath);

        struct PayloadPath
        {
            RestartOwner owner = RestartOwner::FLUID;
            fs::path temporary;
            fs::path finalPath;
        };
        std::vector<PayloadPath> paths;
        for(const RestartOwnerRecord& owner:candidate.manifest.owners)
        {
            if(owner.disposition!=RestartOwnerDisposition::STATE_PAYLOAD)
                continue;
            const RestartPayloadCandidate* payload=
                findPayload(candidate,owner.owner);
            if(payload==nullptr)
                return publishFailure(
                    RestartGenerationStoreStatus::CODEC_REJECTED,
                    "validated codec package lost an owner payload");
            const fs::path finalPath=
                retainedGenerationPath/ownerPayloadFileName(owner.owner);
            const fs::path staged=temporaryPath(finalPath);
            const std::string ownerFile=encodeOwnerFile(owner,*payload);
            const StageResult stage=
                stageFile(staged,ownerFile,false,owner.owner,fault);
            if(!stage.success())
                return publishFailure(stage.status,stage.detail);
            paths.push_back({owner.owner,staged,finalPath});
        }

        for(const PayloadPath& path:paths)
        {
            const StageResult renamed=renameFile(
                path.temporary,path.finalPath,false,path.owner,fault);
            if(!renamed.success())
                return publishFailure(renamed.status,renamed.detail);
            if(payloadFault(
                fault,
                RestartGenerationStoreFaultPoint::PAYLOAD_DIRECTORY_SYNC,
                path.owner))
                return publishFailure(
                    RestartGenerationStoreStatus::DIRECTORY_SYNC_FAILED,
                    "injected payload directory synchronization failure");
            std::string syncDetail;
            if(!generationHandle.sync(syncDetail))
                return publishFailure(
                    RestartGenerationStoreStatus::DIRECTORY_SYNC_FAILED,
                    syncDetail);
        }

        const fs::path manifestPath=retainedGenerationPath/manifestFileName;
        const fs::path stagedManifest=temporaryPath(manifestPath);
        const StageResult manifestStage=
            stageFile(stagedManifest,encoded.bytes,true,RestartOwner::FLUID,fault);
        if(!manifestStage.success())
            return publishFailure(manifestStage.status,manifestStage.detail);

        const RestartCodecDecodeResult decoded=
            RestartManifestCodec::decode(encoded.bytes,expected);
        if(!decoded.success())
            return publishFailure(
                RestartGenerationStoreStatus::MANIFEST_REJECTED,
                "freshly encoded manifest failed exact decode");
        const OwnerFilesResult ownerFiles=validateOwnerFiles(
            retainedGenerationPath,*decoded.candidate,expected);
        if(!ownerFiles.success())
            return publishFailure(ownerFiles.status,ownerFiles.detail);

        const StageResult renamed=renameFile(
            stagedManifest,manifestPath,true,RestartOwner::FLUID,fault);
        if(!renamed.success())
            return publishFailure(renamed.status,renamed.detail);

        const bool injectedDirectoryFailure=
            manifestFault(
                fault,
                RestartGenerationStoreFaultPoint::MANIFEST_DIRECTORY_SYNC)
            ||manifestFault(
                fault,
                RestartGenerationStoreFaultPoint::MANIFEST_ROLLBACK_REMOVE)
            ||manifestFault(
                fault,
                RestartGenerationStoreFaultPoint::
                    MANIFEST_ROLLBACK_DIRECTORY_SYNC);
        std::string syncDetail;
        if(injectedDirectoryFailure
           ||!generationHandle.sync(syncDetail))
        {
            const bool injectedRemoveFailure=manifestFault(
                fault,
                RestartGenerationStoreFaultPoint::MANIFEST_ROLLBACK_REMOVE);
            std::error_code rollbackError;
            const bool removed=!injectedRemoveFailure
                &&fs::remove(manifestPath,rollbackError);
            if(injectedRemoveFailure||rollbackError||!removed)
                return publishFailure(
                    RestartGenerationStoreStatus::COMMIT_STATE_UNCERTAIN,
                    injectedRemoveFailure
                        ?"injected manifest-marker rollback removal failure"
                        :"manifest-marker rollback removal was not confirmed: "
                            +rollbackError.message());

            const bool injectedRollbackSyncFailure=manifestFault(
                fault,
                RestartGenerationStoreFaultPoint::
                    MANIFEST_ROLLBACK_DIRECTORY_SYNC);
            std::string rollbackSyncDetail;
            if(injectedRollbackSyncFailure
               ||!generationHandle.sync(rollbackSyncDetail))
                return publishFailure(
                    RestartGenerationStoreStatus::COMMIT_STATE_UNCERTAIN,
                    injectedRollbackSyncFailure
                        ?"injected rollback directory synchronization failure"
                        :"manifest marker was removed but rollback durability "
                            "is unconfirmed: "+rollbackSyncDetail);

            if(fault.beforeFinalIdentityCheck!=nullptr)
                fault.beforeFinalIdentityCheck(fault.hookContext);
            std::string rollbackIdentityDetail;
            const bool rollbackGenerationStable=
                generationHandle.stillNamedByChild(
                    storeRootHandle,storeRoot,generationName,
                    rollbackIdentityDetail);
            const bool rollbackRootStable=
                storeRootHandle.stillNamedBy(storeRoot,rootDetail);
            if(!rollbackGenerationStable||!rollbackRootStable)
                return publishFailure(
                    RestartGenerationStoreStatus::COMMIT_STATE_UNCERTAIN,
                    "manifest marker rollback completed but the admitted "
                    "root path identity changed before return: "
                        +rollbackIdentityDetail+"; "+rootDetail);
            return publishFailure(
                RestartGenerationStoreStatus::DIRECTORY_SYNC_FAILED,
                injectedDirectoryFailure
                    ?"injected manifest directory synchronization failure"
                    :syncDetail);
        }

        if(fault.beforeFinalIdentityCheck!=nullptr)
            fault.beforeFinalIdentityCheck(fault.hookContext);
        std::string identityDetail;
        const bool generationStable=
            generationHandle.stillNamedByChild(
                storeRootHandle,storeRoot,generationName,identityDetail);
        const bool rootStable=
            storeRootHandle.stillNamedBy(storeRoot,rootDetail);
        if(!generationStable||!rootStable)
            return publishFailure(
                RestartGenerationStoreStatus::COMMIT_STATE_UNCERTAIN,
                "durable manifest was published but the admitted root path "
                "identity changed before return: "+identityDetail+"; "
                    +rootDetail);

        RestartGenerationStorePublishResult result;
        result.status=RestartGenerationStoreStatus::SUCCESS;
        result.detail=
            "owner payloads committed before the durable manifest marker";
        result.manifestPath=visibleGenerationPath/manifestFileName;
        return result;
    }
    catch(const std::exception& error)
    {
        return publishFailure(RestartGenerationStoreStatus::INTERNAL_ERROR,
            std::string("unexpected publish exception: ")+error.what());
    }
    catch(...)
    {
        return publishFailure(RestartGenerationStoreStatus::INTERNAL_ERROR,
            "unexpected non-standard publish exception");
    }
}

RestartGenerationStoreLoadResult RestartGenerationStore::load(
    const fs::path& storeRoot,
    const std::string& generationId,
    const RestartManifestExpectation& expected,
    const RestartGenerationStoreFault& fault) noexcept
{
    try
    {
        const std::lock_guard<std::mutex> lock(generationStoreMutex);
        if(storeRoot.empty())
            return loadFailure(RestartGenerationStoreStatus::INVALID_INPUT,
                "store root is empty");
        if(!safeGenerationId(generationId))
            return loadFailure(
                RestartGenerationStoreStatus::UNSAFE_GENERATION_ID,
                "generation ID is not a portable path token");

        StableDirectoryHandle storeRootHandle;
        std::string rootDetail;
        if(!storeRootHandle.open(storeRoot,rootDetail))
            return loadFailure(
                RestartGenerationStoreStatus::STORE_ROOT_UNAVAILABLE,
                "store root must be an existing non-symlink directory: "
                    +rootDetail);

        const fs::path visibleGenerationPath=
            generationDirectory(storeRoot,generationId);
        const fs::path generationName=visibleGenerationPath.filename();
        const fs::path generationPath=generationDirectory(
            storeRootHandle.anchoredPath(storeRoot),generationId);
        StableDirectoryHandle generationHandle;
        std::string generationDetail;
        if(!generationHandle.openChild(
            storeRootHandle,storeRoot,generationName,generationDetail))
            return loadFailure(
                RestartGenerationStoreStatus::NOT_COMMITTED,
                "generation directory is unavailable: "+generationDetail);
        const fs::path retainedGenerationPath=
            generationHandle.anchoredPath(generationPath);
        const fs::path manifestPath=
            retainedGenerationPath/manifestFileName;

        const FileReadResult manifest=
            readRegularFile(manifestPath,maximumManifestBytes);
        if(!manifest.ok)
        {
            if(manifest.notFound)
                return loadFailure(
                    RestartGenerationStoreStatus::NOT_COMMITTED,
                    "generation has no final manifest commit marker");
            return loadFailure(
                RestartGenerationStoreStatus::MANIFEST_READ_FAILED,
                manifest.detail);
        }
        const RestartCodecDecodeResult decoded=
            RestartManifestCodec::decode(manifest.bytes,expected);
        if(!decoded.success())
            return loadFailure(
                RestartGenerationStoreStatus::MANIFEST_REJECTED,
                std::string("manifest codec rejected disk bytes: ")
                    +restartCodecStatusName(decoded.status)+": "+decoded.detail);
        if(decoded.candidate->manifest.generationId!=generationId)
            return loadFailure(
                RestartGenerationStoreStatus::GENERATION_MISMATCH,
                "manifest generation does not match the requested directory");

        const StageResult layout=
            validateCommittedLayout(retainedGenerationPath,*decoded.candidate);
        if(!layout.success()) return loadFailure(layout.status,layout.detail);
        const OwnerFilesResult ownerFiles=validateOwnerFiles(
            retainedGenerationPath,*decoded.candidate,expected);
        if(!ownerFiles.success())
            return loadFailure(ownerFiles.status,ownerFiles.detail);
        const StageResult finalLayout=
            validateCommittedLayout(retainedGenerationPath,*decoded.candidate);
        if(!finalLayout.success())
            return loadFailure(finalLayout.status,finalLayout.detail);

        if(fault.beforeFinalIdentityCheck!=nullptr)
            fault.beforeFinalIdentityCheck(fault.hookContext);
        const bool generationStable=
            generationHandle.stillNamedByChild(
                storeRootHandle,storeRoot,generationName,generationDetail);
        const bool rootStable=
            storeRootHandle.stillNamedBy(storeRoot,rootDetail);
        if(!generationStable||!rootStable)
            return loadFailure(
                RestartGenerationStoreStatus::STORE_ROOT_CHANGED,
                "retained root or generation identity changed during load: "
                    +generationDetail+"; "+rootDetail);

        RestartGenerationStoreLoadResult result;
        result.status=RestartGenerationStoreStatus::SUCCESS;
        result.detail=
            "committed generation decoded and all seven owners validated";
        result.candidate=std::move(ownerFiles.candidate);
        return result;
    }
    catch(const std::exception& error)
    {
        return loadFailure(RestartGenerationStoreStatus::INTERNAL_ERROR,
            std::string("unexpected load exception: ")+error.what());
    }
    catch(...)
    {
        return loadFailure(RestartGenerationStoreStatus::INTERNAL_ERROR,
            "unexpected non-standard load exception");
    }
}

} // End namespace Foam
