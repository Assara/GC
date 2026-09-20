#pragma once
#include <bit>
#include <array>
#include <cerrno>
#include <cstring>
#include <filesystem>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#if defined(__x86_64__)
#include <nmmintrin.h>
#endif

namespace VectorSpace::serialization {
inline std::uint32_t crc32c_software(std::uint32_t crc,const unsigned char* data,std::size_t size) {
    static const auto table=[] {
        std::array<std::uint32_t,256> t{};
        for(std::uint32_t i=0;i<256;++i){auto c=i;for(int j=0;j<8;++j)c=(c>>1)^((c&1)?0x82f63b78U:0);t[i]=c;}
        return t;
    }();
    for(std::size_t i=0;i<size;++i)crc=(crc>>8)^table[(crc^data[i])&255];return crc;
}
#if defined(__x86_64__)
__attribute__((target("sse4.2"))) inline std::uint32_t crc32c_hardware(std::uint32_t crc,const unsigned char* data,std::size_t size) {
    while(size>=8){std::uint64_t x;std::memcpy(&x,data,8);crc=_mm_crc32_u64(crc,x);data+=8;size-=8;}
    while(size--){crc=_mm_crc32_u8(crc,*data++);}return crc;
}
#endif
inline std::uint32_t crc32c(const void* data,std::size_t size) {
    const auto* bytes=static_cast<const unsigned char*>(data);
#if defined(__x86_64__)
    if(__builtin_cpu_supports("sse4.2"))return ~crc32c_hardware(~0U,bytes,size);
#endif
    return ~crc32c_software(~0U,bytes,size);
}
class mapping {
    int fd_=-1;
    void* data_=MAP_FAILED;
    std::size_t size_=0;
public:
    mapping(const std::filesystem::path& path,std::size_t size=0) {
        fd_=::open(path.c_str(),size ? O_RDWR|O_CREAT|O_TRUNC : O_RDONLY,0600);
        if(fd_<0)throw std::runtime_error("cannot open mmap archive: "+path.string());
        if(size) {
            if(size>std::size_t(std::numeric_limits<off_t>::max()) || ::ftruncate(fd_,size)) {
                ::close(fd_);throw std::runtime_error("cannot size mmap archive");
            }
            // Reserve disk blocks before mapping; detect ENOSPC before a write can SIGBUS.
            const int error=::posix_fallocate(fd_,0,size);
            if(error && error!=EOPNOTSUPP && error!=ENOSYS && error!=EINVAL) {
                ::close(fd_);throw std::runtime_error("cannot reserve mmap archive disk space");
            }
            size_=size;
        } else {
            const auto length=::lseek(fd_,0,SEEK_END);
            if(length<=0){::close(fd_);throw std::runtime_error("empty mmap archive");}
            size_=length;
        }
        data_=::mmap(nullptr,size_,size ? PROT_READ|PROT_WRITE : PROT_READ,MAP_SHARED,fd_,0);
        if(data_==MAP_FAILED){::close(fd_);throw std::runtime_error("cannot map archive");}
    }
    mapping(const mapping&)=delete;
    mapping& operator=(const mapping&)=delete;
    ~mapping(){if(data_!=MAP_FAILED)::munmap(data_,size_);if(fd_>=0)::close(fd_);}
    unsigned char* data(){return static_cast<unsigned char*>(data_);}
    const unsigned char* data() const{return static_cast<const unsigned char*>(data_);}
    std::size_t size() const{return size_;}
    void flush(){if(::msync(data_,size_,MS_SYNC) || ::fsync(fd_))throw std::runtime_error("mmap archive flush failed");}
};
inline void encode_word(unsigned char* out,std::uint64_t value){for(int i=0;i<8;++i)out[i]=value>>(8*i);}
inline std::uint64_t decode_word(const unsigned char* in){std::uint64_t value=0;for(int i=0;i<8;++i)value|=std::uint64_t(in[i])<<(8*i);return value;}

// Small metadata is owned here; large payloads are borrowed until finish().
// The final size is known before mmap, so no remapping or large staging copy.
class archive_output {
    struct piece {std::size_t offset,size;const void* payload;};
    std::filesystem::path path_;
    std::vector<unsigned char> metadata_;
    std::vector<piece> pieces_;
    std::size_t size_=0;
protected:
    void payload(const void* data,std::size_t size) {
        if(size>std::numeric_limits<std::size_t>::max()-size_-8)throw std::length_error("archive too large");
        if(size)pieces_.push_back({0,size,data});size_+=size;
    }
public:
    explicit archive_output(const std::filesystem::path& path):path_(path) {
        word(0x47435041434b3031ULL);word(1);word(std::endian::native==std::endian::little?0:1);word(1); // CRC32C
    }
    void bytes(const void* data,std::size_t size) {
        if(!size)return;
        if(size>std::numeric_limits<std::size_t>::max()-size_-8)throw std::length_error("archive too large");
        const auto offset=metadata_.size();const auto* p=static_cast<const unsigned char*>(data);
        metadata_.insert(metadata_.end(),p,p+size);
        if(!pieces_.empty() && !pieces_.back().payload && pieces_.back().offset+pieces_.back().size==offset)pieces_.back().size+=size;
        else pieces_.push_back({offset,size,nullptr});size_+=size;
    }
    void word(std::uint64_t value){unsigned char data[8];encode_word(data,value);bytes(data,8);}
    void text(const std::string& value){word(value.size());bytes(value.data(),value.size());}
    template<class Function> void section(const std::string& name,Function write) {
        text(name);const auto offset=metadata_.size();word(0);const auto begin=size_;
        write();encode_word(metadata_.data()+offset,size_-begin);
    }
    template<class T> void integers(const std::string& name,std::span<const T> values) {
        static_assert(std::is_integral_v<T>);
        section(name,[&]{word(sizeof(T));word(std::is_signed_v<T>);word(values.size());payload(values.data(),values.size_bytes());});
    }
    template<class K> void block(const std::string& role,std::span<const K> values,std::size_t rows,std::size_t columns) {
        static_assert(std::is_trivially_copyable_v<K>);
        if((columns && rows>values.size()/columns) || rows*columns!=values.size())throw std::invalid_argument("archive block dimensions");
        section(role,[&]{
            text(K::name());word(K::characteristic());word(sizeof(K));word(rows);word(columns);word(1);
            for(std::size_t j=0;j<columns;++j){word(j);word(K::characteristic());word(sizeof(K));word(rows);}
            payload(values.data(),values.size_bytes());
        });
    }
    void finish() {
        const auto temporary=path_.string()+".tmp";
        {mapping file(temporary,size_+8);auto* dest=file.data();
            for(const auto& p:pieces_){std::memcpy(dest,p.payload?p.payload:metadata_.data()+p.offset,p.size);dest+=p.size;}
            encode_word(dest,crc32c(file.data(),size_));file.flush();}
        std::filesystem::rename(temporary,path_);
    }
};
class archive_input {
    mapping file_;
    std::size_t position_=0,end_;
public:
    explicit archive_input(const std::filesystem::path& path):file_(path),end_(file_.size()>=8?file_.size()-8:0) {
        if(file_.size()<40 || decode_word(file_.data()+end_)!=crc32c(file_.data(),end_))
            throw std::runtime_error("truncated or corrupt mmap archive");
        expect(0x47435041434b3031ULL);expect(1);expect(std::endian::native==std::endian::little?0:1);expect(1);
    }
    void bytes(void* data,std::size_t size) {
        if(size>end_-position_)throw std::runtime_error("truncated archive section");
        if(size)std::memcpy(data,file_.data()+position_,size);position_+=size;
    }
    std::uint64_t word(){unsigned char data[8];bytes(data,8);return decode_word(data);}
    void expect(std::uint64_t value){if(word()!=value)throw std::runtime_error("incompatible archive metadata");}
    void text(const std::string& expected) {
        expect(expected.size());if(expected.size()>end_-position_ || std::memcmp(file_.data()+position_,expected.data(),expected.size()))
            throw std::runtime_error("archive field/section mismatch");position_+=expected.size();
    }
    template<class Function> void section(const std::string& name,Function read) {
        text(name);const auto length=word();if(length>end_-position_)throw std::runtime_error("invalid archive section size");
        const auto outer=end_;end_=position_+length;read();
        if(position_!=end_)throw std::runtime_error("unconsumed archive section data");end_=outer;
    }
    template<class T> void integers(const std::string& name,std::span<T> values) {
        static_assert(std::is_integral_v<T>);
        section(name,[&]{expect(sizeof(T));expect(std::is_signed_v<T>);expect(values.size());bytes(values.data(),values.size_bytes());});
    }
    template<class K> void block(const std::string& role,std::span<K> values,std::size_t rows,std::size_t columns) {
        static_assert(std::is_trivially_copyable_v<K>);
        section(role,[&]{
            text(K::name());expect(K::characteristic());expect(sizeof(K));expect(rows);expect(columns);expect(1);
            if((columns && rows>values.size()/columns) || rows*columns!=values.size())throw std::invalid_argument("archive block dimensions");
            for(std::size_t j=0;j<columns;++j){expect(j);expect(K::characteristic());expect(sizeof(K));expect(rows);}
            bytes(values.data(),values.size_bytes());
        });
    }
    void finish(){if(position_!=end_)throw std::runtime_error("trailing archive sections");}
};
}
