#pragma once

#include <array>
#include "GraphGeneration/MappedGraphFile.hpp"

namespace GraphGeneration {

struct GCGraphCounts {
    std::uint64_t total = 0;
    std::uint64_t odd_gc = 0;
    std::uint64_t even_gc = 0;
};

template <typename G>
struct GCGraphRecord {
    G graph;
    bool odd_gc = false;  // Vertex permutation times edge-reversal sign.
    bool even_gc = false; // Odd edges, even vertices.
};

namespace gc_file_detail {
inline constexpr std::array<char, 8> magic{'G', 'C', 'G', 'C', '2', '\r', '\n', '\0'};
inline constexpr std::size_t header_size = 64;
inline constexpr std::uint64_t version = 2; // Odd GC includes edge reversals.

inline void put64(std::byte* to, std::uint64_t value) {
    for (int i = 0; i < 8; ++i) to[i] = std::byte(value >> (8 * i));
}
inline std::uint64_t get64(const std::byte* from) {
    std::uint64_t value = 0;
    for (int i = 0; i < 8; ++i)
        value |= std::uint64_t(std::to_integer<unsigned char>(from[i])) << (8 * i);
    return value;
}

// The output format is separate from the transient and legacy final formats.
class Mapping {
public:
    Mapping(const std::filesystem::path& path, std::size_t create_size = 0) {
        fd_ = ::open(path.c_str(), create_size ? O_RDWR | O_CREAT | O_EXCL : O_RDONLY, 0644);
        if (fd_ == -1) detail::throw_system_error("open " + path.string());
        try {
            if (create_size && ::ftruncate(fd_, static_cast<off_t>(create_size)) == -1)
                detail::throw_system_error("ftruncate GC file");
            struct stat status{};
            if (::fstat(fd_, &status) == -1) detail::throw_system_error("fstat GC file");
            if (status.st_size < static_cast<off_t>(header_size))
                throw std::runtime_error("incomplete GC file header");
            size_ = static_cast<std::size_t>(status.st_size);
            auto* mapping = ::mmap(nullptr, size_, create_size ? PROT_READ | PROT_WRITE : PROT_READ,
                                   create_size ? MAP_SHARED : MAP_PRIVATE, fd_, 0);
            if (mapping == MAP_FAILED) detail::throw_system_error("mmap GC file");
            data_ = static_cast<std::byte*>(mapping);
        } catch (...) { release(); throw; }
    }
    ~Mapping() { release(); }
    Mapping(const Mapping&) = delete;
    Mapping& operator=(const Mapping&) = delete;
    std::byte* data() { return data_; }
    const std::byte* data() const { return data_; }
    std::size_t size() const { return size_; }
    void sync() {
        if (::msync(data_, size_, MS_SYNC) == -1) detail::throw_system_error("msync GC file");
    }
private:
    void release() noexcept {
        if (data_) ::munmap(data_, size_);
        if (fd_ != -1) ::close(fd_);
        data_ = nullptr;
        fd_ = -1;
    }
    int fd_ = -1;
    std::byte* data_ = nullptr;
    std::size_t size_ = 0;
};
} // namespace gc_file_detail

template <typename G>
class MappedGCGraphWriter {
    static_assert(G::N_HAIR == 0 && sizeof(Int) == 1);
    static constexpr std::size_t record_size = G::SIZE + 1;
public:
    MappedGCGraphWriter(const std::filesystem::path& path, std::uint64_t count)
        : mapping_(path, detail::checked_file_size(count, record_size)), expected_(count) {}

    void append(const GCGraphRecord<G>& record) {
        if (finished_ || counts_.total >= expected_) throw std::out_of_range("GC writer is full");
        auto* destination = mapping_.data() + 64 + counts_.total * record_size;
        std::memcpy(destination, record.graph.half_edges.data(), G::SIZE);
        destination[G::SIZE] = std::byte(unsigned(record.odd_gc) | (unsigned(record.even_gc) << 1));
        ++counts_.total;
        counts_.odd_gc += record.odd_gc;
        counts_.even_gc += record.even_gc;
    }

    GCGraphCounts finish() {
        if (counts_.total != expected_) throw std::logic_error("incomplete GC output");
        auto* header = mapping_.data();
        std::memcpy(header, gc_file_detail::magic.data(), 8);
        const std::array<std::uint64_t, 7> fields{gc_file_detail::version, G::N_VERTICES_, G::N_EDGES_,
            counts_.total, counts_.odd_gc, counts_.even_gc, record_size};
        for (std::size_t i = 0; i < fields.size(); ++i)
            gc_file_detail::put64(header + 8 * (i + 1), fields[i]);
        mapping_.sync();
        finished_ = true;
        return counts_;
    }
private:
    gc_file_detail::Mapping mapping_;
    std::uint64_t expected_;
    GCGraphCounts counts_;
    bool finished_ = false;
};

template <typename G>
class MappedGCGraphReader {
    static_assert(G::N_HAIR == 0 && sizeof(Int) == 1);
    static constexpr std::size_t record_size = G::SIZE + 1;
public:
    explicit MappedGCGraphReader(const std::filesystem::path& path) : mapping_(path) {
        const auto* header = mapping_.data();
        if (std::memcmp(header, gc_file_detail::magic.data(), 8) != 0
            || gc_file_detail::get64(header + 8) != gc_file_detail::version
            || gc_file_detail::get64(header + 16) != G::N_VERTICES_
            || gc_file_detail::get64(header + 24) != G::N_EDGES_
            || gc_file_detail::get64(header + 56) != record_size)
            throw std::runtime_error("GC header does not match the requested graph type");
        counts_ = {gc_file_detail::get64(header + 32), gc_file_detail::get64(header + 40),
                   gc_file_detail::get64(header + 48)};
        if (counts_.odd_gc > counts_.total || counts_.even_gc > counts_.total
            || mapping_.size() != detail::checked_file_size(counts_.total, record_size))
            throw std::runtime_error("invalid GC counts or file size");
    }
    const GCGraphCounts& counts() const noexcept { return counts_; }
    std::uint64_t size() const noexcept { return counts_.total; }
    GCGraphRecord<G> operator[](std::uint64_t index) const {
        if (index >= size()) throw std::out_of_range("GC record index");
        const auto* source = mapping_.data() + 64 + index * record_size;
        const auto flags = std::to_integer<unsigned char>(source[G::SIZE]);
        if (flags > 3) throw std::runtime_error("invalid GC survival flags");
        GCGraphRecord<G> record;
        std::memcpy(record.graph.half_edges.data(), source, G::SIZE);
        record.odd_gc = flags & 1;
        record.even_gc = flags & 2;
        return record;
    }
private:
    gc_file_detail::Mapping mapping_;
    GCGraphCounts counts_;
};

} // namespace GraphGeneration
