#pragma once

#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "types.hpp"

namespace GraphGeneration {

struct alignas(64) MappedGraphFileHeader {
	std::uint64_t vertices = 0;
	std::uint64_t edges = 0;
	std::uint64_t graph_count = 0;
	std::uint64_t splittable_graph_count = 0;
	std::uint64_t format_version = 0;
	std::uint64_t record_size_bytes = 0;
	std::uint64_t payload_kind = 0;
	std::uint64_t flags = 0;
};

static_assert(sizeof(MappedGraphFileHeader) == 64);

inline constexpr std::uint64_t MAPPED_GRAPH_FORMAT_VERSION = 2;

enum class MappedGraphPayloadKind : std::uint64_t {
	unspecified = 0,
	transient_graph = 1,
	admissible_graph = 2,
	final_graph_with_metadata = 3,
	support_transient_graph = 4,
	vertex_irreducible_support_transient_graph = 5,
	triconnected_support_transient_graph = 6,
	repairable_support_transient_graph = 7,
	ordered_root_split_support_transient_graph = 8,
	k3_ordered_root_split_support_transient_graph = 9,
	unrooted_support_transient_graph = 10,
	unrooted_support_leaf_transient_graph = 11,
	ordered_unrooted_support_leaf_transient_graph = 12
};

namespace detail {

[[noreturn]] inline void throw_system_error(const std::string& operation) {
	throw std::system_error(errno, std::generic_category(), operation);
}

inline std::size_t checked_file_size(std::uint64_t graph_count, std::size_t graph_size) {
	constexpr std::size_t header_size = sizeof(MappedGraphFileHeader);
	if (graph_count > (std::numeric_limits<std::size_t>::max() - header_size) / graph_size) {
		throw std::overflow_error("mapped graph file is too large");
	}
	return header_size + static_cast<std::size_t>(graph_count) * graph_size;
}

} // namespace detail

inline MappedGraphFileHeader read_mapped_graph_file_header(const std::filesystem::path& path) {
	const int fd = ::open(path.c_str(), O_RDONLY);
	if (fd == -1) {
		detail::throw_system_error("open " + path.string());
	}

	struct stat status {};
	if (::fstat(fd, &status) == -1) {
		const int saved_errno = errno;
		::close(fd);
		errno = saved_errno;
		detail::throw_system_error("fstat " + path.string());
	}
	if (status.st_size < static_cast<off_t>(sizeof(MappedGraphFileHeader))) {
		::close(fd);
		throw std::runtime_error("mapped graph file has no complete header: " + path.string());
	}

	MappedGraphFileHeader header;
	const ssize_t bytes_read = ::pread(fd, &header, sizeof(header), 0);
	if (bytes_read == -1) {
		const int saved_errno = errno;
		::close(fd);
		errno = saved_errno;
		detail::throw_system_error("pread " + path.string());
	}
	::close(fd);
	if (bytes_read != static_cast<ssize_t>(sizeof(header))) {
		throw std::runtime_error(
			"mapped graph file has no complete header: " + path.string()
		);
	}
	return header;
}

template <typename GraphType>
class MappedGraphReader {
	public:
		explicit MappedGraphReader(const std::filesystem::path& path)
		{
			fd_ = ::open(path.c_str(), O_RDONLY);
			if (fd_ == -1) {
				detail::throw_system_error("open " + path.string());
			}

			try {
				struct stat status {};
				if (::fstat(fd_, &status) == -1) {
					detail::throw_system_error("fstat " + path.string());
				}
				mapped_size_ = static_cast<std::size_t>(status.st_size);
				if (mapped_size_ < sizeof(MappedGraphFileHeader)) {
					throw std::runtime_error("mapped graph file has no complete header: " + path.string());
				}

				mapping_ = static_cast<const std::byte*>(::mmap(
					nullptr,
					mapped_size_,
					PROT_READ,
					MAP_PRIVATE,
					fd_,
					0
				));
				if (mapping_ == reinterpret_cast<const std::byte*>(MAP_FAILED)) {
					mapping_ = nullptr;
					detail::throw_system_error("mmap " + path.string());
				}

				std::memcpy(&header_, mapping_, sizeof(header_));
				if (header_.vertices != GraphType::N_VERTICES_ || header_.edges != GraphType::N_EDGES_) {
					throw std::runtime_error("mapped graph dimensions do not match the requested graph type");
				}
				if (header_.format_version == 0) {
					// Files written before the active-prefix field existed remain readable.
					header_.splittable_graph_count = header_.graph_count;
				}
				if (header_.format_version <= 1) {
					if constexpr (GraphType::N_HAIR != 0) {
						throw std::runtime_error(
							"legacy mapped graph files contain hairless graph records"
						);
					}
					header_.record_size_bytes = GraphType::SIZE * sizeof(Int);
					header_.payload_kind = static_cast<std::uint64_t>(
						MappedGraphPayloadKind::unspecified
					);
				} else if (header_.format_version == MAPPED_GRAPH_FORMAT_VERSION) {
					if (header_.record_size_bytes != GraphType::SIZE * sizeof(Int)) {
						throw std::runtime_error(
							"mapped graph record size does not match the requested graph type"
						);
					}
				} else {
					throw std::runtime_error("unsupported mapped graph format version");
				}
				if (header_.splittable_graph_count > header_.graph_count) {
					throw std::runtime_error("splittable graph count exceeds total graph count");
				}
				const std::size_t required_size = detail::checked_file_size(
					header_.graph_count,
					GraphType::SIZE * sizeof(Int)
				);
				if (mapped_size_ != required_size) {
					throw std::runtime_error(
						"mapped graph file size does not match its declared record layout"
					);
				}

				::madvise(const_cast<std::byte*>(mapping_), mapped_size_, MADV_SEQUENTIAL);
			} catch (...) {
				release();
				throw;
			}
		}

		~MappedGraphReader() {
			release();
		}

		MappedGraphReader(const MappedGraphReader&) = delete;
		MappedGraphReader& operator=(const MappedGraphReader&) = delete;

		std::uint64_t size() const {
			return header_.graph_count;
		}

		std::uint64_t splittable_size() const {
			return header_.splittable_graph_count;
		}

		MappedGraphPayloadKind payload_kind() const noexcept {
			return static_cast<MappedGraphPayloadKind>(header_.payload_kind);
		}

		const MappedGraphFileHeader& header() const noexcept {
			return header_;
		}

		GraphType operator[](std::uint64_t index) const {
			if (index >= size()) {
				throw std::out_of_range("mapped graph index");
			}
			GraphType graph;
			std::memcpy(
				graph.half_edges.data(),
				graph_address(index),
				GraphType::SIZE * sizeof(Int)
			);
			return graph;
		}

	private:
		void release() noexcept {
			if (mapping_ != nullptr) {
				::munmap(const_cast<std::byte*>(mapping_), mapped_size_);
				mapping_ = nullptr;
			}
			if (fd_ != -1) {
				::close(fd_);
				fd_ = -1;
			}
		}

		const std::byte* graph_address(std::uint64_t index) const {
			return mapping_ + sizeof(MappedGraphFileHeader)
				+ static_cast<std::size_t>(index) * GraphType::SIZE * sizeof(Int);
		}

		int fd_ = -1;
		const std::byte* mapping_ = nullptr;
		std::size_t mapped_size_ = 0;
		MappedGraphFileHeader header_{};
};

template <typename GraphType>
class MappedGraphWriter {
	public:
		MappedGraphWriter(
			const std::filesystem::path& path,
			std::uint64_t graph_count,
			std::uint64_t splittable_graph_count,
			MappedGraphPayloadKind payload_kind = MappedGraphPayloadKind::unspecified
		) : graph_count_(graph_count) {
			if (splittable_graph_count > graph_count) {
				throw std::invalid_argument("splittable graph count exceeds total graph count");
			}
			if (path.has_parent_path()) {
				std::filesystem::create_directories(path.parent_path());
			}

			mapped_size_ = detail::checked_file_size(
				graph_count,
				GraphType::SIZE * sizeof(Int)
			);
			fd_ = ::open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
			if (fd_ == -1) {
				detail::throw_system_error("open " + path.string());
			}
			try {
				if (::ftruncate(fd_, static_cast<off_t>(mapped_size_)) == -1) {
					detail::throw_system_error("ftruncate " + path.string());
				}

				mapping_ = static_cast<std::byte*>(::mmap(
					nullptr,
					mapped_size_,
					PROT_READ | PROT_WRITE,
					MAP_SHARED,
					fd_,
					0
				));
				if (mapping_ == reinterpret_cast<std::byte*>(MAP_FAILED)) {
					mapping_ = nullptr;
					detail::throw_system_error("mmap " + path.string());
				}

				MappedGraphFileHeader header{
					.vertices = GraphType::N_VERTICES_,
					.edges = GraphType::N_EDGES_,
					.graph_count = graph_count,
					.splittable_graph_count = splittable_graph_count,
					.format_version = MAPPED_GRAPH_FORMAT_VERSION,
					.record_size_bytes = GraphType::SIZE * sizeof(Int),
					.payload_kind = static_cast<std::uint64_t>(payload_kind)
				};
				std::memcpy(mapping_, &header, sizeof(header));
			} catch (...) {
				release();
				throw;
			}
		}

		~MappedGraphWriter() {
			release();
		}

		MappedGraphWriter(const MappedGraphWriter&) = delete;
		MappedGraphWriter& operator=(const MappedGraphWriter&) = delete;

		std::uint64_t size() const {
			return graph_count_;
		}

		// Concurrent writes are safe when every caller owns a distinct index.
		void write(std::uint64_t index, const GraphType& graph) const {
			if (index >= size()) {
				throw std::out_of_range("mapped graph index");
			}
			std::memcpy(
				graph_address(index),
				graph.half_edges.data(),
				GraphType::SIZE * sizeof(Int)
			);
		}

	private:
		void release() noexcept {
			if (mapping_ != nullptr) {
				::munmap(mapping_, mapped_size_);
				mapping_ = nullptr;
			}
			if (fd_ != -1) {
				::close(fd_);
				fd_ = -1;
			}
		}

		std::byte* graph_address(std::uint64_t index) const {
			return mapping_ + sizeof(MappedGraphFileHeader)
				+ static_cast<std::size_t>(index) * GraphType::SIZE * sizeof(Int);
		}

		std::uint64_t graph_count_ = 0;
		int fd_ = -1;
		std::byte* mapping_ = nullptr;
		std::size_t mapped_size_ = 0;
};

} // namespace GraphGeneration
