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
#include <utility>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "GraphGeneration/MappedGraphFile.hpp"
#include "GraphGeneration/SupportTransientGraph.hpp"

namespace GraphGeneration {

// Keep the unrooted support frontier separate from every earlier rooted,
// one-hair, and unrestricted-support transient format and from final records.
// The payload is deliberately exposed as the common enum type so callers can
// compare it with MappedGraphFileHeader::payload_kind without knowing the
// on-disk integer.
inline constexpr MappedGraphPayloadKind MAPPED_SUPPORT_TRANSIENT_PAYLOAD_KIND
	= MappedGraphPayloadKind::ordered_unrooted_support_leaf_transient_graph;

template <typename SupportGraph>
inline constexpr std::size_t SUPPORT_TRANSIENT_RECORD_SIZE
	= (SupportGraph::SUPPORT_BIT_COUNT + 7) / 8;

namespace detail {

inline std::size_t checked_support_transient_file_size(
	std::uint64_t graph_count,
	std::size_t record_size
) {
	constexpr std::size_t header_size = sizeof(MappedGraphFileHeader);
	if (record_size != 0
		&& graph_count
			> (std::numeric_limits<std::size_t>::max() - header_size)
				/ record_size) {
		throw std::overflow_error("mapped support transient file is too large");
	}
	return header_size + static_cast<std::size_t>(graph_count) * record_size;
}

template <typename SupportGraph>
SupportGraph decode_support_transient_record(const std::byte* source) {
	typename SupportGraph::support_word_array words{};
	for (std::size_t byte_index = 0;
		 byte_index < SUPPORT_TRANSIENT_RECORD_SIZE<SupportGraph>;
		 ++byte_index) {
		const std::size_t bit = byte_index * 8;
		words[bit / 64] |= static_cast<std::uint64_t>(
			std::to_integer<std::uint8_t>(source[byte_index])
		) << (bit % 64);
	}
	return SupportGraph(std::move(words));
}

template <typename SupportGraph>
void encode_support_transient_record(
	std::byte* destination,
	const SupportGraph& graph
) {
	graph.validate();
	constexpr std::size_t record_size
		= SUPPORT_TRANSIENT_RECORD_SIZE<SupportGraph>;
	const auto& words = graph.support_words();
	for (std::size_t byte_index = 0; byte_index < record_size; ++byte_index) {
		const std::size_t bit = byte_index * 8;
		destination[byte_index] = static_cast<std::byte>(
			words[bit / 64] >> (bit % 64)
		);
	}
}

} // namespace detail

template <typename SupportGraph>
class MappedSupportTransientReader {
	public:
		static constexpr std::size_t RECORD_SIZE_BYTES
			= SUPPORT_TRANSIENT_RECORD_SIZE<SupportGraph>;

		explicit MappedSupportTransientReader(const std::filesystem::path& path) {
			fd_ = ::open(path.c_str(), O_RDONLY);
			if (fd_ == -1) {
				detail::throw_system_error("open " + path.string());
			}

			try {
				struct stat status {};
				if (::fstat(fd_, &status) == -1) {
					detail::throw_system_error("fstat " + path.string());
				}
				if (status.st_size
					< static_cast<off_t>(sizeof(MappedGraphFileHeader))) {
					throw std::runtime_error(
						"mapped support transient file has no complete header: "
						+ path.string()
					);
				}
				if (static_cast<std::uintmax_t>(status.st_size)
					> std::numeric_limits<std::size_t>::max()) {
					throw std::overflow_error(
						"mapped support transient file is too large to map"
					);
				}
				mapped_size_ = static_cast<std::size_t>(status.st_size);

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
				if (header_.vertices != SupportGraph::N_VERTICES_
					|| header_.edges != SupportGraph::N_EDGES_
					|| header_.format_version != MAPPED_GRAPH_FORMAT_VERSION
					|| header_.record_size_bytes != RECORD_SIZE_BYTES
					|| header_.payload_kind != static_cast<std::uint64_t>(
						MAPPED_SUPPORT_TRANSIENT_PAYLOAD_KIND
					)
					|| header_.splittable_graph_count != header_.graph_count
					|| header_.flags != 0) {
					throw std::runtime_error(
						"mapped support transient header does not match "
						"the requested record type"
					);
				}

				const std::size_t required_size
					= detail::checked_support_transient_file_size(
						header_.graph_count,
						RECORD_SIZE_BYTES
					);
				if (mapped_size_ != required_size) {
					throw std::runtime_error(
						"mapped support transient file size does not match "
						"its declared record layout"
					);
				}

				::madvise(
					const_cast<std::byte*>(mapping_),
					mapped_size_,
					MADV_SEQUENTIAL
				);
			} catch (...) {
				release();
				throw;
			}
		}

		~MappedSupportTransientReader() {
			release();
		}

		MappedSupportTransientReader(const MappedSupportTransientReader&) = delete;
		MappedSupportTransientReader& operator=(
			const MappedSupportTransientReader&
		) = delete;

		std::uint64_t size() const noexcept {
			return header_.graph_count;
		}

		std::uint64_t splittable_size() const noexcept {
			return header_.splittable_graph_count;
		}

		MappedGraphPayloadKind payload_kind() const noexcept {
			return static_cast<MappedGraphPayloadKind>(header_.payload_kind);
		}

		const MappedGraphFileHeader& header() const noexcept {
			return header_;
		}

		SupportGraph operator[](std::uint64_t index) const {
			if (index >= size()) {
				throw std::out_of_range("mapped support transient index");
			}
			return detail::decode_support_transient_record<SupportGraph>(
				record_address(index)
			);
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

		const std::byte* record_address(std::uint64_t index) const noexcept {
			return mapping_ + sizeof(MappedGraphFileHeader)
				+ static_cast<std::size_t>(index) * RECORD_SIZE_BYTES;
		}

		int fd_ = -1;
		const std::byte* mapping_ = nullptr;
		std::size_t mapped_size_ = 0;
		MappedGraphFileHeader header_{};
};

template <typename SupportGraph>
class MappedSupportTransientWriter {
	public:
		static constexpr std::size_t RECORD_SIZE_BYTES
			= SUPPORT_TRANSIENT_RECORD_SIZE<SupportGraph>;

		MappedSupportTransientWriter(
			const std::filesystem::path& path,
			std::uint64_t graph_count
		) : graph_count_(graph_count) {
			if (path.has_parent_path()) {
				std::filesystem::create_directories(path.parent_path());
			}

			mapped_size_ = detail::checked_support_transient_file_size(
				graph_count,
				RECORD_SIZE_BYTES
			);
			if (mapped_size_
				> static_cast<std::uintmax_t>(std::numeric_limits<off_t>::max())) {
				throw std::overflow_error(
					"mapped support transient file is too large to create"
				);
			}

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

				const MappedGraphFileHeader header{
					.vertices = SupportGraph::N_VERTICES_,
					.edges = SupportGraph::N_EDGES_,
					.graph_count = graph_count,
					.splittable_graph_count = graph_count,
					.format_version = MAPPED_GRAPH_FORMAT_VERSION,
					.record_size_bytes = RECORD_SIZE_BYTES,
					.payload_kind = static_cast<std::uint64_t>(
						MAPPED_SUPPORT_TRANSIENT_PAYLOAD_KIND
					),
					.flags = 0
				};
				std::memcpy(mapping_, &header, sizeof(header));
			} catch (...) {
				release();
				throw;
			}
		}

		~MappedSupportTransientWriter() {
			release();
		}

		MappedSupportTransientWriter(const MappedSupportTransientWriter&) = delete;
		MappedSupportTransientWriter& operator=(
			const MappedSupportTransientWriter&
		) = delete;

		std::uint64_t size() const noexcept {
			return graph_count_;
		}

		// Concurrent writes are safe when every caller owns a distinct index.
		void write(std::uint64_t index, const SupportGraph& graph) const {
			if (index >= size()) {
				throw std::out_of_range("mapped support transient index");
			}
			detail::encode_support_transient_record(
				record_address(index),
				graph
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

		std::byte* record_address(std::uint64_t index) const noexcept {
			return mapping_ + sizeof(MappedGraphFileHeader)
				+ static_cast<std::size_t>(index) * RECORD_SIZE_BYTES;
		}

		std::uint64_t graph_count_ = 0;
		int fd_ = -1;
		std::byte* mapping_ = nullptr;
		std::size_t mapped_size_ = 0;
};

} // namespace GraphGeneration
