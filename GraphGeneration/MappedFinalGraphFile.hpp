#pragma once

#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "GraphGeneration/FinalCanonicalization.hpp"
#include "GraphGeneration/MappedGraphFile.hpp"

namespace GraphGeneration {

template <typename GraphType>
struct final_graph_record {
	static_assert(GraphType::N_HAIR == 0,
		"final graph records require hairless graphs");

	GraphType graph{};
	std::uint64_t automorphism_order = 0;
	final_graph_survival survival = final_graph_survival::none;

	constexpr bool survives_odd_edges() const noexcept {
		return survives(survival, final_graph_survival::odd_edges);
	}

	constexpr bool survives_even_edges_odd_vertices() const noexcept {
		return survives(
			survival,
			final_graph_survival::even_edges_odd_vertices
		);
	}
};

template <typename GraphType>
inline constexpr std::size_t FINAL_GRAPH_RECORD_SIZE
	= GraphType::SIZE * sizeof(Int)
	+ sizeof(std::uint64_t)
	+ sizeof(std::uint8_t);

template <typename GraphType>
class MappedFinalGraphReader {
	public:
		explicit MappedFinalGraphReader(const std::filesystem::path& path) {
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
					throw std::runtime_error(
						"mapped final graph file has no complete header: " + path.string()
					);
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
				if (header_.format_version != MAPPED_GRAPH_FORMAT_VERSION
					|| header_.vertices != GraphType::N_VERTICES_
					|| header_.edges != GraphType::N_EDGES_
					|| header_.record_size_bytes != FINAL_GRAPH_RECORD_SIZE<GraphType>
					|| header_.payload_kind != static_cast<std::uint64_t>(
						MappedGraphPayloadKind::final_graph_with_metadata
					)
					|| header_.splittable_graph_count != 0) {
					throw std::runtime_error(
						"mapped final graph header does not match the requested record type"
					);
				}

				const std::size_t required_size = detail::checked_file_size(
					header_.graph_count,
					FINAL_GRAPH_RECORD_SIZE<GraphType>
				);
				if (mapped_size_ != required_size) {
					throw std::runtime_error(
						"mapped final graph file size does not match its declared layout"
					);
				}
				::madvise(const_cast<std::byte*>(mapping_), mapped_size_, MADV_SEQUENTIAL);
			} catch (...) {
				release();
				throw;
			}
		}

		~MappedFinalGraphReader() {
			release();
		}

		MappedFinalGraphReader(const MappedFinalGraphReader&) = delete;
		MappedFinalGraphReader& operator=(const MappedFinalGraphReader&) = delete;

		std::uint64_t size() const noexcept {
			return header_.graph_count;
		}

		const MappedGraphFileHeader& header() const noexcept {
			return header_;
		}

		final_graph_record<GraphType> operator[](std::uint64_t index) const {
			if (index >= size()) {
				throw std::out_of_range("mapped final graph index");
			}
			const std::byte* address = record_address(index);
			final_graph_record<GraphType> result;
			std::memcpy(
				result.graph.half_edges.data(),
				address,
				GraphType::SIZE * sizeof(Int)
			);
			address += GraphType::SIZE * sizeof(Int);
			std::memcpy(
				&result.automorphism_order,
				address,
				sizeof(result.automorphism_order)
			);
			address += sizeof(result.automorphism_order);
			std::uint8_t survival = 0;
			std::memcpy(&survival, address, sizeof(survival));
			result.survival = static_cast<final_graph_survival>(survival);
			return result;
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
				+ static_cast<std::size_t>(index) * FINAL_GRAPH_RECORD_SIZE<GraphType>;
		}

		int fd_ = -1;
		const std::byte* mapping_ = nullptr;
		std::size_t mapped_size_ = 0;
		MappedGraphFileHeader header_{};
};

template <typename GraphType>
class MappedFinalGraphWriter {
	public:
		MappedFinalGraphWriter(
			const std::filesystem::path& path,
			std::uint64_t graph_count
		) : graph_count_(graph_count) {
			if (path.has_parent_path()) {
				std::filesystem::create_directories(path.parent_path());
			}
			mapped_size_ = detail::checked_file_size(
				graph_count,
				FINAL_GRAPH_RECORD_SIZE<GraphType>
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

				const MappedGraphFileHeader header{
					.vertices = GraphType::N_VERTICES_,
					.edges = GraphType::N_EDGES_,
					.graph_count = graph_count,
					.splittable_graph_count = 0,
					.format_version = MAPPED_GRAPH_FORMAT_VERSION,
					.record_size_bytes = FINAL_GRAPH_RECORD_SIZE<GraphType>,
					.payload_kind = static_cast<std::uint64_t>(
						MappedGraphPayloadKind::final_graph_with_metadata
					)
				};
				std::memcpy(mapping_, &header, sizeof(header));
			} catch (...) {
				release();
				throw;
			}
		}

		~MappedFinalGraphWriter() {
			release();
		}

		MappedFinalGraphWriter(const MappedFinalGraphWriter&) = delete;
		MappedFinalGraphWriter& operator=(const MappedFinalGraphWriter&) = delete;

		std::uint64_t size() const noexcept {
			return graph_count_;
		}

		void write(
			std::uint64_t index,
			const final_graph_record<GraphType>& record
		) const {
			if (index >= size()) {
				throw std::out_of_range("mapped final graph index");
			}
			std::byte* address = record_address(index);
			std::memcpy(
				address,
				record.graph.half_edges.data(),
				GraphType::SIZE * sizeof(Int)
			);
			address += GraphType::SIZE * sizeof(Int);
			std::memcpy(
				address,
				&record.automorphism_order,
				sizeof(record.automorphism_order)
			);
			address += sizeof(record.automorphism_order);
			const std::uint8_t survival
				= static_cast<std::uint8_t>(record.survival);
			std::memcpy(address, &survival, sizeof(survival));
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
				+ static_cast<std::size_t>(index) * FINAL_GRAPH_RECORD_SIZE<GraphType>;
		}

		std::uint64_t graph_count_ = 0;
		int fd_ = -1;
		std::byte* mapping_ = nullptr;
		std::size_t mapped_size_ = 0;
};

} // namespace GraphGeneration
