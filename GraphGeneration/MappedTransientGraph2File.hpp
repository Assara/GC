#pragma once

#include "GraphGeneration/MappedGraphFile.hpp"

namespace GraphGeneration {

// New payload identifier; the legacy mmap implementation remains unchanged.
inline constexpr auto MAPPED_TRANSIENT_GRAPH2_PAYLOAD_KIND
	= static_cast<MappedGraphPayloadKind>(13);

template <typename GraphType>
class MappedTransientGraph2Writer : public MappedGraphWriter<GraphType> {
	static_assert(GraphType::N_HAIR == 0 && GraphType::SIZE > 0);
	static_assert(sizeof(Int) == 1);
public:
	MappedTransientGraph2Writer(const std::filesystem::path& path, std::uint64_t count)
		: MappedGraphWriter<GraphType>(path, count, count,
			MAPPED_TRANSIENT_GRAPH2_PAYLOAD_KIND) {}
};

template <typename GraphType>
class MappedTransientGraph2Reader : public MappedGraphReader<GraphType> {
	static_assert(GraphType::N_HAIR == 0 && GraphType::SIZE > 0);
	static_assert(sizeof(Int) == 1);
public:
	explicit MappedTransientGraph2Reader(const std::filesystem::path& path)
		: MappedGraphReader<GraphType>(path) {
		const auto& header = this->header();
		if (header.format_version != MAPPED_GRAPH_FORMAT_VERSION
			|| this->payload_kind() != MAPPED_TRANSIENT_GRAPH2_PAYLOAD_KIND
			|| header.splittable_graph_count != header.graph_count || header.flags != 0)
			throw std::runtime_error("not a transient_graph2 stage file");
	}
};

} // namespace GraphGeneration
