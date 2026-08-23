#include <array>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

#include <unistd.h>

#include "GraphGeneration/MappedSupportTransientFile.hpp"
#include "GraphGeneration/SupportTransientGraph.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
	if (!condition) {
		std::cerr << "FAIL: " << message << '\n';
		++failures;
	}
}

template <typename Function>
void check_throws(Function&& function, const std::string& message) {
	try {
		std::invoke(std::forward<Function>(function));
		check(false, message);
	} catch (const std::exception&) {
		// Expected.
	}
}

class temporary_directory {
	public:
		temporary_directory() {
			std::array<char, 64> pattern{};
			const char text[] = "/tmp/gc-mapped-support-test-XXXXXX";
			static_assert(sizeof(text) <= pattern.size());
			std::memcpy(pattern.data(), text, sizeof(text));
			char* const result = ::mkdtemp(pattern.data());
			if (result == nullptr) {
				throw std::system_error(
					errno, std::generic_category(), "mkdtemp"
				);
			}
			path_ = result;
		}

		~temporary_directory() {
			std::error_code ignored;
			std::filesystem::remove_all(path_, ignored);
		}

		temporary_directory(const temporary_directory&) = delete;
		temporary_directory& operator=(const temporary_directory&) = delete;

		std::filesystem::path file(const std::string& name) const {
			return path_ / name;
		}

	private:
		std::filesystem::path path_;
};

using Support = GraphGeneration::support_transient_graph<12, 20>;
using Reader = GraphGeneration::MappedSupportTransientReader<Support>;
using Writer = GraphGeneration::MappedSupportTransientWriter<Support>;
using Header = GraphGeneration::MappedGraphFileHeader;

constexpr std::size_t packed_record_size = (Support::SUPPORT_BIT_COUNT + 7) / 8;
static_assert(Support::SUPPORT_BIT_COUNT == 66);
static_assert(packed_record_size == 9);

Support make_support(
	std::initializer_list<std::pair<Int, Int>> edges
) {
	Support result;
	for (const auto [first, second] : edges) {
		result.set_support_edge(first, second);
	}
	result.validate();
	return result;
}

std::vector<std::byte> expected_bytes(const Support& graph) {
	std::vector<std::byte> result(packed_record_size);
	std::size_t bit = 0;
	for (Int first = 0; first < Support::N_VERTICES_; ++first) {
		for (Int second = static_cast<Int>(first + 1);
			 second < Support::N_VERTICES_;
			 ++second, ++bit) {
			if (graph.has_support_edge(first, second)) {
				result[bit / 8] |= std::byte{1} << (bit % 8);
			}
		}
	}
	return result;
}

void write_records(
	const std::filesystem::path& path,
	const std::vector<Support>& records
) {
	Writer writer(path, records.size());
	for (std::size_t i = 0; i < records.size(); ++i) {
		writer.write(i, records[i]);
	}
}

template <typename Mutator>
void mutate_header(const std::filesystem::path& path, Mutator&& mutator) {
	Header header = GraphGeneration::read_mapped_graph_file_header(path);
	std::invoke(std::forward<Mutator>(mutator), header);
	std::fstream file(
		path, std::ios::binary | std::ios::in | std::ios::out
	);
	if (!file) {
		throw std::runtime_error("cannot open mapped test file for mutation");
	}
	file.write(reinterpret_cast<const char*>(&header), sizeof(header));
	if (!file) {
		throw std::runtime_error("cannot mutate mapped test header");
	}
}

void overwrite_payload_byte(
	const std::filesystem::path& path,
	std::size_t byte_index,
	std::byte value
) {
	std::fstream file(
		path, std::ios::binary | std::ios::in | std::ios::out
	);
	if (!file) {
		throw std::runtime_error("cannot open mapped test payload for mutation");
	}
	file.seekp(static_cast<std::streamoff>(sizeof(Header) + byte_index));
	const char raw = static_cast<char>(value);
	file.write(&raw, 1);
	if (!file) {
		throw std::runtime_error("cannot mutate mapped test payload");
	}
}

template <typename Mutation>
void malformed_header_case(
	const temporary_directory& directory,
	const std::string& name,
	Mutation&& mutation
) {
	const auto path = directory.file(name + ".gcg");
	write_records(path, {make_support({{0, 1}, {1, 11}})});
	mutate_header(path, std::forward<Mutation>(mutation));
	check_throws([&] { Reader input(path); }, "reject " + name);
}

void test_round_trip_and_packing(const temporary_directory& directory) {
	const std::vector<Support> expected{
		make_support({{0, 1}, {0, 11}, {1, 2}}),
		make_support({{2, 7}, {3, 9}, {10, 11}}),
		make_support({})
	};
	const auto path = directory.file("round-trip.gcg");
	write_records(path, expected);

	const Header header = GraphGeneration::read_mapped_graph_file_header(path);
	check(header.vertices == Support::N_VERTICES_, "header vertex count");
	check(header.edges == Support::N_EDGES_, "header expanded edge count");
	check(header.graph_count == expected.size(), "header graph count");
	check(
		header.splittable_graph_count == expected.size(),
		"all support records are an active frontier"
	);
	check(
		header.format_version == GraphGeneration::MAPPED_GRAPH_FORMAT_VERSION,
		"support file format version"
	);
	check(
		header.record_size_bytes == packed_record_size,
		"support records use packed bytes rather than uint64 word padding"
	);
	check(
		header.payload_kind == static_cast<std::uint64_t>(
			GraphGeneration::MAPPED_SUPPORT_TRANSIENT_PAYLOAD_KIND
		),
		"support transient payload kind"
	);
	check(header.flags == 0, "support file flags are zero");
	check(
		std::filesystem::file_size(path)
			== sizeof(Header) + expected.size() * packed_record_size,
		"support mapped file has the exact packed size"
	);

	Reader input(path);
	check(input.size() == expected.size(), "reader record count");
	check(
		input.splittable_size() == expected.size(),
		"reader active-frontier count"
	);
	check(
		input.payload_kind()
			== GraphGeneration::MAPPED_SUPPORT_TRANSIENT_PAYLOAD_KIND,
		"reader support payload kind"
	);
	for (std::size_t i = 0; i < expected.size(); ++i) {
		check(input[i] == expected[i], "packed support round trip record "
			+ std::to_string(i));
	}
	check_throws(
		[&] { (void)input[expected.size()]; },
		"reader rejects an out-of-range index"
	);

	std::ifstream raw(path, std::ios::binary);
	raw.seekg(static_cast<std::streamoff>(sizeof(Header)));
	for (std::size_t i = 0; i < expected.size(); ++i) {
		std::vector<std::byte> actual(packed_record_size);
		raw.read(
			reinterpret_cast<char*>(actual.data()),
			static_cast<std::streamsize>(actual.size())
		);
		check(actual == expected_bytes(expected[i]),
			"on-disk bit packing record " + std::to_string(i));
	}

	const auto range_path = directory.file("writer-range.gcg");
	Writer writer(range_path, 1);
	check_throws(
		[&] { writer.write(1, expected.front()); },
		"writer rejects an out-of-range index"
	);
}

void test_empty_file(const temporary_directory& directory) {
	const auto path = directory.file("empty.gcg");
	write_records(path, {});
	Reader input(path);
	check(input.size() == 0, "zero-record support file round trip");
	check(
		std::filesystem::file_size(path) == sizeof(Header),
		"zero-record support file contains only its header"
	);
}

void test_zero_stride_rose(const temporary_directory& directory) {
	using Rose = GraphGeneration::support_transient_graph<1, 20>;
	using RoseReader = GraphGeneration::MappedSupportTransientReader<Rose>;
	using RoseWriter = GraphGeneration::MappedSupportTransientWriter<Rose>;
	static_assert(Rose::SUPPORT_BIT_COUNT == 0);
	static_assert(GraphGeneration::SUPPORT_TRANSIENT_RECORD_SIZE<Rose> == 0);

	const auto path = directory.file("rose.gcg");
	{
		RoseWriter writer(path, 1);
		writer.write(0, Rose::rose());
	}
	RoseReader input(path);
	check(input.size() == 1, "zero-stride rose retains its record count");
	check(input[0] == Rose::rose(), "zero-stride rose round trip");
	check(
		input.header().record_size_bytes == 0,
		"one-vertex support record has zero payload bytes"
	);
	check(
		std::filesystem::file_size(path) == sizeof(Header),
		"one-vertex rose file contains only its header"
	);
}

void test_malformed_headers(const temporary_directory& directory) {
	malformed_header_case(directory, "legacy-version", [](Header& header) {
		header.format_version = 1;
	});
	malformed_header_case(directory, "future-version", [](Header& header) {
		header.format_version = GraphGeneration::MAPPED_GRAPH_FORMAT_VERSION + 1;
	});
	malformed_header_case(directory, "wrong-stride", [](Header& header) {
		++header.record_size_bytes;
	});
	malformed_header_case(directory, "wrong-kind", [](Header& header) {
		header.payload_kind = static_cast<std::uint64_t>(
			GraphGeneration::MappedGraphPayloadKind::transient_graph
		);
	});
	malformed_header_case(directory, "inactive-record", [](Header& header) {
		header.splittable_graph_count = 0;
	});
	malformed_header_case(directory, "unknown-flags", [](Header& header) {
		header.flags = 1;
	});
	malformed_header_case(directory, "wrong-vertices", [](Header& header) {
		++header.vertices;
	});
	malformed_header_case(directory, "wrong-edges", [](Header& header) {
		++header.edges;
	});
}

void test_malformed_file_sizes(const temporary_directory& directory) {
	const auto incomplete = directory.file("incomplete-header.gcg");
	write_records(incomplete, {make_support({{0, 1}})});
	std::filesystem::resize_file(incomplete, sizeof(Header) - 1);
	check_throws(
		[&] { Reader input(incomplete); }, "reject incomplete header"
	);

	const auto short_file = directory.file("short-payload.gcg");
	write_records(short_file, {make_support({{0, 1}})});
	std::filesystem::resize_file(
		short_file, sizeof(Header) + packed_record_size - 1
	);
	check_throws(
		[&] { Reader input(short_file); }, "reject short payload"
	);

	const auto long_file = directory.file("long-payload.gcg");
	write_records(long_file, {make_support({{0, 1}})});
	std::filesystem::resize_file(
		long_file, sizeof(Header) + packed_record_size + 1
	);
	check_throws(
		[&] { Reader input(long_file); }, "reject trailing payload bytes"
	);
}

void test_malformed_records(const temporary_directory& directory) {
	const auto padding = directory.file("nonzero-padding.gcg");
	write_records(padding, {make_support({{0, 1}})});
	// With 66 support bits, only bits zero and one of the ninth byte are used.
	overwrite_payload_byte(padding, packed_record_size - 1, std::byte{0x80});
	check_throws([&] {
		Reader input(padding);
		(void)input[0];
	}, "reject nonzero record padding bits");

	const auto too_many_edges = directory.file("too-many-support-edges.gcg");
	write_records(too_many_edges, {make_support({})});
	// The first 21 packed bits now describe 21 support edges although E=20.
	for (std::size_t byte = 0; byte < 2; ++byte) {
		overwrite_payload_byte(too_many_edges, byte, std::byte{0xff});
	}
	overwrite_payload_byte(too_many_edges, 2, std::byte{0x1f});
	check_throws([&] {
		Reader input(too_many_edges);
		(void)input[0];
	}, "reject a support record with more than E edges");
}

} // namespace

int main() {
	try {
		temporary_directory directory;
		test_round_trip_and_packing(directory);
		test_empty_file(directory);
		test_zero_stride_rose(directory);
		test_malformed_headers(directory);
		test_malformed_file_sizes(directory);
		test_malformed_records(directory);
	} catch (const std::exception& error) {
		std::cerr << "unexpected exception: " << error.what() << '\n';
		return 1;
	}

	if (failures != 0) {
		std::cerr << failures << " mapped support transient test(s) failed\n";
		return 1;
	}
	std::cout << "all mapped support transient tests passed\n";
	return 0;
}
