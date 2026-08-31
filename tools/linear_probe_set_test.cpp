#include <bit>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <unordered_set>

#include "LinearProbeSet.hpp"
#include "GraphGeneration/SupportTransientGraph.hpp"
#include "graph.hpp"

namespace {

int failures = 0;

void check(bool condition, const char* message) {
	if (!condition) {
		std::cerr << "FAILED: " << message << '\n';
		++failures;
	}
}

struct collision_value {
	int value = 0;

	std::size_t hash() const noexcept { return 3; }
	bool empty() const noexcept { return value == 0; }
	bool operator==(const collision_value&) const noexcept = default;
};

static_assert(
	sizeof(linear_probe_set<collision_value>)
		== sizeof(std::unique_ptr<collision_value[]>) + 2 * sizeof(std::size_t),
	"linear_probe_set must contain only its array pointer, size, and capacity"
);

void test_collisions_duplicates_and_iteration() {
	linear_probe_set<collision_value> set;
	check(set.capacity() == 0 && set.data() == nullptr,
		"empty set has no allocation");
	check(set.insert(collision_value{1}), "first colliding value inserts");
	check(set.insert(collision_value{2}), "second colliding value inserts");
	check(set.insert(collision_value{3}), "wrapped colliding value inserts");
	check(!set.insert(collision_value{1}), "duplicate is rejected");
	check(set.contains(collision_value{1}), "first value found");
	check(set.contains(collision_value{2}), "second value found");
	check(!set.contains(collision_value{7}), "missing collision value not found");

	int sum = 0;
	for (std::size_t index = 0; index < set.capacity(); ++index) {
		if (!set.data()[index].empty()) sum += set.data()[index].value;
	}
	check(sum == 6, "direct array scan skips empty slots");

	check(!set.insert(collision_value{}), "empty sentinel is ignored");
}

void test_growth_and_reserve() {
	linear_probe_set<collision_value> set;
	for (int value = 1; value <= 1000; ++value) {
		check(set.insert(collision_value{value}), "growing insertion succeeds");
		check(std::has_single_bit(set.capacity()), "growth keeps power-of-two capacity");
	}
	check(set.size() == 1000, "growth preserves size");
	for (int value = 1; value <= 1000; ++value) {
		check(set.contains(collision_value{value}),
			"growth preserves elements");
	}
	const std::size_t old_capacity = set.capacity();
	set.reserve(10000);
	check(set.size() == 1000, "reserve preserves element count");
	check(set.capacity() > old_capacity,
		"reserve increases capacity for the requested elements");
	check(std::has_single_bit(set.capacity()),
		"reserved capacity remains a power of two");
	for (int value = 1; value <= 1000; ++value) {
		check(set.contains(collision_value{value}),
			"reserve preserves elements");
	}

	linear_probe_set<collision_value> empty_set;
	empty_set.reserve(100);
	check(empty_set.size() == 0 && empty_set.capacity() >= 100,
		"an empty set can reserve storage");
}

void test_graph_sentinels_and_hash() {
	using MultiVertex = Graph<3, 2, 0, 0, 0, 0, fieldType>;
	MultiVertex empty_graph;
	check(empty_graph.empty(), "multi-vertex graph uses all-zero empty sentinel");

	MultiVertex path;
	path.setEdge(0, 0, 1);
	path.setEdge(1, 1, 2);
	check(!path.empty(), "valid multi-vertex graph is not empty");
	check(path.hash() == std::hash<MultiVertex>{}(path),
		"graph member and std hashes agree");

	linear_probe_set<MultiVertex> graphs;
	check(graphs.insert(path), "graph inserts into linear set");
	check(graphs.contains(path), "graph is found in linear set");

	using OneVertex = Graph<1, 3, 0, 0, 0, 0, fieldType>;
	OneVertex empty_one_vertex;
	check(empty_one_vertex.empty(), "one-vertex graph uses all-one empty sentinel");
	OneVertex rose;
	rose.half_edges.fill(0);
	check(!rose.empty(), "one-vertex rose is not the empty sentinel");
	linear_probe_set<OneVertex> roses;
	check(roses.insert(rose) && roses.contains(rose),
		"one-vertex rose can be stored");
}

void test_transient_graph_sentinel() {
	using Transient = GraphGeneration::support_transient_graph<3, 5, fieldType>;
	Transient empty_transient;
	check(empty_transient.empty(), "empty transient support is an empty slot");
	const Transient triangle = Transient::triangle();
	check(!triangle.empty(), "connected transient support is not empty");
	check(triangle.hash() == triangle.hash_value(),
		"transient member hash matches its established hash");
	linear_probe_set<Transient> transients;
	check(transients.insert(triangle) && transients.contains(triangle),
		"transient graph can be stored");
}

} // namespace

int main() {
	test_collisions_duplicates_and_iteration();
	test_growth_and_reserve();
	test_graph_sentinels_and_hash();
	test_transient_graph_sentinel();
	if (failures != 0) {
		std::cerr << failures << " linear probing hash set tests failed\n";
		return 1;
	}
	std::cout << "all linear probing hash set tests passed\n";
	return 0;
}
