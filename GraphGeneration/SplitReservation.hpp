#pragma once

#include <algorithm>
#include <bit>
#include <charconv>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include "LinearProbeSet.hpp"

namespace GraphGeneration {

// MemAvailable already excludes the resident parents. Leave half the remaining
// headroom for the OS, worker scratch space and later growth. Also honour finite
// cgroup-v2 limits at every ancestor of the current process's group.
inline std::size_t split_reservation_budget() {
    std::ifstream memory("/proc/meminfo");
    std::size_t available = 0;
    for (std::string key, rest; memory >> key;) {
        if (key == "MemAvailable:") { memory >> available; break; }
        std::getline(memory, rest);
    }
    available *= 1024;
    std::ifstream groups("/proc/self/cgroup");
    for (std::string line; std::getline(groups, line);) {
        if (!line.starts_with("0::/")) continue;
        const std::filesystem::path root("/sys/fs/cgroup");
        auto group = root / std::filesystem::path(line.substr(3)).lexically_normal().relative_path();
        for (;;) {
            std::ifstream limit_file(group / "memory.max"), used_file(group / "memory.current");
            std::size_t limit, used;
            if ((limit_file >> limit) && (used_file >> used))
                available = std::min(available, limit > used ? limit - used : 0);
            if (group == root) break;
            group = group.parent_path();
        }
        break;
    }
    auto budget = available / 2;
    // Optional reservation cap, not a process-wide memory limit. Useful when
    // sharing the machine and for exercising the step-down on small runs.
    if (const char* value = std::getenv("GC_TRIANGLE_RESERVE_MIB")) {
        std::string_view input(value);
        std::size_t mib = 0;
        const auto [end, error] = std::from_chars(input.data(), input.data() + input.size(), mib);
        constexpr std::size_t unit = 1024 * 1024;
        if (error != std::errc{} || end != input.data() + input.size()
            || mib > std::numeric_limits<std::size_t>::max() / unit)
            throw std::runtime_error("GC_TRIANGLE_RESERVE_MIB must be a nonnegative integer");
        budget = std::min(budget, mib * unit);
    }
    return budget;
}

// Account for actual power-of-two table allocations, not just graph payload.
// Work from the byte budget first, so even an enormous estimate cannot overflow.
template <linear_probe_value Entry>
std::size_t step_down_split_reservation(std::size_t estimate, std::size_t buckets,
                                      std::size_t budget) {
    if (!buckets) return 0;
    using Set = linear_probe_set<Entry>;
    const auto slots = std::bit_floor(budget / buckets / sizeof(Entry));
    if (slots < Set::MINIMUM_CAPACITY) return 0;
    const auto supported = (slots / Set::LOAD_DENOMINATOR) * Set::LOAD_NUMERATOR * buckets;
    while (estimate > supported) estimate /= 2;
    return estimate;
}

} // namespace GraphGeneration
