#include <cassert>
#include <iostream>
#include <limits>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include "GraphGeneration/SplitReservation.hpp"

struct Entry {
    std::size_t value = 0;
    bool empty() const noexcept { return value == 0; }
    std::size_t hash() const noexcept { return value; }
    bool operator==(const Entry&) const noexcept = default;
};

int main() {
    using GraphGeneration::step_down_split_reservation;
    constexpr std::size_t buckets = 64;
    constexpr std::size_t bytes_per_slot_row = buckets * sizeof(Entry);
    assert(step_down_split_reservation<Entry>(1000, buckets, 0) == 0);
    assert(step_down_split_reservation<Entry>(1000, buckets, 16 * bytes_per_slot_row - 1) == 0);
    // 16 slots hold 14 graphs per bucket: 896 graphs across 64 buckets.
    assert(step_down_split_reservation<Entry>(896, buckets, 16 * bytes_per_slot_row) == 896);
    assert(step_down_split_reservation<Entry>(897, buckets, 16 * bytes_per_slot_row) == 448);
    assert(step_down_split_reservation<Entry>(0, buckets, 100000) == 0);
    for (auto budget : {16 * bytes_per_slot_row, 31 * bytes_per_slot_row, 32 * bytes_per_slot_row}) {
        for (auto requested : {std::size_t{1}, std::size_t{1000}, std::numeric_limits<std::size_t>::max()}) {
            const auto n = step_down_split_reservation<Entry>(requested, buckets, budget);
            linear_probe_set<Entry> table;
            assert(table.try_reserve(n / buckets + (n % buckets != 0)));
            assert(table.capacity() * bytes_per_slot_row <= budget);
        }
    }
    linear_probe_set<Entry> table;
    assert(table.insert(Entry{1}));
    const auto capacity = table.capacity();
    const pid_t child = fork();
    assert(child >= 0);
    if (child == 0) {
        // Force a genuine allocation failure without consuming system memory.
        rlimit limit;
        if (getrlimit(RLIMIT_AS, &limit)) _exit(2);
        limit.rlim_cur = 0;
        if (setrlimit(RLIMIT_AS, &limit)) _exit(3);
        if (table.try_reserve(1000000)) _exit(4);
        if (table.size() != 1 || table.capacity() != capacity || !table.contains(Entry{1})) _exit(5);
        if (!table.try_reserve(1)) _exit(6);
        _exit(0);
    }
    int status = 0;
    assert(waitpid(child, &status, 0) == child);
    assert(WIFEXITED(status) && WEXITSTATUS(status) == 0);
    std::cout << "Reservation sizing, overflow and allocation-failure checks passed\n";
}
