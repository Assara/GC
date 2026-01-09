#include <chrono>
#include <cstdint>
#include <iostream>

struct timer_accum {
    using clock = std::chrono::steady_clock;

    double seconds = 0.0;
    std::uint64_t calls = 0;

    struct guard {
        timer_accum& t;
        clock::time_point start;
        explicit guard(timer_accum& t_) : t(t_), start(clock::now()) {}
        ~guard() {
            t.seconds += std::chrono::duration<double>(
                clock::now() - start
            ).count();
            ++t.calls;
        }
    };

    void print(const char* name) const {
        std::cout << name << ": " << seconds << " s"
                  << " (" << calls << " calls";
        if (calls) std::cout << ", avg " << (seconds / calls) << " s";
        std::cout << ")\n";
    }
};
