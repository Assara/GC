#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace GraphGeneration {

// A sizing hint, not a graph-count bound. Approximate unlabelled simple graphs
// of minimum degree three by weighted degree sequences and endpoint pairings.
// Only the ratio between neighbouring stages is used, anchored to actual size.
inline std::size_t estimate_next_split_stage(int loop, int vertices, std::size_t size) {
    if (!size || vertices >= 2 * (loop - 1)) return 0;
    // Dense, tiny stages have large symmetry corrections. Avoid extrapolating
    // there; normal table growth is cheap for these stages.
    if (size < 1024) return size;
    const auto log_count = [loop](int v) {
        const int edges = loop + v - 1, surplus = 2 * edges - 3 * v;
        std::vector<long double> weights(surplus + 1), coefficients(surplus + 1), previous;
        // Factor out 1/3! per vertex; k is the degree surplus above three.
        weights[0] = coefficients[0] = 1;
        for (int k = 1; k <= surplus; ++k)
            weights[k] = k + 3 < v ? weights[k - 1] / (k + 3) : 0;
        for (int vertex = 0; vertex < v; ++vertex) {
            previous = coefficients;
            std::fill(coefficients.begin(), coefficients.end(), 0);
            for (int total = 0; total <= surplus; ++total)
                for (int k = 0; k <= total; ++k)
                    coefficients[total] += previous[total - k] * weights[k];
        }
        long double second_moment = 0;
        for (int k = 0; k <= surplus; ++k)
            second_moment += (k + 3) * (k + 2) * weights[k] * previous[surplus - k];
        const auto excess = v * second_moment / coefficients[surplus] / (2 * edges);
        // Pairings / (vertex permutations * degree factorials), with a sparse
        // simplicity correction for loops and parallel edges.
        return std::lgamma(2.L * edges + 1) - edges * std::log(2.L)
            - std::lgamma(edges + 1.L) - std::lgamma(v + 1.L) - v * std::log(6.L)
            + std::log(coefficients[surplus]) - excess / 2 - excess * excess / 4;
    };
    const auto predicted = std::ceil(size * std::exp(log_count(vertices + 1) - log_count(vertices)));
    if (!std::isfinite(predicted) || predicted >= std::numeric_limits<std::size_t>::max()) return 0;
    return static_cast<std::size_t>(predicted);
}

// Correct half the previous multiplicative error, without feeding corrected
// estimates back into the model. Tiny-stage fallback values are not evidence.
class SplitStageEstimator {
    std::size_t previous_model_estimate_ = 0;
public:
    std::size_t next(int loop, int vertices, std::size_t size) {
        const auto model = estimate_next_split_stage(loop, vertices, size);
        const auto previous = previous_model_estimate_;
        previous_model_estimate_ = size >= 1024 ? model : 0;
        if (!previous || size < 1024 || !model) return model;
        const auto error = std::clamp(static_cast<long double>(size) / previous, 0.5L, 2.L);
        const auto corrected = std::ceil(model * std::sqrt(error));
        if (!std::isfinite(corrected) || corrected >= std::numeric_limits<std::size_t>::max()) return model;
        return static_cast<std::size_t>(corrected);
    }
};

} // namespace GraphGeneration
