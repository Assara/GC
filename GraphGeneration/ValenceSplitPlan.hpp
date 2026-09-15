#pragma once

#include <algorithm>
#include <map>
#include <vector>
#include "types.hpp"

namespace GraphGeneration {

using ValenceArray = std::vector<Int>;
struct ValenceSplit {
    int vertex;
    Int moved;
    bool operator==(const ValenceSplit&) const = default;
};

// Incidence zero stays on the old vertex. Keep both unequal side sizes:
// which side contains that incidence depends on the particular graph.
inline std::map<ValenceArray, std::vector<ValenceSplit>>
valence_split_plans(const ValenceArray& parent) {
    std::map<ValenceArray, std::vector<ValenceSplit>> plans;
    for (int vertex = int(parent.size()) - 1; vertex >= 0; --vertex) {
        const int degree = parent[vertex];
        for (int moved = 2; moved <= degree - 2; ++moved) {
            auto child = parent;
            child[vertex] = degree - moved + 1;
            child.push_back(moved + 1);
            std::sort(child.begin(), child.end());
            plans[child].push_back({vertex, Int(moved)});
        }
    }
    return plans;
}

} // namespace GraphGeneration
