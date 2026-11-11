#include "lightfem/reordering.hpp"

#include <algorithm>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace lightfem {

std::vector<int> reverse_cuthill_mckee(const Model& model) {
    if (!model.is_indexed()) throw std::logic_error("model must be indexed before RCM ordering");

    std::unordered_map<int, std::unordered_set<int>> adjacency;
    adjacency.reserve(model.nodes.size());
    for (const Node& node : model.nodes) adjacency.emplace(node.id, std::unordered_set<int>{});
    for (const Q4Element& element : model.elements) {
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = i + 1; j < 4; ++j) {
                adjacency[element.node_ids[i]].insert(element.node_ids[j]);
                adjacency[element.node_ids[j]].insert(element.node_ids[i]);
            }
        }
    }
    for (const TrussElement& truss : model.trusses) {
        adjacency[truss.node_ids[0]].insert(truss.node_ids[1]);
        adjacency[truss.node_ids[1]].insert(truss.node_ids[0]);
    }

    const auto by_degree_then_id = [&adjacency](const int left, const int right) {
        const std::size_t left_degree = adjacency.at(left).size();
        const std::size_t right_degree = adjacency.at(right).size();
        return left_degree != right_degree ? left_degree < right_degree : left < right;
    };

    std::vector<int> seeds;
    seeds.reserve(model.nodes.size());
    for (const Node& node : model.nodes) seeds.push_back(node.id);
    std::sort(seeds.begin(), seeds.end(), by_degree_then_id);

    std::unordered_set<int> visited;
    visited.reserve(model.nodes.size());
    std::vector<int> order;
    order.reserve(model.nodes.size());
    for (const int seed : seeds) {
        if (!visited.insert(seed).second) continue;
        std::queue<int> pending;
        pending.push(seed);
        while (!pending.empty()) {
            const int current = pending.front();
            pending.pop();
            order.push_back(current);
            std::vector<int> neighbours(adjacency.at(current).begin(), adjacency.at(current).end());
            std::sort(neighbours.begin(), neighbours.end(), by_degree_then_id);
            for (const int neighbour : neighbours) {
                if (visited.insert(neighbour).second) pending.push(neighbour);
            }
        }
    }
    std::reverse(order.begin(), order.end());
    return order;
}

}  // namespace lightfem
