#include "lightfem/model.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace lightfem {

namespace {

template <typename T>
void require_unique_ids(const std::vector<T>& entities, const char* label) {
    for (std::size_t i = 1; i < entities.size(); ++i) {
        if (entities[i - 1].id == entities[i].id) {
            throw std::logic_error(std::string("duplicate ") + label + " ID: " +
                                   std::to_string(entities[i].id));
        }
    }
}

}  // namespace

void Model::sort_and_index() {
    std::sort(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) { return a.id < b.id; });
    std::sort(materials.begin(), materials.end(),
              [](const Material& a, const Material& b) { return a.id < b.id; });
    std::sort(elements.begin(), elements.end(),
              [](const Q4Element& a, const Q4Element& b) { return a.id < b.id; });
    std::sort(trusses.begin(), trusses.end(),
              [](const TrussElement& a, const TrussElement& b) { return a.id < b.id; });

    require_unique_ids(nodes, "node");
    require_unique_ids(materials, "material");
    require_unique_ids(elements, "element");
    require_unique_ids(trusses, "truss");

    rebuild_indices();
}

void Model::index_with_node_order(const std::vector<int>& node_order) {
    if (node_order.size() != nodes.size()) {
        throw std::invalid_argument("node reordering must contain every model node exactly once");
    }
    require_unique_ids(nodes, "node");
    std::sort(materials.begin(), materials.end(),
              [](const Material& a, const Material& b) { return a.id < b.id; });
    std::sort(elements.begin(), elements.end(),
              [](const Q4Element& a, const Q4Element& b) { return a.id < b.id; });
    std::sort(trusses.begin(), trusses.end(),
              [](const TrussElement& a, const TrussElement& b) { return a.id < b.id; });
    require_unique_ids(materials, "material");
    require_unique_ids(elements, "element");
    require_unique_ids(trusses, "truss");

    std::unordered_map<int, Node> by_id;
    by_id.reserve(nodes.size());
    for (const Node& node : nodes) by_id.emplace(node.id, node);
    std::vector<Node> reordered;
    reordered.reserve(nodes.size());
    for (const int node_id : node_order) {
        const auto it = by_id.find(node_id);
        if (it == by_id.end()) {
            throw std::invalid_argument("node reordering references an unknown or duplicate node ID");
        }
        reordered.push_back(it->second);
        by_id.erase(it);
    }
    if (!by_id.empty()) throw std::invalid_argument("node reordering omits one or more model nodes");
    nodes = std::move(reordered);
    rebuild_indices();
}

void Model::rebuild_indices() {

    node_index_.clear();
    material_index_.clear();
    element_index_.clear();

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        nodes[i].dof_x = 2 * i;
        nodes[i].dof_y = 2 * i + 1;
        node_index_.emplace(nodes[i].id, i);
    }
    for (std::size_t i = 0; i < materials.size(); ++i) {
        material_index_.emplace(materials[i].id, i);
    }
    for (std::size_t i = 0; i < elements.size(); ++i) {
        element_index_.emplace(elements[i].id, i);
    }
    indexed_ = true;
}

const Node* Model::find_node(const int id) const noexcept {
    if (indexed_) {
        const auto it = node_index_.find(id);
        return it == node_index_.end() ? nullptr : &nodes[it->second];
    }
    const auto it = std::find_if(nodes.begin(), nodes.end(), [id](const Node& node) { return node.id == id; });
    return it == nodes.end() ? nullptr : &*it;
}

const Material* Model::find_material(const int id) const noexcept {
    if (indexed_) {
        const auto it = material_index_.find(id);
        return it == material_index_.end() ? nullptr : &materials[it->second];
    }
    const auto it = std::find_if(materials.begin(), materials.end(),
                                 [id](const Material& material) { return material.id == id; });
    return it == materials.end() ? nullptr : &*it;
}

const Q4Element* Model::find_element(const int id) const noexcept {
    if (indexed_) {
        const auto it = element_index_.find(id);
        return it == element_index_.end() ? nullptr : &elements[it->second];
    }
    const auto it = std::find_if(elements.begin(), elements.end(),
                                 [id](const Q4Element& element) { return element.id == id; });
    return it == elements.end() ? nullptr : &*it;
}

ModelSummary Model::summary() const noexcept {
    return ModelSummary{nodes.size(), elements.size(), materials.size(), boundaries.size(), loads.size(),
                        body_forces.size(), edge_loads.size(), 2 * nodes.size(), trusses.size()};
}

const char* to_string(const PlaneMode mode) noexcept {
    switch (mode) {
        case PlaneMode::PlaneStress:
            return "plane_stress";
        case PlaneMode::PlaneStrain:
            return "plane_strain";
    }
    return "unknown";
}

}  // namespace lightfem
