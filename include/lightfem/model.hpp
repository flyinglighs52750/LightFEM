#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace lightfem {

enum class PlaneMode {
    PlaneStress,
    PlaneStrain,
};

struct Vec2 {
    double x{0.0};
    double y{0.0};
};

struct Node {
    int id{0};
    Vec2 position{};
    std::size_t dof_x{std::numeric_limits<std::size_t>::max()};
    std::size_t dof_y{std::numeric_limits<std::size_t>::max()};
    std::size_t source_line{0};
};

struct Material {
    int id{0};
    double young_modulus{0.0};
    double poisson_ratio{0.0};
    double thickness{0.0};
    std::size_t source_line{0};
};

struct Q4Element {
    int id{0};
    std::array<int, 4> node_ids{{0, 0, 0, 0}};
    int material_id{0};
    std::size_t source_line{0};
};

// Two-node axial bar. Resists only the member axial force: its end nodes are
// free to rotate, which is the classical truss (pin-jointed) assumption.
struct TrussElement {
    int id{0};
    std::array<int, 2> node_ids{{0, 0}};
    int material_id{0};
    double area{0.0};  // cross-section area (m^2)
    std::size_t source_line{0};
};

struct DirichletBC {
    int node_id{0};
    std::optional<double> ux{};
    std::optional<double> uy{};
    std::size_t source_line{0};
};

struct NodalLoad {
    int node_id{0};
    double fx{0.0};
    double fy{0.0};
    std::size_t source_line{0};
};

// Constant volumetric force density applied over one Q4 element (force / volume).
struct BodyForce {
    int element_id{0};
    double bx{0.0};
    double by{0.0};
    std::size_t source_line{0};
};

enum class EdgeLoadType {
    Traction,
    Pressure,
};

// Q4 local edge numbering: 1=n1-n2, 2=n2-n3, 3=n3-n4, 4=n4-n1.
// Traction components are force / area; positive pressure acts inward.
struct EdgeLoad {
    int element_id{0};
    int edge{0};
    EdgeLoadType type{EdgeLoadType::Traction};
    double first{0.0};  // tx or pressure
    double second{0.0}; // ty (unused for pressure)
    std::size_t source_line{0};
};

struct SolverSettings {
    PlaneMode plane_mode{PlaneMode::PlaneStress};
    std::string solver{"pcg"};
    double tolerance{1.0e-10};
    std::size_t max_iterations{5000};
    std::string reordering{"none"};
    std::string assembly{"serial"};
    std::size_t threads{0};  // 0: OpenMP runtime default; serial ignores this value.
    // Iterative-phase choices. "csr" and "jacobi" reproduce a textbook
    // preconditioned CG; the alternatives change how many bytes a product moves
    // and how much of the matrix one preconditioner application captures.
    std::string storage{"csr"};       // csr | bsr
    std::string preconditioner{"jacobi"};  // jacobi | block
    std::string subdomain_split{"rcb"};    // rcb | chunk
    std::size_t subdomain_nodes{64};       // nodes per subdomain when preconditioner=block
    std::string output_file{"results.vtk"};
    double deformation_scale{1.0};
    std::string csv_prefix{};  // Empty: derive from VTK output file stem.
};

struct ModelSummary {
    std::size_t nodes{0};
    std::size_t elements{0};
    std::size_t materials{0};
    std::size_t boundary_records{0};
    std::size_t loads{0};
    std::size_t body_forces{0};
    std::size_t edge_loads{0};
    std::size_t dofs{0};
    std::size_t trusses{0};
};

class Model {
public:
    std::vector<Node> nodes;
    std::vector<Material> materials;
    std::vector<Q4Element> elements;
    std::vector<TrussElement> trusses;
    std::vector<DirichletBC> boundaries;
    std::vector<NodalLoad> loads;
    std::vector<BodyForce> body_forces;
    std::vector<EdgeLoad> edge_loads;
    SolverSettings settings;
    std::string source_name;

    // Call only after semantic validation succeeds. The method sorts entities
    // by ID, assigns two displacement DOFs per node, and builds lookup tables.
    void sort_and_index();

    // Reorders nodes according to a full ID permutation, then assigns DOFs in
    // that order. Materials and elements retain deterministic ID ordering.
    void index_with_node_order(const std::vector<int>& node_order);

    [[nodiscard]] const Node* find_node(int id) const noexcept;
    [[nodiscard]] const Material* find_material(int id) const noexcept;
    [[nodiscard]] const Q4Element* find_element(int id) const noexcept;
    [[nodiscard]] ModelSummary summary() const noexcept;
    [[nodiscard]] bool is_indexed() const noexcept { return indexed_; }

private:
    std::unordered_map<int, std::size_t> node_index_;
    std::unordered_map<int, std::size_t> material_index_;
    std::unordered_map<int, std::size_t> element_index_;
    bool indexed_{false};

    void rebuild_indices();
};

[[nodiscard]] const char* to_string(PlaneMode mode) noexcept;

}  // namespace lightfem
