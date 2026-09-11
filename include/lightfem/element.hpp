#pragma once

#include "lightfem/model.hpp"

#include <array>
#include <cstddef>

namespace lightfem {

using ElementMatrix = std::array<double, 64>;

struct ElementField {
    std::array<double, 3> strain{{0.0, 0.0, 0.0}};  // exx, eyy, gamma_xy
    std::array<double, 3> stress{{0.0, 0.0, 0.0}};  // sxx, syy, tau_xy
    double von_mises{0.0};
};

struct GaussPointField {
    double xi{0.0};
    double eta{0.0};
    ElementField field{};
};

// Returns [ux1, uy1, ux2, uy2, ux3, uy3, ux4, uy4]. Model must be indexed.
[[nodiscard]] std::array<std::size_t, 8> q4_dof_indices(const Model& model,
                                                         const Q4Element& element);

// Computes the 8 x 8 stiffness matrix with 2 x 2 Gauss integration.
[[nodiscard]] ElementMatrix q4_stiffness(const std::array<Vec2, 4>& coordinates,
                                         const Material& material, PlaneMode plane_mode);

// Consistent 8-entry nodal loads for constant body-force density and uniform
// edge traction/pressure. Q4 coordinates must be counter-clockwise.
[[nodiscard]] std::array<double, 8> q4_body_force(const std::array<Vec2, 4>& coordinates,
                                                    const Material& material, double bx, double by);
[[nodiscard]] std::array<double, 8> q4_edge_traction(const std::array<Vec2, 4>& coordinates,
                                                       const Material& material, int edge, double tx,
                                                       double ty);
[[nodiscard]] std::array<double, 8> q4_edge_pressure(const std::array<Vec2, 4>& coordinates,
                                                       const Material& material, int edge, double pressure);

// Evaluates strain and stress at the element centre (xi=eta=0).
[[nodiscard]] ElementField q4_center_field(const std::array<Vec2, 4>& coordinates,
                                            const Material& material, PlaneMode plane_mode,
                                            const std::array<double, 8>& displacement);

// Fields at (-g,-g), (g,-g), (g,g), (-g,g), g=1/sqrt(3).
[[nodiscard]] std::array<GaussPointField, 4> q4_gauss_point_fields(
    const std::array<Vec2, 4>& coordinates, const Material& material, PlaneMode plane_mode,
    const std::array<double, 8>& displacement);

}  // namespace lightfem
