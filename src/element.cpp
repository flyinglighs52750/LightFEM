#include "lightfem/element.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lightfem {

namespace {

struct Kinematics {
    std::array<double, 24> b{};  // 3 x 8 row-major strain-displacement matrix
    double det_jacobian{0.0};
};

std::array<double, 9> constitutive_matrix(const Material& material, const PlaneMode plane_mode) {
    std::array<double, 9> d{};
    const double e = material.young_modulus;
    const double nu = material.poisson_ratio;
    if (plane_mode == PlaneMode::PlaneStress) {
        const double factor = e / (1.0 - nu * nu);
        d = {{factor, factor * nu, 0.0,
              factor * nu, factor, 0.0,
              0.0, 0.0, factor * (1.0 - nu) * 0.5}};
    } else {
        const double factor = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
        d = {{factor * (1.0 - nu), factor * nu, 0.0,
              factor * nu, factor * (1.0 - nu), 0.0,
              0.0, 0.0, factor * (1.0 - 2.0 * nu) * 0.5}};
    }
    return d;
}

Kinematics q4_kinematics(const std::array<Vec2, 4>& coordinates, const double xi, const double eta) {
    const std::array<double, 4> dxi{{-0.25 * (1.0 - eta), +0.25 * (1.0 - eta),
                                     +0.25 * (1.0 + eta), -0.25 * (1.0 + eta)}};
    const std::array<double, 4> deta{{-0.25 * (1.0 - xi), -0.25 * (1.0 + xi),
                                      +0.25 * (1.0 + xi), +0.25 * (1.0 - xi)}};

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        dx_dxi += dxi[i] * coordinates[i].x;
        dx_deta += deta[i] * coordinates[i].x;
        dy_dxi += dxi[i] * coordinates[i].y;
        dy_deta += deta[i] * coordinates[i].y;
    }

    const double determinant = dx_dxi * dy_deta - dx_deta * dy_dxi;
    if (!(determinant > 1.0e-14) || !std::isfinite(determinant)) {
        throw std::runtime_error("Q4 stiffness evaluation received a non-positive Jacobian");
    }
    const double inv00 = dy_deta / determinant;
    const double inv01 = -dx_deta / determinant;
    const double inv10 = -dy_dxi / determinant;
    const double inv11 = dx_dxi / determinant;

    Kinematics result;
    result.det_jacobian = determinant;
    for (std::size_t i = 0; i < 4; ++i) {
        const double dnx = inv00 * dxi[i] + inv10 * deta[i];
        const double dny = inv01 * dxi[i] + inv11 * deta[i];
        result.b[2 * i] = dnx;
        result.b[2 * i + 1] = 0.0;
        result.b[8 + 2 * i] = 0.0;
        result.b[8 + 2 * i + 1] = dny;
        result.b[16 + 2 * i] = dny;
        result.b[16 + 2 * i + 1] = dnx;
    }
    return result;
}

std::array<double, 3> multiply_3x8(const std::array<double, 24>& b,
                                    const std::array<double, 8>& displacement) {
    std::array<double, 3> result{};
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 8; ++column) {
            result[row] += b[8 * row + column] * displacement[column];
        }
    }
    return result;
}

std::array<double, 3> multiply_3x3(const std::array<double, 9>& matrix,
                                    const std::array<double, 3>& vector) {
    std::array<double, 3> result{};
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            result[row] += matrix[3 * row + column] * vector[column];
        }
    }
    return result;
}

std::array<double, 4> q4_shape_functions(const double xi, const double eta) {
    return {{0.25 * (1.0 - xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 - eta),
             0.25 * (1.0 + xi) * (1.0 + eta), 0.25 * (1.0 - xi) * (1.0 + eta)}};
}

std::pair<std::size_t, std::size_t> edge_nodes(const int edge) {
    switch (edge) {
        case 1: return {0, 1};
        case 2: return {1, 2};
        case 3: return {2, 3};
        case 4: return {3, 0};
        default: throw std::invalid_argument("Q4 edge must be in [1, 4]");
    }
}

}  // namespace

std::array<std::size_t, 8> q4_dof_indices(const Model& model, const Q4Element& element) {
    if (!model.is_indexed()) {
        throw std::logic_error("model must be indexed before Q4 DOFs are requested");
    }
    std::array<std::size_t, 8> dofs{};
    for (std::size_t i = 0; i < 4; ++i) {
        const Node* node = model.find_node(element.node_ids[i]);
        if (node == nullptr) {
            throw std::logic_error("validated Q4 element references an unknown node");
        }
        dofs[2 * i] = node->dof_x;
        dofs[2 * i + 1] = node->dof_y;
    }
    return dofs;
}

ElementMatrix q4_stiffness(const std::array<Vec2, 4>& coordinates, const Material& material,
                           const PlaneMode plane_mode) {
    constexpr double gauss = 0.57735026918962576451;
    const std::array<std::array<double, 2>, 4> points{{
        {{-gauss, -gauss}}, {{gauss, -gauss}}, {{gauss, gauss}}, {{-gauss, gauss}},
    }};
    const std::array<double, 9> d = constitutive_matrix(material, plane_mode);
    ElementMatrix stiffness{};

    for (const auto& point : points) {
        const Kinematics kinematics = q4_kinematics(coordinates, point[0], point[1]);
        std::array<double, 24> db{};
        for (std::size_t row = 0; row < 3; ++row) {
            for (std::size_t column = 0; column < 8; ++column) {
                for (std::size_t inner = 0; inner < 3; ++inner) {
                    db[8 * row + column] += d[3 * row + inner] * kinematics.b[8 * inner + column];
                }
            }
        }
        const double weight = kinematics.det_jacobian * material.thickness;
        for (std::size_t row = 0; row < 8; ++row) {
            for (std::size_t column = 0; column < 8; ++column) {
                double contribution = 0.0;
                for (std::size_t strain = 0; strain < 3; ++strain) {
                    contribution += kinematics.b[8 * strain + row] * db[8 * strain + column];
                }
                stiffness[8 * row + column] += contribution * weight;
            }
        }
    }
    return stiffness;
}

std::array<double, 8> q4_body_force(const std::array<Vec2, 4>& coordinates, const Material& material,
                                     const double bx, const double by) {
    constexpr double gauss = 0.57735026918962576451;
    const std::array<std::array<double, 2>, 4> points{{
        {{-gauss, -gauss}}, {{gauss, -gauss}}, {{gauss, gauss}}, {{-gauss, gauss}},
    }};
    std::array<double, 8> result{};
    for (const auto& point : points) {
        const Kinematics kinematics = q4_kinematics(coordinates, point[0], point[1]);
        const auto shape = q4_shape_functions(point[0], point[1]);
        const double weight = kinematics.det_jacobian * material.thickness;
        for (std::size_t node = 0; node < 4; ++node) {
            result[2 * node] += shape[node] * bx * weight;
            result[2 * node + 1] += shape[node] * by * weight;
        }
    }
    return result;
}

std::array<double, 8> q4_edge_traction(const std::array<Vec2, 4>& coordinates,
                                        const Material& material, const int edge, const double tx,
                                        const double ty) {
    const auto [first, second] = edge_nodes(edge);
    const double dx = coordinates[second].x - coordinates[first].x;
    const double dy = coordinates[second].y - coordinates[first].y;
    const double length = std::hypot(dx, dy);
    if (!(length > 1.0e-14)) throw std::runtime_error("Q4 edge traction received a zero-length edge");
    const double endpoint_weight = 0.5 * length * material.thickness;
    std::array<double, 8> result{};
    result[2 * first] = endpoint_weight * tx;
    result[2 * first + 1] = endpoint_weight * ty;
    result[2 * second] = endpoint_weight * tx;
    result[2 * second + 1] = endpoint_weight * ty;
    return result;
}

std::array<double, 8> q4_edge_pressure(const std::array<Vec2, 4>& coordinates,
                                        const Material& material, const int edge,
                                        const double pressure) {
    const auto [first, second] = edge_nodes(edge);
    const double dx = coordinates[second].x - coordinates[first].x;
    const double dy = coordinates[second].y - coordinates[first].y;
    const double length = std::hypot(dx, dy);
    if (!(length > 1.0e-14)) throw std::runtime_error("Q4 edge pressure received a zero-length edge");
    // CCW connectivity gives outward normal (dy, -dx) / length; positive p acts inward.
    return q4_edge_traction(coordinates, material, edge, -pressure * dy / length,
                             pressure * dx / length);
}

ElementField q4_center_field(const std::array<Vec2, 4>& coordinates, const Material& material,
                             const PlaneMode plane_mode,
                             const std::array<double, 8>& displacement) {
    const Kinematics kinematics = q4_kinematics(coordinates, 0.0, 0.0);
    ElementField field;
    field.strain = multiply_3x8(kinematics.b, displacement);
    field.stress = multiply_3x3(constitutive_matrix(material, plane_mode), field.strain);
    const double sx = field.stress[0];
    const double sy = field.stress[1];
    const double txy = field.stress[2];
    field.von_mises = std::sqrt(std::max(0.0, sx * sx - sx * sy + sy * sy + 3.0 * txy * txy));
    return field;
}

std::array<GaussPointField, 4> q4_gauss_point_fields(
    const std::array<Vec2, 4>& coordinates, const Material& material, const PlaneMode plane_mode,
    const std::array<double, 8>& displacement) {
    constexpr double gauss = 0.57735026918962576451;
    const std::array<std::array<double, 2>, 4> points{{
        {{-gauss, -gauss}}, {{gauss, -gauss}}, {{gauss, gauss}}, {{-gauss, gauss}},
    }};
    const auto constitutive = constitutive_matrix(material, plane_mode);
    std::array<GaussPointField, 4> result{};
    for (std::size_t i = 0; i < points.size(); ++i) {
        const Kinematics kinematics = q4_kinematics(coordinates, points[i][0], points[i][1]);
        result[i].xi = points[i][0];
        result[i].eta = points[i][1];
        result[i].field.strain = multiply_3x8(kinematics.b, displacement);
        result[i].field.stress = multiply_3x3(constitutive, result[i].field.strain);
        const double sx = result[i].field.stress[0];
        const double sy = result[i].field.stress[1];
        const double txy = result[i].field.stress[2];
        result[i].field.von_mises =
            std::sqrt(std::max(0.0, sx * sx - sx * sy + sy * sy + 3.0 * txy * txy));
    }
    return result;
}

}  // namespace lightfem
