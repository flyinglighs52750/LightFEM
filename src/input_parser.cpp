#include "lightfem/input_parser.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <unordered_set>
#include <utility>
#include <vector>

namespace lightfem {

namespace {

enum class Section {
    None,
    Settings,
    Nodes,
    Materials,
    Elements,
    Boundary,
    Loads,
    BodyForces,
    EdgeLoads,
};

std::string trim(std::string value) {
    const auto first = std::find_if_not(value.begin(), value.end(),
                                        [](const unsigned char c) { return std::isspace(c) != 0; });
    if (first == value.end()) {
        return {};
    }
    const auto last = std::find_if_not(value.rbegin(), value.rend(),
                                       [](const unsigned char c) { return std::isspace(c) != 0; })
                          .base();
    return std::string(first, last);
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

std::vector<std::string> tokenize(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> tokens;
    for (std::string token; stream >> token;) {
        tokens.push_back(std::move(token));
    }
    return tokens;
}

[[noreturn]] void fail(const std::string& source, const std::size_t line, const std::string& message) {
    throw ParseError(source, line, message);
}

int parse_int(const std::string& token, const std::string& source, const std::size_t line,
              const char* field) {
    try {
        std::size_t consumed = 0;
        const long long value = std::stoll(token, &consumed);
        if (consumed != token.size() || value < std::numeric_limits<int>::min() ||
            value > std::numeric_limits<int>::max()) {
            fail(source, line, std::string("invalid integer for ") + field + ": " + token);
        }
        return static_cast<int>(value);
    } catch (const std::invalid_argument&) {
        fail(source, line, std::string("invalid integer for ") + field + ": " + token);
    } catch (const std::out_of_range&) {
        fail(source, line, std::string("integer out of range for ") + field + ": " + token);
    }
}

std::size_t parse_size(const std::string& token, const std::string& source, const std::size_t line,
                       const char* field) {
    try {
        std::size_t consumed = 0;
        const unsigned long long value = std::stoull(token, &consumed);
        if (consumed != token.size() || value > std::numeric_limits<std::size_t>::max()) {
            fail(source, line, std::string("invalid non-negative integer for ") + field + ": " + token);
        }
        return static_cast<std::size_t>(value);
    } catch (const std::invalid_argument&) {
        fail(source, line, std::string("invalid non-negative integer for ") + field + ": " + token);
    } catch (const std::out_of_range&) {
        fail(source, line, std::string("integer out of range for ") + field + ": " + token);
    }
}

double parse_double(const std::string& token, const std::string& source, const std::size_t line,
                    const char* field) {
    try {
        std::size_t consumed = 0;
        const double value = std::stod(token, &consumed);
        if (consumed != token.size() || !std::isfinite(value)) {
            fail(source, line, std::string("invalid finite number for ") + field + ": " + token);
        }
        return value;
    } catch (const std::invalid_argument&) {
        fail(source, line, std::string("invalid number for ") + field + ": " + token);
    } catch (const std::out_of_range&) {
        fail(source, line, std::string("number out of range for ") + field + ": " + token);
    }
}

std::optional<double> parse_constraint(const std::string& token, const std::string& source,
                                       const std::size_t line, const char* field) {
    const std::string normalized = lower(token);
    if (normalized == "free" || normalized == "none" || normalized == "-") {
        return std::nullopt;
    }
    return parse_double(token, source, line, field);
}

Section parse_section(const std::string& name, const std::string& source, const std::size_t line) {
    const std::string section = lower(trim(name));
    if (section == "settings") return Section::Settings;
    if (section == "nodes") return Section::Nodes;
    if (section == "materials") return Section::Materials;
    if (section == "elements") return Section::Elements;
    if (section == "boundary" || section == "boundaries" || section == "dirichlet") {
        return Section::Boundary;
    }
    if (section == "loads" || section == "forces") return Section::Loads;
    if (section == "body_forces" || section == "body_loads") return Section::BodyForces;
    if (section == "edge_loads" || section == "surface_loads" || section == "tractions" ||
        section == "pressures") return Section::EdgeLoads;
    fail(source, line, "unknown section [" + section + "]");
}

void require_count(const std::vector<std::string>& tokens, const std::size_t expected,
                   const std::string& source, const std::size_t line, const char* format) {
    if (tokens.size() != expected) {
        fail(source, line, std::string("expected ") + std::to_string(expected) +
                               " fields (" + format + "), got " + std::to_string(tokens.size()));
    }
}

}  // namespace

ParseError::ParseError(std::string source, const std::size_t line, std::string message)
    : std::runtime_error(source + ":" + std::to_string(line) + ": " + message),
      source_(std::move(source)),
      line_(line) {}

Model InputParser::parse_file(const std::string& path) const {
    std::ifstream input(path);
    if (!input) {
        throw ParseError(path, 0, "cannot open input file");
    }
    return parse(input, path);
}

Model InputParser::parse(std::istream& input, std::string source_name) const {
    Model model;
    model.source_name = source_name;
    Section section = Section::None;
    std::unordered_set<std::string> seen_settings;

    std::string raw;
    std::size_t line_number = 0;
    while (std::getline(input, raw)) {
        ++line_number;
        const auto comment = raw.find('#');
        if (comment != std::string::npos) {
            raw.erase(comment);
        }
        const std::string line = trim(raw);
        if (line.empty()) {
            continue;
        }

        if (line.front() == '[') {
            if (line.back() != ']') {
                fail(source_name, line_number, "section header must end with ']'");
            }
            section = parse_section(line.substr(1, line.size() - 2), source_name, line_number);
            continue;
        }
        if (section == Section::None) {
            fail(source_name, line_number, "data appears before the first section header");
        }

        const auto tokens = tokenize(line);
        switch (section) {
            case Section::Settings: {
                require_count(tokens, 2, source_name, line_number, "key value");
                const std::string key = lower(tokens[0]);
                if (!seen_settings.insert(key).second) {
                    fail(source_name, line_number, "duplicate setting: " + key);
                }
                if (key == "plane_mode") {
                    const std::string mode = lower(tokens[1]);
                    if (mode == "plane_stress") {
                        model.settings.plane_mode = PlaneMode::PlaneStress;
                    } else if (mode == "plane_strain") {
                        model.settings.plane_mode = PlaneMode::PlaneStrain;
                    } else {
                        fail(source_name, line_number, "plane_mode must be plane_stress or plane_strain");
                    }
                } else if (key == "solver") {
                    model.settings.solver = lower(tokens[1]);
                } else if (key == "tolerance") {
                    model.settings.tolerance = parse_double(tokens[1], source_name, line_number, "tolerance");
                } else if (key == "max_iterations") {
                    model.settings.max_iterations =
                        parse_size(tokens[1], source_name, line_number, "max_iterations");
                } else if (key == "reordering") {
                    model.settings.reordering = lower(tokens[1]);
                } else if (key == "assembly") {
                    model.settings.assembly = lower(tokens[1]);
                } else if (key == "threads") {
                    model.settings.threads = parse_size(tokens[1], source_name, line_number, "threads");
                } else if (key == "output") {
                    model.settings.output_file = tokens[1];
                } else if (key == "deformation_scale") {
                    model.settings.deformation_scale =
                        parse_double(tokens[1], source_name, line_number, "deformation_scale");
                } else if (key == "csv_prefix") {
                    model.settings.csv_prefix = tokens[1];
                } else {
                    fail(source_name, line_number, "unknown setting: " + key);
                }
                break;
            }
            case Section::Nodes:
                require_count(tokens, 3, source_name, line_number, "id x y");
                model.nodes.push_back(Node{parse_int(tokens[0], source_name, line_number, "node id"),
                                           Vec2{parse_double(tokens[1], source_name, line_number, "x"),
                                                parse_double(tokens[2], source_name, line_number, "y")},
                                           std::numeric_limits<std::size_t>::max(),
                                           std::numeric_limits<std::size_t>::max(), line_number});
                break;
            case Section::Materials:
                require_count(tokens, 4, source_name, line_number, "id E nu thickness");
                model.materials.push_back(
                    Material{parse_int(tokens[0], source_name, line_number, "material id"),
                             parse_double(tokens[1], source_name, line_number, "E"),
                             parse_double(tokens[2], source_name, line_number, "nu"),
                             parse_double(tokens[3], source_name, line_number, "thickness"), line_number});
                break;
            case Section::Elements:
                require_count(tokens, 6, source_name, line_number, "id n1 n2 n3 n4 material_id");
                model.elements.push_back(Q4Element{
                    parse_int(tokens[0], source_name, line_number, "element id"),
                    {parse_int(tokens[1], source_name, line_number, "n1"),
                     parse_int(tokens[2], source_name, line_number, "n2"),
                     parse_int(tokens[3], source_name, line_number, "n3"),
                     parse_int(tokens[4], source_name, line_number, "n4")},
                    parse_int(tokens[5], source_name, line_number, "material id"), line_number});
                break;
            case Section::Boundary:
                if (tokens.size() == 2) {
                    const std::string type = lower(tokens[1]);
                    const int node_id = parse_int(tokens[0], source_name, line_number, "boundary node id");
                    if (type == "fixed" || type == "pinned") {
                        model.boundaries.push_back(DirichletBC{node_id, 0.0, 0.0, line_number});
                    } else if (type == "roller_x" || type == "symmetry_x") {
                        model.boundaries.push_back(DirichletBC{node_id, 0.0, std::nullopt, line_number});
                    } else if (type == "roller_y" || type == "symmetry_y") {
                        model.boundaries.push_back(DirichletBC{node_id, std::nullopt, 0.0, line_number});
                    } else {
                        fail(source_name, line_number,
                             "boundary type must be fixed, pinned, roller_x, roller_y, symmetry_x, or symmetry_y");
                    }
                } else {
                    require_count(tokens, 3, source_name, line_number, "node_id ux uy or node_id boundary_type");
                    model.boundaries.push_back(
                        DirichletBC{parse_int(tokens[0], source_name, line_number, "boundary node id"),
                                    parse_constraint(tokens[1], source_name, line_number, "ux"),
                                    parse_constraint(tokens[2], source_name, line_number, "uy"), line_number});
                }
                break;
            case Section::Loads:
                require_count(tokens, 3, source_name, line_number, "node_id fx fy");
                model.loads.push_back(NodalLoad{
                    parse_int(tokens[0], source_name, line_number, "load node id"),
                    parse_double(tokens[1], source_name, line_number, "fx"),
                    parse_double(tokens[2], source_name, line_number, "fy"), line_number});
                break;
            case Section::BodyForces:
                require_count(tokens, 3, source_name, line_number, "element_id bx by");
                model.body_forces.push_back(BodyForce{
                    parse_int(tokens[0], source_name, line_number, "element id"),
                    parse_double(tokens[1], source_name, line_number, "bx"),
                    parse_double(tokens[2], source_name, line_number, "by"), line_number});
                break;
            case Section::EdgeLoads: {
                if (tokens.size() < 4) {
                    fail(source_name, line_number, "expected edge load: element_id edge traction tx ty or element_id edge pressure p");
                }
                const int element_id = parse_int(tokens[0], source_name, line_number, "element id");
                const int edge = parse_int(tokens[1], source_name, line_number, "edge");
                const std::string type = lower(tokens[2]);
                if (type == "traction") {
                    require_count(tokens, 5, source_name, line_number, "element_id edge traction tx ty");
                    model.edge_loads.push_back(EdgeLoad{element_id, edge, EdgeLoadType::Traction,
                                                        parse_double(tokens[3], source_name, line_number, "tx"),
                                                        parse_double(tokens[4], source_name, line_number, "ty"), line_number});
                } else if (type == "pressure") {
                    require_count(tokens, 4, source_name, line_number, "element_id edge pressure p");
                    model.edge_loads.push_back(EdgeLoad{element_id, edge, EdgeLoadType::Pressure,
                                                        parse_double(tokens[3], source_name, line_number, "pressure"),
                                                        0.0, line_number});
                } else {
                    fail(source_name, line_number, "edge load type must be traction or pressure");
                }
                break;
            }
            case Section::None:
                break;
        }
    }

    return model;
}

}  // namespace lightfem
