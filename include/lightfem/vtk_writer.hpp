#pragma once

#include "lightfem/model.hpp"
#include "lightfem/solver.hpp"

#include <string>

namespace lightfem {

// Writes an ASCII legacy VTK unstructured grid using scaled deformed coordinates,
// plus nodal smoothed stresses and centre/Gauss-point element fields.
void write_vtk(const Model& model, const AnalysisResult& result, const std::string& path);

struct CsvExportPaths {
    std::string nodes;
    std::string gauss_points;
};

// Writes <prefix>_nodes.csv and <prefix>_gauss_points.csv.
[[nodiscard]] CsvExportPaths write_csv_results(const Model& model, const AnalysisResult& result,
                                                const std::string& prefix);

}  // namespace lightfem
