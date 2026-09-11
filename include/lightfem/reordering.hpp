#pragma once

#include "lightfem/model.hpp"

#include <vector>

namespace lightfem {

// Returns a node-ID permutation from the reverse Cuthill-McKee algorithm.
// The model must be indexed and semantically validated.
[[nodiscard]] std::vector<int> reverse_cuthill_mckee(const Model& model);

}  // namespace lightfem
