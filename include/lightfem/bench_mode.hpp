#pragma once

#include <cstdlib>

namespace lightfem {

// Measurement switch for the end-to-end scaling study. Each variable selects the
// pre-optimization implementation of one stage, so a single binary can run both
// sides of an A/B comparison back to back. Interleaving the two runs on the same
// machine removes the build-to-build and thermal drift that otherwise swamps the
// difference at small mesh sizes. Normal runs leave every variable unset, and
// the optimized path is the default.
inline bool legacy_stage(const char* name) noexcept {
    const char* value = std::getenv(name);
    return value != nullptr && value[0] != '\0' && value[0] != '0';
}

}  // namespace lightfem
