#pragma once

#include "lightfem/model.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace lightfem {

enum class Severity {
    Error,
    Warning,
};

struct Diagnostic {
    Severity severity{Severity::Error};
    std::string code;
    std::string message;
    std::size_t source_line{0};
};

class ValidationReport {
public:
    void error(std::string code, std::string message, std::size_t source_line = 0);
    void warning(std::string code, std::string message, std::size_t source_line = 0);

    [[nodiscard]] bool ok() const noexcept;
    [[nodiscard]] std::size_t error_count() const noexcept;
    [[nodiscard]] std::size_t warning_count() const noexcept;
    [[nodiscard]] const std::vector<Diagnostic>& diagnostics() const noexcept { return diagnostics_; }

private:
    std::vector<Diagnostic> diagnostics_;
};

class ModelValidator {
public:
    [[nodiscard]] ValidationReport validate(const Model& model) const;
};

[[nodiscard]] const char* to_string(Severity severity) noexcept;

}  // namespace lightfem

