#pragma once

#include "lightfem/model.hpp"

#include <cstddef>
#include <istream>
#include <stdexcept>
#include <string>

namespace lightfem {

class ParseError : public std::runtime_error {
public:
    ParseError(std::string source, std::size_t line, std::string message);

    [[nodiscard]] const std::string& source() const noexcept { return source_; }
    [[nodiscard]] std::size_t line() const noexcept { return line_; }

private:
    std::string source_;
    std::size_t line_{0};
};

class InputParser {
public:
    [[nodiscard]] Model parse_file(const std::string& path) const;
    [[nodiscard]] Model parse(std::istream& input, std::string source_name = "<stream>") const;
};

}  // namespace lightfem

