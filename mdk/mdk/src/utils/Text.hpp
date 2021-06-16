#pragma once
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace mdk {
    std::string_view view(std::string const &s, int i, int j);
    char view(std::string const& s, int i);

    std::string_view ltrim(std::string_view s);
    std::string_view rtrim(std::string_view s);
    std::string_view trim(std::string_view s);

    std::string line(std::istream& is);
    void skipLine(std::istream& is);
    std::stringstream lineStream(std::istream& is);
}