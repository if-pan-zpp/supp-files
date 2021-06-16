#include "Text.hpp"
#include <algorithm>
using namespace mdk;
using namespace std;

string_view mdk::view(string const& s, int i, int j) {
    return string_view(s.c_str() + (i-1), j-i+1);
}

char mdk::view(string const& s, int i) {
    return s[i - 1];
}

string_view mdk::ltrim(string_view s) {
    auto beg = find_if(s.begin(), s.end(), [](unsigned char c) -> auto {
        return !std::isspace(c);
    });
    auto end = s.end();

    return string_view(beg, end - beg);
}

string_view mdk::rtrim(string_view s) {
    auto beg = s.begin();
    auto end = find_if(s.rbegin(), s.rend(), [](unsigned char c) -> auto {
        return !std::isspace(c);
    }).base();

    return string_view(beg, end - beg);
}

string_view mdk::trim(string_view s) {
    return rtrim(ltrim(s));
}

string mdk::line(istream& is) {
    string s;
    getline(is, s);
    return s;
}

void mdk::skipLine(istream& is) {
    string line;
    getline(is, line);
}

stringstream mdk::lineStream(istream& is) {
    string line;
    getline(is, line);
    stringstream ss;
    ss << line;
    return ss;
}
