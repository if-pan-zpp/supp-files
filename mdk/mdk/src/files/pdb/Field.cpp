#include "utils/Text.hpp"
#include "Field.hpp"
#include <stdexcept>
#include <cstring>
using namespace mdk::pdb;
using namespace std;

void Integer::read(const string &s) {
    *v = stoi((string)view(s, i, j)) + offset;
}

/**
 * This is an analogue of \p snprintf, i.e. it takes a format string and
 * outputs the result into a \p std::string_view.
 * @tparam Args Format argument types.
 * @param sv Destination \p std::string_view.
 * @param fmt Format string.
 * @param args Format arguments.
 */
template<typename... Args>
void format(string_view sv, const char *fmt, Args const&... args) {
    string buf(81, ' ');
    int n = snprintf(buf.data(), sv.size()+1, fmt, args...);
    memcpy((char*)sv.data(), buf.data(), n);
}

void Integer::write(string &s) const {
    auto sv = view(s, i, j);
    format(sv, "%*d", sv.size(), *v - offset);
}

void Real::read(const string &s) {
    *v = stod((string)view(s, i, j)) * scalar;
}

void Real::write(string &s) const {
    auto sv = view(s, i, j);
    format(sv, "%*.*f", n, m, *v / scalar);
}

void String::read(const string &s) {
    if (mode & Exact) *v = view(s, i, j);
    else *v = trim(view(s, i, j));
}

void String::write(string &s) const {
    auto sv = view(s, i, j);
    if (mode & Trim) *v = trim(*v);

    /** Here we determine how much to offset the string when using
     * \p Left or \p Right mode.
     */
    int size = min(sv.size(), v->size());
    int offset = (mode & Left) ? 0 : max<int>(0, sv.size() - v->size());
    memcpy((char*)sv.data() + offset, v->data(), size);
}

void Char::read(const string &s) {
    *v = view(s, i);
}

void Char::write(string &s) const {
    s[i-1] = *v;
}

void SymOp::read(const string &s) {
    *v = trim(view(s, i, j));
}

void SymOp::write(string &s) const {
    auto sv = view(s, i, j);
    format(sv, "%*s", (int)v->size(), v->data());
}

void Literal::read(const string &s) {
    if (view(s, i, j) != lit)
        throw runtime_error("literal not satisfied");
}

void Literal::write(string &s) const {
    auto sv = view(s, i, j);
    format(sv, "%*s", (int)lit.size(), lit.data());
}
