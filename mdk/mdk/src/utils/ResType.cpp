#include "utils/ResType.hpp"
#include "utils/AminoAcid.hpp"
using namespace mdk;
using namespace std;

ResType::operator ResTypeIdx const&() const {
    return code;
}

ResType::ResType(const std::string &name) {
    (int8_t&)code = (int8_t)AminoAcid(name);
}

ResType::operator int8_t() const {
    return (int8_t)code;
}

ResType::operator std::string() const {
    return (string)AminoAcid((int8_t)code);
}

double ResType::mass() const {
    return AminoAcid((int8_t)code).info().mass;
}
