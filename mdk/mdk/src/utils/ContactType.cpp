#include "utils/ContactType.hpp"
#include <unordered_map>
using namespace mdk;
using namespace std;

ContactType::operator ContactTypeIdx const&() const {
    return code;
}

ContactType::ContactType(const std::string &name) {
    static unordered_map<string, ContactTypeIdx> conversions = {
        { "NAT",    ContactTypeIdx::NAT },
        { "NAT_BB", ContactTypeIdx::NAT_BB },
        { "NAT_BS", ContactTypeIdx::NAT_BS },
        { "NAT_SB", ContactTypeIdx::NAT_SB },
        { "NAT_SS", ContactTypeIdx::NAT_SS },
        { "SSBOND", ContactTypeIdx::SSBOND }
    };
    code = conversions.at(name);
}

ContactType::operator int8_t() const {
    return (int8_t)code;
}

ContactType::operator std::string() const {
    static unordered_map<ContactTypeIdx, string> conversions = {
        {ContactTypeIdx::NAT,    "NAT" },
        {ContactTypeIdx::NAT_BB, "NAT_BB" },
        {ContactTypeIdx::NAT_BS, "NAT_BS" },
        {ContactTypeIdx::NAT_SB, "NAT_SB" },
        {ContactTypeIdx::NAT_SS, "NAT_SS" },
        {ContactTypeIdx::SSBOND, "SSBOND" }
    };
    return conversions.at(code);
}
