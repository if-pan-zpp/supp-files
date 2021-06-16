#pragma once
#include <vector>
#include <cstdint>
#include <string>

namespace mdk {
    /**
     * An underlying contact type index.
     */
    enum class ContactTypeIdx: int8_t {
        NAT, NAT_BB, NAT_BS, NAT_SB, NAT_SS, SSBOND
    };

    /**
     * An object representing a contact type, along with a number of
     * facilities, conversions from/to other types etc.
     */
    class ContactType {
    public:
        ContactType() = default;

        constexpr explicit ContactType(ContactTypeIdx code): code {code} {};
        explicit operator ContactTypeIdx const&() const;

        constexpr ContactType(int8_t x): ContactType((ContactTypeIdx)x) {};
        operator int8_t() const;

        explicit ContactType(std::string const& name);
        explicit operator std::string() const;

    private:
        ContactTypeIdx code;
        friend struct std::hash<ContactType>;
    };
}

namespace std {
    /**
     * A \p std::hash instantiation for \p ContactType, allows us to use
     * \p ContactType as keys in STL maps.
     */
    template<>
    struct hash<mdk::ContactType> {
        size_t operator()(mdk::ContactType const &resType) const {
            return std::hash<mdk::ContactTypeIdx>()(resType.code);
        }
    };
}