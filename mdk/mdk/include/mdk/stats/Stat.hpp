#pragma once
#include <cstdint>

namespace mdk {
    /**
     * An object containing the coordination stats for a given
     * residue. As opposed to the Fortran version, where the numbers
     * represent the number of contacts formed, here it's the number
     * of slots remaining. This simplifies the checking procedure, as
     * we can simply compare the coordination numbers with zero. We also
     * provide facilities for manipulating these stats, such as addition etc.
     * This object can also naturally represent a stat difference.
     */
    struct Stat {
        /**
         * Number of remaining backbone contacts that can be formed.
         */
        int8_t backbone = 0;

        /**
         * Number of remaining sidechain contacts that can be formed.
         */
        int8_t sidechain = 0;

        /**
         * Number of remaining sidechain contacts that can be formed with
         * hydrophobic residues.
         */
        int8_t hydrophobicSS = 0;

        /**
         * Number of remaining sidechain contacts that can be formed with
         * polar residues.
         */
        int8_t polarSS = 0;

        /**
         * Checks whether a stat is valid; in particular, one can check whether
         * the formation of a contact will exceed the number of available slots
         * this way.
         * @return Whether a stat is valid.
         */
        inline bool valid() const {
            return backbone >= 0 && sidechain >= 0 &&
                   hydrophobicSS >= 0 && polarSS >= 0;
        }

        inline Stat& operator+=(Stat const& other) {
            backbone += other.backbone;
            sidechain += other.sidechain;
            hydrophobicSS += other.hydrophobicSS;
            polarSS += other.polarSS;
            return *this;
        }

        inline Stat operator+(Stat const& other) const {
            auto res = *this;
            res += other;
            return res;
        }

        inline Stat operator-() const {
            auto res = *this;
            res.backbone = -res.backbone;
            res.sidechain = -res.sidechain;
            res.polarSS = -res.polarSS;
            res.hydrophobicSS = -res.hydrophobicSS;
            return res;
        }

        inline Stat& operator-=(Stat const& other) {
            *this += -other;
            return *this;
        }

        inline Stat operator-(Stat const& other) const {
            auto res = *this;
            res -= other;
            return res;
        }
    };
}