#pragma once
#include <vector>
#include <tuple>

/* Include this file in order to provide support for using tuples and pairs
 * as hash map keys. Taken from https://stackoverflow.com/a/55113454/13180615. */

template<class T>
class HashMixer {
public:
    T const &value;
    explicit HashMixer(T const &value): value {value} {};

    size_t operator,(size_t seed) const {
        seed ^= std::hash<T>()(value) + 0x9e3779b9 + (seed << 6u) + (seed >> 2u);
        return seed;
    };
};

namespace std {
    template<typename... Types>
    struct hash<tuple<Types...>> {
        size_t operator()(tuple<Types...> const &t) const {
            return std::apply([](Types const&... tx) -> size_t {
                return (HashMixer(tx), ..., 0);
            }, t);
        }
    };

    template<typename T1, typename T2>
    struct hash<pair<T1, T2>> {
        size_t operator()(pair<T1, T2> const &t) const {
            return hash<tuple<T1, T2>>()(t);
        }
    };
}
