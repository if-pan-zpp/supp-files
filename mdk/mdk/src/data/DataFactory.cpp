#include "data/DataFactory.hpp"
#include "data/Chains.hpp"
#include "data/Charges.hpp"
#include "data/Masses.hpp"
#include "data/Types.hpp"
using namespace mdk;

template<>
Chains DataFactory::create<Chains>() const {
    return Chains(*model);
}

template<>
Charges DataFactory::create<Charges>() const {
    return Charges(*model, *params);
}

template<>
Masses DataFactory::create<Masses>() const {
    return Masses(*model);
}

template<>
Types DataFactory::create<Types>() const {
    return Types(*model);
}

template<>
Model DataFactory::create<Model>() const {
    return *model;
}

template<>
param::Parameters DataFactory::create<param::Parameters>() const {
    return *params;
}