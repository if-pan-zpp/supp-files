#pragma once
#include "../model/Model.hpp"
#include "../files/param/Parameters.hpp"
#include <memory>
#include <typeindex>
#include <unordered_map>

namespace mdk {
    /**
     * A generic provider of data classes derived from a model and a parameter
     * file. This class is provided chiefly for convenience.
     */
    class DataFactory {
    public:
        DataFactory(Model const& model, param::Parameters const& params):
            model(&model), params(&params) {};

        /** Extract a data class of type Data, if such is defined; otherwise,
         * one will most likely get a linker error, as relevant template
         * specialzations of /p create are placed in a .cpp file.
         * @tparam Data Type of data class to retrieve
         * @return A const reference to a stored data class.
         */
        template<typename Data>
        Data const& data() {
            auto idx = std::type_index(typeid(Data));
            // If a data class has not been emplaced in the savedData,
            // we add it.
            if (savedData.find(idx) == savedData.end()) {
                savedData.emplace(idx, std::make_shared<Data>(create<Data>()));
            }
            return *(Data const*)savedData.at(idx).get();
        }

    private:
        Model const *model;
        param::Parameters const *params;

        /**
         * A realization of a heterogeneous container. It holds values based on
         * the types of elements, i.e. at most one value of each type can be
         * stored. Using std::shared_ptr<void> allows for safe destruction of
         * stored entities.
         */
        std::unordered_map<std::type_index, std::shared_ptr<void>> savedData;

        /**
         * Underlying data class constructor.
         * @tparam Data Data class to construct
         * @return Constructed data class.
         */
        template<typename Data>
        Data create() const;
    };
}