#include "data/Chains.hpp"
using namespace mdk;

Chains::Chains(const Model &model) {
    this->model = &model;
    chainIdx = Integers(model.n, -1);
    isTerminal = Bytes(model.n, false);

    for (int chIdx = 0; chIdx < (int)model.chains.size(); ++chIdx) {
        auto const& chain = model.chains[chIdx];
        for (int i = chain.start; i < chain.end; ++i) {
            chainIdx[i] = chIdx;
        }
        isTerminal[chain.start] = isTerminal[chain.end - 1] = true;
    }

    pairs = tuples(2);
    triples = tuples(3);
    quads = tuples(4);

    nativeTriples = nativeTuples(3);
    nativeQuads = nativeTuples(4);
}

Bytes Chains::tuples(int k) const {
    Bytes res(model->n, false);
    for (auto const& chain: model->chains) {
        auto start = chain.start, end = chain.end;
        for (int i = start + (k-2); i+1 < end; ++i) {
            res[i] = true;
        }
    }
    return res;
}

Bytes Chains::nativeTuples(int k) const {
    Bytes res(model->n, false);
    for (auto const& chain: model->chains) {
        for (auto const& spIdx: chain.structuredParts) {
            auto const& sp = model->structuredParts[spIdx];
            auto start = chain.start + sp.off;
            auto end = chain.start + sp.off + sp.len;
            for (int i = start + (k-2); i+1 < end; ++i) {
                res[i] = true;
            }
        }
    }
    return res;
}
