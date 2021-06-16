#include "files/cmap/LegacyParser.hpp"
#include "utils/Units.hpp"
#include <sstream>

using namespace mdk::cmap;
using namespace std;

ContactMap LegacyParser::read(std::istream &is) {
    ContactMap cmap;

    int numOfContacts, numOfResidues;
    is >> cmap.offset;
    is >> numOfResidues;
    is >> numOfContacts;

    cmap.contacts = vector<ContactMap::Contact>(numOfContacts);
    for (int i = 0; i < numOfContacts; ++i) {
        auto& contact = cmap.contacts[i];
        is >> contact.res[0] >> contact.res[1] >> contact.dist0;
        --contact.res[0];
        --contact.res[1];
        contact.dist0 *= f77unit * pow(2.0, 1./6.);
    }

    cmap.len = numOfResidues;
    cmap.angle = vector<double>(numOfResidues);
    cmap.dihedral = vector<double>(numOfResidues);

    for (int i = 0; i < numOfResidues; ++i) {
        is >> cmap.angle[i] >> cmap.dihedral[i];
        cmap.angle[i] *= rad;
        cmap.dihedral[i] *= rad;
    }

    return cmap;
}

std::ostream &LegacyParser::write(ostream &os, ContactMap const& cmap) {
    os << cmap.offset << endl;
    os << cmap.len << endl;
    os << cmap.contacts.size() << endl;

    for (auto& contact: cmap.contacts) {
        os << contact.res[0] + 1 << "\t";
        os << contact.res[1] + 1 << "\t";
        os << contact.dist0 / f77unit << endl;
    }

    for (int i = 0; i < cmap.len; ++i) {
        os << cmap.angle[i] / rad << "\t";
        os << cmap.dihedral[i] / rad << endl;
    }

    return os;
}
