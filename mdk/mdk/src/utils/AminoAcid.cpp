#include "utils/AminoAcid.hpp"
#include "utils/Units.hpp"
#include <algorithm>
using namespace mdk;
using namespace std;

std::vector<AminoAcidIdx> AminoAcid::types() {
    return {
        AminoAcidIdx::ALA, AminoAcidIdx::ARG, AminoAcidIdx::ASN, AminoAcidIdx::ASP, AminoAcidIdx::CYS,
        AminoAcidIdx::GLU, AminoAcidIdx::GLN, AminoAcidIdx::GLY, AminoAcidIdx::HIS, AminoAcidIdx::ILE,
        AminoAcidIdx::LEU, AminoAcidIdx::LYS, AminoAcidIdx::MET, AminoAcidIdx::PHE, AminoAcidIdx::PRO,
        AminoAcidIdx::SER, AminoAcidIdx::THR, AminoAcidIdx::TRP, AminoAcidIdx::TYR, AminoAcidIdx::VAL
    };
}

std::vector<AminoAcid> AminoAcid::aminoAcids() {
    std::vector<AminoAcid> _aminoAcids;
    for (auto code: codes()) {
        _aminoAcids.emplace_back(code);
    }
    return _aminoAcids;
}

AminoAcid::operator AminoAcidIdx const&() const {
    return type;
}

AminoAcid::operator int8_t() const {
    return (int8_t)type;
}

AminoAcid::AminoAcid(char code) {
    static const auto _codes = codes();
    type = (AminoAcidIdx)_codes.find(code);
}

std::string AminoAcid::codes() {
    return "ARNDCEQGHILKMFPSTWYV";
}

AminoAcid::operator char const&() const {
    static const auto _codes = codes();
    return _codes[(int8_t)type];
}

bool AminoAcid::isProper(std::string const& name) {
    static const auto _names = names();
    return find(_names.begin(), _names.end(), name) != _names.end();
}

AminoAcid::AminoAcid(std::string const& name) {
    static const auto _names = names();
    auto iter = find(_names.begin(), _names.end(), name);
    if (iter != _names.end()) {
        static const auto _types = types();
        auto idx = distance(_names.begin(), iter);
        type = _types[idx];
    }
}

std::vector<std::string> AminoAcid::names() {
    return {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    };
}

AminoAcid::operator string() const {
    static const auto _names = names();
    return _names[(int8_t)type];
}

bool AminoAcid::operator==(const AminoAcid &other) const {
    return type == other.type;
}

bool AminoAcid::operator!=(const AminoAcid &other) const {
    return type != other.type;
}

bool AminoAcid::operator<(const AminoAcid &other) const {
    return type < other.type;
}

bool AminoAcid::operator<=(const AminoAcid &other) const {
    return type <= other.type;
}

bool AminoAcid::operator>(const AminoAcid &other) const {
    return type > other.type;
}

bool AminoAcid::operator>=(const AminoAcid &other) const {
    return type >= other.type;
}

using ResidueData = std::unordered_map<AminoAcid, AminoAcidInfo>;
static ResidueData createResData() {
    ResidueData types;

    for (auto const& acid: AminoAcid::aminoAcids()) {
        auto& acidInfo = types[acid];
        acidInfo.atomInfo["N"].radius = 1.64;
        acidInfo.atomInfo["CA"].radius = 1.88;
        acidInfo.atomInfo["C"].radius = 1.61;
        acidInfo.atomInfo["O"].radius = 1.42;
        acidInfo.atomInfo["OXT"].radius = 1.42; // ??

        for (auto& [type, entry]: acidInfo.atomInfo) {
            entry.inBackbone = true;
        }
    }

    auto& PRO = types[(AminoAcid)"PRO"];
    PRO.mass = 97.052763875;
    PRO.atomInfo["CB"].radius = 1.88;
    PRO.atomInfo["CG"].radius = 1.88;
    PRO.atomInfo["CD"].radius = 1.88;

    auto& GLN = types[(AminoAcid)"GLN"];
    GLN.mass = 128.058577540;
    GLN.atomInfo["CB"].radius = 1.88;
    GLN.atomInfo["CB"].radius = 1.88;
    GLN.atomInfo["CG"].radius = 1.88;
    GLN.atomInfo["CD"].radius = 1.61;
    GLN.atomInfo["OE1"].radius = 1.42;
    GLN.atomInfo["NE2"].radius = 1.64;

    auto& CYS = types[(AminoAcid)"CYS"];
    CYS.mass = 103.009184505;
    CYS.atomInfo["CB"].radius = 1.88;
    CYS.atomInfo["SG"].radius = 1.77;

    auto& VAL = types[(AminoAcid)"VAL"];
    VAL.mass = 99.068413945;
    VAL.atomInfo["CB"].radius = 1.88;
    VAL.atomInfo["CG1"].radius = 1.88;
    VAL.atomInfo["CG2"].radius = 1.88;

    auto& PHE = types[(AminoAcid)"PHE"];
    PHE.mass = 147.068413945;
    PHE.atomInfo["CB"].radius = 1.88;
    PHE.atomInfo["CG"].radius = 1.88;
    PHE.atomInfo["CD1"].radius = 1.61;
    PHE.atomInfo["CD2"].radius = 1.76;
    PHE.atomInfo["CE1"].radius = 1.76;
    PHE.atomInfo["CE2"].radius = 1.76;
    PHE.atomInfo["CZ"].radius = 1.76;

    auto& MET = types[(AminoAcid)"MET"];
    MET.mass = 131.040484645;
    MET.atomInfo["CB"].radius = 1.88;
    MET.atomInfo["CG"].radius = 1.88;
    MET.atomInfo["SD"].radius = 1.77;
    MET.atomInfo["CE"].radius = 1.88;

    auto& ILE = types[(AminoAcid)"ILE"];
    ILE.mass = 113.084064015;
    ILE.atomInfo["CB"].radius = 1.88;
    ILE.atomInfo["CG1"].radius = 1.88;
    ILE.atomInfo["CG2"].radius = 1.88;
    ILE.atomInfo["CD1"].radius = 1.88;

    auto& ASP = types[(AminoAcid)"ASP"];
    ASP.mass = 115.026943065;
    ASP.atomInfo["CB"].radius = 1.88;
    ASP.atomInfo["CG"].radius = 1.61;
    ASP.atomInfo["OD1"].radius = 1.46;
    ASP.atomInfo["OD2"].radius = 1.42;

    auto& GLU = types[(AminoAcid)"GLU"];
    GLU.mass = 129.042593135;
    GLU.atomInfo["CB"].radius = 1.88;
    GLU.atomInfo["CG"].radius = 1.88;
    GLU.atomInfo["CD"].radius = 1.61;
    GLU.atomInfo["OE1"].radius = 1.46;
    GLU.atomInfo["OE2"].radius = 1.42;

    auto& LYS = types[(AminoAcid)"LYS"];
    LYS.mass = 128.094963050;
    LYS.atomInfo["CB"].radius = 1.88;
    LYS.atomInfo["CG"].radius = 1.88;
    LYS.atomInfo["CD"].radius = 1.88;
    LYS.atomInfo["CE"].radius = 1.88;
    LYS.atomInfo["NZ"].radius = 1.64;

    auto& ARG = types[(AminoAcid)"ARG"];
    ARG.mass = 156.101111050;
    ARG.atomInfo["CB"].radius = 1.88;
    ARG.atomInfo["CG"].radius = 1.88;
    ARG.atomInfo["CD"].radius = 1.88;
    ARG.atomInfo["NE"].radius = 1.64;
    ARG.atomInfo["CZ"].radius = 1.61;
    ARG.atomInfo["NH1"].radius = 1.64;
    ARG.atomInfo["NH2"].radius = 1.64;

    auto& SER = types[(AminoAcid)"SER"];
    SER.mass = 87.032028435;
    SER.atomInfo["CB"].radius = 1.88;
    SER.atomInfo["OG"].radius = 1.46;

    auto& THR = types[(AminoAcid)"THR"];
    THR.mass = 101.047678505;
    THR.atomInfo["CB"].radius = 1.88;
    THR.atomInfo["OG1"].radius = 1.46;
    THR.atomInfo["CG2"].radius = 1.88;

    auto& TYR = types[(AminoAcid)"TYR"];
    TYR.mass = 163.063328575;
    TYR.atomInfo["CB"].radius = 1.88;
    TYR.atomInfo["CG"].radius = 1.61;
    TYR.atomInfo["CD1"].radius = 1.76;
    TYR.atomInfo["CD2"].radius = 1.76;
    TYR.atomInfo["CE1"].radius = 1.76;
    TYR.atomInfo["CE2"].radius = 1.76;
    TYR.atomInfo["CZ"].radius = 1.61;
    TYR.atomInfo["OH"].radius = 1.46;

    auto& HIS = types[(AminoAcid)"HIS"];
    HIS.mass = 137.058911875;
    HIS.atomInfo["CB"].radius = 1.88;
    HIS.atomInfo["CG"].radius = 1.61;
    HIS.atomInfo["ND1"].radius = 1.64;
    HIS.atomInfo["CD2"].radius = 1.76;
    HIS.atomInfo["CE1"].radius = 1.76;
    HIS.atomInfo["NE2"].radius = 1.64;

    auto& ASN = types[(AminoAcid)"ASN"];
    ASN.mass = 114.042927470;
    ASN.atomInfo["CB"].radius = 1.88;
    ASN.atomInfo["CG"].radius = 1.61;
    ASN.atomInfo["OD1"].radius = 1.42;
    ASN.atomInfo["ND2"].radius = 1.64;

    auto& TRP = types[(AminoAcid)"TRP"];
    TRP.mass = 186.079312980;
    TRP.atomInfo["CB"].radius = 1.88;
    TRP.atomInfo["CG"].radius = 1.61;
    TRP.atomInfo["CD1"].radius = 1.76;
    TRP.atomInfo["NE1"].radius = 1.61;
    TRP.atomInfo["CE2"].radius = 1.64;
    TRP.atomInfo["CD2"].radius = 1.61;
    TRP.atomInfo["CE3"].radius = 1.76;
    TRP.atomInfo["CZ3"].radius = 1.76;
    TRP.atomInfo["CZ2"].radius = 1.76;
    TRP.atomInfo["CH2"].radius = 1.76;

    auto& ALA = types[(AminoAcid)"ALA"];
    ALA.mass = 71.037113805;
    ALA.atomInfo["CB"].radius = 1.88;

    auto& LEU = types[(AminoAcid)"LEU"];
    LEU.mass = 113.084064015;
    LEU.atomInfo["CB"].radius = 1.88;
    LEU.atomInfo["CG"].radius = 1.88;
    LEU.atomInfo["CD1"].radius = 1.88;
    LEU.atomInfo["CD2"].radius = 1.88;

    auto& GLY = types[(AminoAcid)"GLY"];
    GLY.mass = 57.021463735;

    for (auto& [acid, acidInfo]: types) {
        acidInfo.mass *= au;
        for (auto& [atomType, atomInfo]: acidInfo.atomInfo) {
            acidInfo.heavyAtoms.insert(atomType);
            atomInfo.radius *= angstrom;
        }
    }

    return types;
}

AminoAcidInfo const &AminoAcid::info() const {
    static auto allInfo = createResData();
    return allInfo.at(*this);
}
