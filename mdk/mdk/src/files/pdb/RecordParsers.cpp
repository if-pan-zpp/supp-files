#include "RecordParsers.hpp"
#include "utils/Units.hpp"
#include "Field.hpp"
using namespace mdk::pdb;
using namespace std;

std::ostream& RecordParser::write(std::ostream& os, Record const& other) {
    if (record.index() == other.index()) {
        record = other;

        std::string line(80, ' ');
        for (auto const& field: fields) {
            field->write(line);
        }
        os << line;
    }
    return os;
}

Record RecordParser::tryParse(const string &s) const {
    try {
        for (auto const& field: fields) {
            field->read(s);
        }
        return record;
    }
    catch (exception const& e) {
        return std::monostate();
    }
}

RemarkParser::RemarkParser() {
    record = Remark();
    auto& remark = get<Remark>(record);

    fields = {
        make_shared<Literal>(1, 6, "REMARK"),
        make_shared<Integer>(8, 10, remark.number),
        make_shared<String>(12, 80, remark.text)
    };
}

AtomParser::AtomParser() {
    record = Atom();
    auto& atom = get<Atom>(record);

    fields = {
        make_shared<Literal>(1, 6, "ATOM  "),
        make_shared<Integer>(7, 11, atom.serialNum, -1),
        make_shared<String>(13, 16, atom.atomName),
        make_shared<Char>(17, atom.altLocation),
        make_shared<String>(18, 20, atom.residueName),
        make_shared<Char>(22, atom.chainID),
        make_shared<Integer>(23, 26, atom.residueSeqNum, -1),
        make_shared<Char>(27, atom.insertionCode),
        make_shared<Real>(31, 38, 8, 3, atom.pos.x(), angstrom),
        make_shared<Real>(39, 46, 8, 3, atom.pos.y(), angstrom),
        make_shared<Real>(47, 54, 8, 3, atom.pos.z(), angstrom),
        make_shared<Real>(55, 60, 6, 2, atom.occupancy),
        make_shared<Real>(61, 66, 6, 2, atom.tempFactor),
        make_shared<String>(77, 78, atom.element, Trim | Right),
        make_shared<String>(79, 80, atom.charge)
    };
}

SSBondParser::SSBondParser() {
    record = SSBond();
    auto& ssbond = get<SSBond>(record);

    fields = {
        make_shared<Literal>(1, 6, "SSBOND"),
        make_shared<Integer>(8, 10, ssbond.serialNum, -1),
        make_shared<Literal>(12, 14, "CYS"),
        make_shared<Char>(16, ssbond.res[0].chainId),
        make_shared<Integer>(18, 21, ssbond.res[0].residueSeqNum, -1),
        make_shared<Char>(22, ssbond.res[0].insertionCode),
        make_shared<Char>(30, ssbond.res[1].chainId),
        make_shared<Integer>(32, 35, ssbond.res[1].residueSeqNum, -1),
        make_shared<Char>(36, ssbond.res[1].insertionCode),
        make_shared<SymOp>(60, 65, ssbond.res[0].symmetryOp),
        make_shared<SymOp>(67, 72, ssbond.res[1].symmetryOp),
        make_shared<Real>(74, 78, 5, 2, ssbond.dist0)
    };
}

Cryst1Parser::Cryst1Parser() {
    record = Cryst1();
    auto& cryst1 = get<Cryst1>(record);

    fields = {
        make_shared<Literal>(1, 6, "CRYST1"),
        make_shared<Real>(7, 15, 9, 3, cryst1.size.x(), angstrom),
        make_shared<Real>(16, 24, 9, 3, cryst1.size.y(), angstrom),
        make_shared<Real>(25, 33, 9, 3, cryst1.size.z(), angstrom),
        make_shared<Real>(34, 40, 7, 2, cryst1.angles.x(), degree),
        make_shared<Real>(41, 47, 7, 2, cryst1.angles.y(), degree),
        make_shared<Real>(48, 54, 7, 2, cryst1.angles.z(), degree),
        make_shared<String>(56, 66, cryst1.spaceGroup),
        make_shared<Integer>(67, 70, cryst1.z)
    };
}

LinkParser::LinkParser() {
    record = Link();
    auto& link = get<Link>(record);

    fields = {
        make_shared<Literal>(1, 6, "LINK  "),
        make_shared<String>(13, 16, link.res[0].atomName),
        make_shared<Char>(17, link.res[0].altLocation),
        make_shared<String>(18, 20, link.res[0].residueName),
        make_shared<Char>(22, link.res[0].chainId),
        make_shared<Integer>(23, 26, link.res[0].residueSeqNum),
        make_shared<Char>(27, link.res[0].insertionCode),
        make_shared<SymOp>(60, 65, link.res[0].symmetryOp),
        make_shared<String>(43, 46, link.res[1].atomName),
        make_shared<Char>(47, link.res[1].altLocation),
        make_shared<String>(48, 50, link.res[1].residueName),
        make_shared<Char>(52, link.res[1].chainId),
        make_shared<Integer>(53, 56, link.res[1].residueSeqNum),
        make_shared<Char>(57, link.res[1].insertionCode),
        make_shared<SymOp>(67, 72, link.res[1].symmetryOp),
        make_shared<Real>(74, 78, 5, 2, link.linkLength, angstrom),
    };
}

ModelParser::ModelParser() {
    record = Model();
    auto& model = get<Model>(record);

    fields = {
        make_unique<Literal>(1, 6, "MODEL "),
        make_unique<Integer>(11, 14, model.serialNum)
    };
}

EndmdlParser::EndmdlParser() {
    record = Endmdl();
    fields = {
        make_unique<Literal>(1, 6, "ENDMDL")
    };
}

EndParser::EndParser() {
    record = End();
    fields = {
        make_unique<Literal>(1, 6, "END   ")
    };
}

TerParser::TerParser() {
    record = Ter();
    auto& ter = get<Ter>(record);

    fields = {
        make_unique<Literal>(1, 6, "TER   "),
        make_unique<Integer>(7, 11, ter.serialNum, -1),
        make_unique<String>(18, 20, ter.residueName),
        make_unique<Char>(22, ter.chainId),
        make_unique<Integer>(23, 26, ter.residueSeqNum, -1),
        make_unique<Char>(27, ter.insertionCode)
    };
}