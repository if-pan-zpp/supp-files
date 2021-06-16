#pragma once
#include "../../utils/Units.hpp"
#include <Eigen/Core>
#include <string>
#include <variant>

/**
 * A namespace containing PODs of PDB records. For reference, see PDB format
 * report.
 */
namespace mdk::pdb::records {
    class Remark {
    public:
        int number = 6;
        std::string text = "";
    };

    class Atom {
    public:
        int serialNum = 0; // Zero-indexed here.
        std::string atomName;
        char altLocation = ' ';
        std::string residueName;
        char chainID = 'A';
        int residueSeqNum = 0; // Zero-indexed here.
        char insertionCode = ' ';
        Eigen::Vector3d pos;
        double occupancy = 1.0;
        double tempFactor = 0.0;
        std::string element;
        std::string charge;
    };

    class SSBond {
    public:
        int serialNum = 0;

        struct PerResidueData {
            char chainId = 'A';
            int residueSeqNum = 0; // Zero-indexed here.
            char insertionCode = ' ';
            std::string symmetryOp;
        };
        PerResidueData res[2];

        double dist0 = 0.0;
    };

    class Cryst1 {
    public:
        Eigen::Vector3d size;
        Eigen::Vector3d angles = { 90*degree, 90*degree, 90*degree };
        std::string spaceGroup;
        int z = 0;
    };

    class Ter {
    public:
        int serialNum = 0; // Zero-indexed here.
        std::string residueName;
        char chainId = 'A';
        int residueSeqNum = 0; // Zero-indexed here.
        char insertionCode = ' ';
    };

    class Link {
    public:
        struct PerResidueData {
            std::string atomName;
            char altLocation = ' ';
            std::string residueName;
            char chainId = 'A';
            int residueSeqNum;
            char insertionCode = ' ';
            std::string symmetryOp;
        };
        PerResidueData res[2];
        double linkLength;
    };

    class Model {
    public:
        int serialNum = 0; // Zero-indexed here.
    };

    class Endmdl {
    public:
    };

    class End {
    public:
    };

    /**
     * A variant type denoting either a record of some type, or no record
     * \p std::monostate in the case of a parsing failure.
     */
    using Record = std::variant<std::monostate,
        Remark, Atom, SSBond, Cryst1, End, Link,
        Model, Endmdl, Ter>;
}