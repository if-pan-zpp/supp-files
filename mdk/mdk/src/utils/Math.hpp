#pragma once
#include <Eigen/Dense>

namespace mdk {
    using VRef = Eigen::Vector3d const&;

    inline double angle(VRef v1, VRef v2, VRef v3) {
        auto u1 = v2 - v1, u2 = v2 - v3;
        return acos(u1.dot(u2) / (u1.norm() * u2.norm()));
    }

    inline double dihedral(VRef v1, VRef v2, VRef v3, VRef v4) {
        auto u1 = v2 - v1, u2 = v3 - v2, u3 = v4 - v3;
        auto d1 = u1.cross(u2), d2 = u2.cross(u3);

        if (d1.isZero() || d2.isZero()) return 0.0;
        double phi = acos(d1.dot(d2) / (d1.norm() * d2.norm()));
        if (d1.dot(u3) < 0.) phi = -phi;
        return phi;
    }

    double asphericity(Vectors const& v) {
        Vector mean = v.vectorwise().mean();

        Eigen::Matrix3d mit;
        double crg;
        for (size_t i = 0; i < v.size(); ++i) {
            Vector diff = v.col(i) - mean;
            crg += diff.squaredNorm();

            Vector sq_diff = diff.array() * diff.array();
            mit(0,0) += sq_diff(1) + sq_diff(2);
            mit(1,1) += sq_diff(2) + sq_diff(0);
            mit(2,2) += sq_diff(0) + sq_diff(1);
            mit(0,1) -= diff(0) * diff(1);
            mit(1,2) -= diff(1) * diff(2);
            mit(0,2) -= diff(0) * diff(2);
        }
        mit(1,0) = mit(0,1);
        mit(2,1) = mit(1,2);
        mit(2,0) = mit(0,2);
        crg = sqrt(crg / v.size());

        auto eigenvaluesC = mit.eigenvalues();
        double eigenvalues[3] = {
            sqrt(eigenvaluesC(0).real() / v.size()),
            sqrt(eigenvaluesC(1).real() / v.size()),
            sqrt(eigenvaluesC(2).real() / v.size())
        };

        std::sort(eigenvalues, eigenvalues + 3);
        return eigenvalues[1] - 0.5 * (eigenvalues[0] + eigenvalues[2]);
    }
}
