#include "model/Model.hpp"
#include "utils/AminoAcid.hpp"
#include <Eigen/Geometry>
using namespace std;
using namespace mdk;

Model::Residue& Model::addResidue(Chain *chain) {
    auto& res = residues.emplace_back();
    res.idx = residues.size()-1;
    if (chain) {
        res.chainIdx = chain->idx;
        if (chain->start < 0) chain->start = res.idx;
        chain->end = res.idx+1;
    }
    else {
        res.chainIdx = -1;
    }
    ++n;
    return res;
}

Model::Chain& Model::addChain() {
    auto& chain = chains.emplace_back();
    chain.idx = chains.size()-1;
    chain.start = chain.end = -1;
    return chain;
}

Model::Contact &Model::addContact() {
    auto& cont = contacts.emplace_back();
    cont.idx = contacts.size()-1;
    return cont;
}

Model::StructuredPart& Model::addSP() {
    auto& sp = structuredParts.emplace_back();
    sp.idx = structuredParts.size()-1;
    return sp;
}

void Model::morphIntoLine() {
    double tether0 = 3.8*angstrom;

    Eigen::Vector3d dir = { 0.0, 0.0, tether0 };
    for (auto& res: residues) {
        res.r = res.idx * dir;
    }
}

static Eigen::Vector3d sampleBox(Eigen::AlignedBox3d const& box, Random& rand) {
    Eigen::Vector3d v;
    for (int i = 0; i < 3; ++i) {
        v[i] = rand.uniform(box.min()[i], box.max()[i]);
    }
    return v;
}

static Eigen::Vector3d sampleSphere(Random& rand) {
    Eigen::Vector3d v;
    for (int i = 0; i < 3; ++i) {
        v[i] = rand.normal();
    }
    return v.normalized();
}

/**
 * Find a vector perpendicular to \p v.
 * @param v
 * @return
 */
static Eigen::Vector3d findPerp(Eigen::Vector3d const& v) {
    auto perp1 = v.cross(Eigen::Vector3d::UnitX());
    auto perp2 = v.cross(Eigen::Vector3d::UnitY());
    auto& perp = perp1.norm() > perp2.norm() ? perp1 : perp2;
    return perp.normalized();
}

void Model::morphIntoSAW(Random& rand, bool useTop, double density,
        double minDist) {
    Eigen::Vector3d minCorner, maxCorner;
    Eigen::AlignedBox3d box;

    double tether0 = 3.8*angstrom;
    auto pairs = nonlocalPairs();
    auto state0 = residues;

    if (density > 0) {
        auto vol = residues.size() / (density * atom);
        auto a = 0.5 * pow(vol, 1.0/3.0);
        minCorner = -a * Eigen::Vector3d { 1.0, 1.0, 1.0 };
        maxCorner = -minCorner;
    }
    else {
        minCorner = maxCorner = Eigen::Vector3d { 0.0, 0.0, 0.0 };
    }
    box = Eigen::AlignedBox3d(minCorner, maxCorner);
    if (useTop) top.setCell(box.max() - box.min());

    for (int attempt = 0; attempt < 9000; ++attempt) {
        bool failed = false;

        for (auto const& chain: chains) {
            auto pos = sampleBox(box, rand);
            auto dir = sampleSphere(rand);

            for (int i = chain.start; i < chain.end; ++i) {
                residues[i].r = pos;
                pos += tether0 * dir;

                double spread = rand.uniform(M_PI/12.0, M_PI/2.0);
                auto spreadRot = Eigen::AngleAxisd(spread, findPerp(dir));

                double around = rand.uniform(-M_PI, M_PI);
                auto aroundRot = Eigen::AngleAxisd(around, dir);

                dir = aroundRot * spreadRot * dir;
                dir.normalize();
            }
        }

        for (auto& res: residues) {
            if (res.chainIdx < 0)
                res.r = sampleBox(box, rand);
        }

        for (auto const& [res1, res2]: pairs) {
            Eigen::Vector3d dx = res1->r - res2->r;
            if (useTop) top.fix(dx);

            if (dx.norm() < minDist) {
                failed = true;
                break;
            }
        }

        if (!failed) return;
    }

    residues = state0;
}

Vectors genSAW(int n, double bondLen, Random& rand) {
    vector<double> phi(n);   // around
    vector<double> theta(n); // spread

    phi[0] = 3.0 * M_PI / 2.;
    phi[1] = M_PI;
    theta[0] = 0.;
    theta[1] = rand.uniform(0.0, M_PI / 3.0);
    for (int i = 2; i + 1 < n; ++i) {
        phi[i] = rand.uniform(0.0, 2.0 * M_PI);
        theta[i] = rand.uniform(0.0, M_PI / 3.0);
    }

    double dir_theta = M_PI / 2.0 - acos(rand.uniform(-1.0, 1.0));
    double dir_phi = rand.uniform(0.0, 2.0 * M_PI); 

    auto dir = Eigen::Vector3d{cos(dir_theta)*cos(dir_phi),
                               -cos(dir_theta)*sin(dir_phi),
                               sin(dir_theta)};

    auto rot_matrix = Eigen::AngleAxisd::Identity();

    Vector pos = Vector::Zero();

    Vectors res(n);
    res[0] = pos;
    for(int i = 1; i < n; i++) {
        rot_matrix = rot_matrix *
            Eigen::AngleAxisd(phi[i - 1], Eigen::Vector3d::UnitX()) *
            Eigen::AngleAxisd(theta[i - 1], Eigen::Vector3d::UnitZ());
        pos += rot_matrix * dir * bondLen;
        res[i] = pos;
    }

    auto flip = Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitX());
    for (int i = 0; i < n; ++i) res[i] = flip * res[i];
    return res;
}

// SAW implementation that mirrors cg.f:confstart
// Assumes that each residue is in some chain
// and chains vector is sorted
void Model::legacyMorphIntoSAW(Random& rand, bool useTop, double density,
                               double minDist, bool nativeBondLen) {

    Eigen::Vector3d minCorner, maxCorner;
    Eigen::AlignedBox3d box;

    double tether0 = (float) 3.8*angstrom;
    if (nativeBondLen) {
        double avgBondLen = 0.0;
        for (Chain const& chain : chains) {
            for (int i = chain.start; i + 1 < chain.end; ++i) {
                avgBondLen += (residues[i + 1].r - residues[i].r).norm();
            }
        }
        avgBondLen /= (int) residues.size() - 1;
        tether0 = avgBondLen;
    }
    
    auto pairs = nonlocalPairs();
    auto state0 = residues;

    if (density > 0) {
        auto vol = residues.size() / (density * atom);
        auto a = 0.5 * pow(vol, 1.0/3.0);
        minCorner = -a * Eigen::Vector3d { 1.0, 1.0, 1.0 };
        maxCorner = -minCorner;
    }
    else {
        minCorner = maxCorner = Eigen::Vector3d { 0.0, 0.0, 0.0 };
    }
    box = Eigen::AlignedBox3d(minCorner, maxCorner);
    if (useTop) top.setCell(box.max() - box.min());

    for (auto& res: residues) {
        // this SAW implementation doesn't support residues without chains
        assert (res.chainIdx >= 0);
    }

    for (auto const& chain: chains) {
        int attemptsLeft = 9000;
        bool success = false;
        while (not success && attemptsLeft > 0) {
            attemptsLeft--;

            Vectors saw = genSAW(chain.end - chain.start, tether0, rand);
            auto displacement = sampleBox(box, rand);
            for (int i = chain.start; i < chain.end; ++i) {
                residues[i].r = saw[i - chain.start] + displacement;
            }
            success = true;

            double minDistSquared = minDist * minDist;
            for (int i = 0; i < chain.end && success; ++i) {
                for (int j = i + 3; j < chain.end; ++j) {
                    Eigen::Vector3d vec = residues[i].r - residues[j].r;
                    if (useTop) top.fix(vec);

                    if (vec.squaredNorm() < minDistSquared) {
                        success = false;
                        break;
                    }
                }
            }
        } 

        assert (success /* morphing into SAW failed */);
    }
}

void Model::exactLegacySAW(Random &rand, bool useTop, double density,
    double minDist, double bond) {

    auto minDistSq = pow(minDist, 2.0);

    Eigen::AlignedBox3d box;
    if (density > 0) {
        auto startvolume = residues.size() / density;
        auto startboxsize = 0.5 * pow(startvolume, 1.0/3.0);
        box.min() = Vector::Constant(-startboxsize);
        box.max() = -box.min();

        if (useTop) {
            top.setCell(box.max() - box.min());
        }
    }
    else {
        if (useTop) throw;

        box.min() = box.max() = Vector::Zero();
    }

    for (auto& chain: chains) {
        int n = chain.end - chain.start;
        for (int kter = 0; kter < 9000; ++kter) {
            std::vector<double> phi(n), theta(n);

            auto pi = acos(-1.0);

            phi[0] = pi / 2.0;
            theta[0] = 0.0;

            phi[1] = 0.0;
            theta[1] = rand.uniform(0.0, M_PI/3.0);

            for (int i = 2; i < n - 1; ++i) {
                phi[i] = rand.uniform(-M_PI, M_PI);
                theta[i] = rand.uniform(0.0, M_PI/3.0);
            }

            std::vector<Eigen::Matrix3d> T(n);
            for (int i = 0; i < n - 1; ++i) {
                auto ct = cos(theta[i]), st = sin(theta[i]),
                     cp = cos(phi[i]), sp = sin(phi[i]);

                T[i] << ct,      st,       0.0,
                        st * cp, -ct * cp,  sp,
                        st * sp, -ct * sp, -cp;
            }

            auto theta0 = acos(-rand.uniform(-1.0, 1.0));
            auto phi0 = rand.uniform(0.0, 2.0 * pi);

            std::vector<Vector> R(n);
            for (int i = 0; i < n - 1; ++i) {
                Vector r;
                r << bond * sin(theta0) * cos(phi0),
                     bond * sin(theta0) * sin(phi0),
                     bond * cos(theta0);

                for (int j = i; j >= 0; --j) {
                    R[i] = T[j] * r;
                    r = R[i];
                }
            }

            Vector ran = sampleBox(box, rand);
            for (int i = 0; i < n; ++i) {
                Vector r = ran;
                for (int j = 0; j < i; ++j) {
                    r += R[j];
                }
                residues[i].r = r;
            }

            bool nonIntersecting = true;
            for (int i = 0; i < n - 3 && nonIntersecting; ++i) {
                auto r1 = residues[i].r;
                for (int j = i + 3; j < n && nonIntersecting; ++j) {
                    auto r2 = residues[j].r;
                    Vector dx = !useTop ? (r2 - r1) : top(r2 - r1);
                    if (dx.squaredNorm() < minDistSq) {
                        nonIntersecting = false;
                    }
                }
            }

            if (nonIntersecting) break;
        }
    }
}

void Model::initVelocity(Random& rand, double temperature, bool useMass) {
    int n = residues.size();
    vector<double> masses(n, 1.0 * f77mass);

    if (useMass) {
        assert(0);
    }

    for (int i = 0; i < n; ++i) {
        Vector& v = residues[i].v;
        for (int d = 0; d < 3; ++d) {
            v(d) = rand.uniform(-1.0, 1.0);
        }
        v.normalize();
    }

    // Shifting velocities so that total momentum is 0

    double totMass = 0.0;
    Vector totMomentum = Vector::Zero();

    for (int i = 0; i < n; ++i) {
        totMass += masses[i];
        totMomentum += masses[i] * residues[i].v;
    }

    for (int i = 0; i < n; ++i) {
        residues[i].v -= totMomentum * (1.0 / totMass);
    }

    // Scaling velocities to desired temperature

    double avgKinEnergy = 0.0;
    for (int i = 0; i < n; ++i) {
        avgKinEnergy += masses[i] * residues[i].v.squaredNorm() / 2.0;
    }
    avgKinEnergy /= n;

    double ratio = 1.5 * temperature / avgKinEnergy;
    for (int i = 0; i < n; ++i) {
        residues[i].v *= sqrt(ratio);
    }
}

Model::StructuredPart& Model::addContactMap(cmap::ContactMap const& contactMap) {
    auto& sp = addSP();
    sp.off = contactMap.offset;
    sp.len = contactMap.len;
    sp.angle = contactMap.angle;
    sp.dihedral = contactMap.dihedral;
    return sp;
}

void Model::addCMapContacts(cmap::ContactMap const& contactMap, Chain &chain) {
    for (auto const& cmapCont: contactMap.contacts) {
        auto& modelCont = addContact();
        modelCont.type = ContactType(ContactTypeIdx::NAT);
        modelCont.dist0 = cmapCont.dist0;
        for (int i = 0; i < 2; ++i) {
            auto resIdx = chain.start + contactMap.offset + cmapCont.res[i];
            modelCont.res[i] = resIdx;
        }
    }

    auto sp = addContactMap(contactMap);
    chain.structuredParts.push_back(sp.idx);
}

vector<pair<Model::Residue*, Model::Residue*>> Model::nonlocalPairs() {
    vector<pair<Residue*, Residue*>> pairs;
    for (auto& res1: residues) {
        for (auto& res2: residues) {
            if (res1.chainIdx == res2.chainIdx && res1.chainIdx >= 0)
                continue;
            if (abs(res1.idx - res2.idx) < 3)
                continue;

            pairs.emplace_back(&res1, &res2);
        }
    }

    return pairs;
}

void Model::useAverageMasses() {
    double sum = 0.0;
    for (auto& acid: AminoAcid::aminoAcids())
        sum += acid.info().mass;
    sum /= AminoAcid::N;

    for (auto& res: residues)
        res.mass = sum;
}
