#pragma once
#include "NonlocalForce.hpp"
#include "../kernels/LennardJones.hpp"
#include "../kernels/SidechainLJ.hpp"
#include "../data/Types.hpp"
#include "../data/Chains.hpp"
#include "../utils/AminoAcid.hpp"

namespace mdk {
    /**
     * An object representing either of the three lambda functions. It is
     * theoretically a Gaussian function, but we implement it here in terms of
     * computationally faster functions (either cosine or an "algebraic version"
     * thereof, with a cutoff at the far ends).
     */
    class LambdaPeak {
    public:
        /**
         * The "width" of the lambda function.
         */
        double alpha;

        /**
         * The "center" of the lambda function.
         */
        double psi0;

        /**
         * Whether we should use the cosine version of the function.
         */
        bool cosineVersion;

        /**
         * Check whether a value $psi$ is in the support of the function.
         * @param psi Value to test.
         * @return true if $lambda(psi) != 0$, false otherwise.
         */
        bool supp(double psi) const;

        /**
         * Compute the value of the function.
         * @param psi Argument, i.e. the improper dihedral angle
         * @param L Value of the lambda function at \p psi,
         * @param dL_dpsi Derivative of L wrt \p psi
         */
        void eval(double psi, double& L, double& dL_dpsi) const;
    };

    /**
     * The "pseudo-improper-dihedral" custom potential. A more computationally
     * expensive but more accurate and Hamiltonian alternative to the
     * quasi-adiabatic potential.
     */
    class PseudoImproperDihedral: public NonlocalForce {
    private:
        /**
         * Derive the angle data (specifically the angles and derivatives)
         * pertaining to the pair.
         * @param pair Pair of residues.
         * @param psi An array, where psi[0] is the improper dihedral angle
         * between i_1-1, i_1, i_1+1 and i_2, and psi[1] is the improper dihedral
         * angle between i_2-1, i_2, i_2+1 and i_1.
         * @param dpsi_dr An array, where dpsi_dr[0][0..2] are the derivatives
         * d psi_1/d r_{i_1-1 .. i_1+1} and dspi_dr[0][3..6] are the derivatives
         * d psi_1/d r_{i_2-1 .. i_2+1}; dpsi_dr[1][0..3] is defined as for 0,
         * but with psi_2 instead of psi_1 (order of variables is not changed).
         */
        void deriveAngles(vl::PairInfo const& pair, double psi[2],
            Vector dpsi_dr[2][6]) const;

        /**
         * A const pointer to the types of the residues.
         */
        Types const* types = nullptr;

        /**
         * A const pointer to the chain data of the model.
         */
        Chains const* seqs = nullptr;

    public:
        /**
         * Lambda function corresponding to the bb+ part of the potential.
         */
        LambdaPeak bb_pos;

        /**
         * Lennard-Jones potential corresponding to the bb+ part of the
         * potential.
         */
        LennardJones bb_pos_lj;

        /**
         * Lambda function corresponding to the bb- part of the potential.
         */
        LambdaPeak bb_neg;

        /**
         * Lennard-Jones potential corresponding to the bb- part of the
         * potential.
         */
        LennardJones bb_neg_lj;

        /**
         * Lambda function corresponding to the ss part of the potential.
         */
        LambdaPeak ss;

        /**
         * An array of sidechain L-J potentials for each pair of types of
         * amino acids interacting.
         */
        SidechainLJ ss_ljs[AminoAcid::N][AminoAcid::N];

        /**
         * Bind the object to the simulation. Here we also fetch from it the
         * types and chains.
         * @param simulation Simulation to bind to.
         */
        void bind(Simulation& simulation) override;

        /**
         * Asynchronous part of the force computation.
         * @param dynamics Dynamics object to add potential energy and
         * forces to.
         */
        void asyncPart(Dynamics &dynamics) override;

        /**
         * An action to be performed when the Verlet list is updated. We
         * essentially copy the residues but exclude the pairs within bond
         * distance of less than 4, or at the ends of their chains.
         */
        void vlUpdateHook() override;

    protected:
        /**
         * Generate spec for the Verlet list.
         * @return Generated spec.
         */
        vl::Spec spec() const override;

    private:
        /**
         * Local copy of the Verlet list with some pairs potentially excluded
         * (see \p vlUpdateHook for details).
         */
        Pairs pairs;
    };
}
