#pragma once
#include <Eigen/Core>

namespace mdk {
    using Vector = Eigen::Vector3d;
    using VRef = Vector const&;

    using VectorBase = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>;

    /**
     * This class is a sort of a decorator over standard Eigen::Matrix3Xd. The
     * chief difference is that this class offers an interface more appropriate
     * for dealing with it as though it were a list of vectors (for example, one
     * needs not use m.col(i) to get t'th vector in the matrix. Also, now size()
     * corresponds to the number of vectors and not the elements of the array.
     * The array itself is not (easily) resizeable.
     */
    class Vectors: public VectorBase {
    public:
        Vectors() = default;

        /**
         * Create an unitialized list of n vectors.
         * @param n Number of vectors
         */
        explicit Vectors(int n):
            VectorBase(3, n) {};

        /**
         * Create an initialized list of n vectors.
         * @param n Number of vectors
         * @param init Initial value
         */
        Vectors(int n, Vector const& init) {
            resize(3, n);
            colwise() = init;
        }

        /**
         * Returns an iterator over constituent vectors.
         * @return An interator over constituent vectors.
         */
        auto vectorwise() {
            return colwise();
        }

        /**
         * Returns a const iterator over constituent vectors.
         * @return A const interator over constituent vectors.
         */
        auto vectorwise() const {
            return colwise();
        }

        /**
         *
         * @return A number of vectors.
         */
        int size() const {
            return cols();
        }

        /**
         * Access to vector.
         * @param i Index of a vector to access
         * @return A slice of a matrix corresponding to i'th vector.
         */
        inline auto operator[](int i) {
            return col(i);
        }

        /**
         * Const access to vector.
         * @param i Index of a vector to access
         * @return A const slice of a matrix corresponding to i'th vector.
         */
        inline auto operator[](int i) const {
            return col(i);
        }
    };
}
