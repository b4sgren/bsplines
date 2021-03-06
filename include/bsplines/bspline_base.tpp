#include <eigen3/Eigen/Dense>

#include "bspline_base.hpp"

namespace spline {


template<typename T, int K>
BSplineBase<T,K>::~BSplineBase() {}

template<typename T, int K>
typename BSplineBase<T,K>::VectorXT BSplineBase<T,K>::operator()(double t) {
    int i{getI(t)};

    T u = (t - knot_pts_(i))/(knot_pts_(i+1) - knot_pts_(i));
    VectorXT u_vec = VectorXT(K+1);
    for(int j{0}; j != K+1; ++j)
        u_vec(j) = pow(u, j);

    MatrixXT V = ctrl_pts_.block(i-K, 0, K+1, ctrl_pts_.cols());

    return u_vec.transpose() * M_[i] * V;
}

template<typename T, int K>
typename BSplineBase<T,K>::VectorXT BSplineBase<T,K>::evaluateDerivative
                                                        (double t, int n) {
    int i{getI(t)};

    T u = (t - knot_pts_(i))/(knot_pts_(i+1) - knot_pts_(i));
    VectorXT u_vec = VectorXT(K+1);
    for(int j{0}; j != K+1; ++j) {
        if (j-n >= 0)
            u_vec(j) = factorial(j)/factorial(std::max(0, j-n)) * pow(u, j-n);
        else
            u_vec(j) = T(0.0);
    }

    MatrixXT V = ctrl_pts_.block(i-K, 0, K+1, ctrl_pts_.cols());

    return u_vec.transpose() * M_[i] * V;
}

template<typename T, int K>
bool BSplineBase<T,K>::addControlPoint(const Eigen::Ref<const VectorXT> &pt) {
    if (pt.rows() != ctrl_pts_.cols())
        return false;

    ctrl_pts_.conservativeResize(ctrl_pts_.rows()+1, ctrl_pts_.cols());
    ctrl_pts_.row(ctrl_pts_.rows()-1) = pt.transpose();

    knot_pts_.conservativeResize(knot_pts_.rows()+1);
    knot_pts_.tail(K+1) = VectorXT::Ones(K+1) * (ctrl_pts_.rows()-K);

    M_.resize(knot_pts_.rows());
    bool success = updateBasisMatrices();

    return success;
}

template<typename T, int K>
double BSplineBase<T,K>::factorial(int n) {
    double val{1.0};
    for (int i{1}; i != n+1; ++i)
        val *= i;

    return val;
}

template<typename T, int K>
int BSplineBase<T,K>::getI(double t) {
    // int i = static_cast<int>(std::floor((t - tr_.start)/dt_)); // Generalized version
    int i = static_cast<int>(t) + K; // When start = 0 and dt = 1
    if (i == ctrl_pts_.rows())
        i -= 1;
    return i;
}

template<typename T, int K>
int BSplineBase<T,K>::degree() const {
    return K;
}

template<typename T, int K>
typename BSplineBase<T,K>::MatrixXT BSplineBase<T,K>::getControlPoints() const {
    return ctrl_pts_;
}

template<typename T, int K>
typename BSplineBase<T,K>::VectorXT BSplineBase<T,K>::getKnotPoints() const {
    return knot_pts_;
}

}  // namespace spline
