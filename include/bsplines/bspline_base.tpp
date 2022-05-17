#include <eigen3/Eigen/Dense>

#include "bspline_base.hpp"

namespace spline {


template<typename T, int K>
BSplineBase<T,K>::~BSplineBase() {}

template<typename T, int K>
typename BSplineBase<T,K>::VectorXT BSplineBase<T,K>::operator()(double t) {
    int i{0};
    if (t == knot_pts_(ctrl_pts_.rows()+K)) {
        i = ctrl_pts_.rows()-1;
    } else {
        for(; i != knot_pts_.rows()-1; ++i) {
            if (t >= knot_pts_(i) && t < knot_pts_(i+1))
                break;
        }
    }

    T u = (t - knot_pts_(i))/(knot_pts_(i+1) - knot_pts_(i));
    VectorXT u_vec = VectorXT(K+1);
    for(int j{0}; j != K+1; ++j)
        u_vec(j) = pow(u, j);

    MatrixXT V = ctrl_pts_.block(i-K, 0, K+1, ctrl_pts_.cols());

    return u_vec.transpose() * M_[i] * V;
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
