#include <eigen3/Eigen/Dense>

#include "bspline_base.hpp"

namespace spline {


template<typename T, int K>
BSplineBase<T,K>::~BSplineBase() {}


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
