#ifndef BSPLINES_BSPLINE_BASE_HPP_
#define BSPLINES_BSPLINE_BASE_HPP_

#include <eigen3/Eigen/Dense>
#include <vector>

#include <iostream>

namespace spline {

template <typename T, int K>
class BSplineBase {
 public:
    using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    BSplineBase() {}
    // Doesn't like this in the tpp file
    explicit BSplineBase(const MatrixXT ctrl_pts) : ctrl_pts_(ctrl_pts) {
        int cnt = ctrl_pts_.rows();
        knot_pts_ = VectorXT(ctrl_pts_.rows() + K + 1);
        knot_pts_.setZero();
        knot_pts_.tail(K) = VectorXT::Ones(K) * (cnt-K);
        VectorXT temp = VectorXT::LinSpaced(cnt-K+1, 0, cnt-K);
        knot_pts_.segment(K, cnt-K+1) = temp;
    }
    virtual ~BSplineBase();

    VectorXT operator()(double t);
    VectorXT evaluateDerivative(double t, int n);

    int degree() const;
    MatrixXT getControlPoints() const;
    VectorXT getKnotPoints() const;

 protected:
    virtual void computeBasisMatrices() = 0;

    double factorial(int i);

 private:

 protected:
    MatrixXT ctrl_pts_;
    VectorXT knot_pts_;
    std::vector<MatrixXT> M_;  // Vector of basis matrices
};

}  // namespace spline

#include "bspline_base.tpp"

#endif  // BSPLINES_BSPLINE_BASE_HPP_
