#ifndef BSPLINES_BSPLINE3_HPP_
#define BSPLINES_BSPLINE3_HPP_

#include "bspline_base.hpp"

namespace spline {

template<typename T>
class BSpline3 : public BSplineBase<T, 3> {
 public:
    BSpline3() : BSplineBase<T, 3>() {}
    explicit BSpline3(const Eigen::MatrixXd &ctrl_pts) :
                                                BSplineBase<T, 3>(ctrl_pts) {}
    ~BSpline3();

 protected:
};

}  // namespace spline

#include "bspline3.tpp"
#endif  //BSPLINES_BSPLINE3_HPP_
