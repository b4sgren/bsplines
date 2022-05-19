#ifndef BSPLINES_BSPLINE3_HPP_
#define BSPLINES_BSPLINE3_HPP_

/*
Code specifics for a 3rd order BSpline. 3rd order seems to be the most useful
for what I want. If another order is needed the class can be modeled after this
one. All the generalized code that does not depend on the spline order is in
the bspline_base class
*/

#include "bspline_base.hpp"

#include <cmath>

namespace spline {

template<typename T>
class BSpline3 : public BSplineBase<T, 3> {
 public:
    BSpline3() : BSplineBase<T, 3>() {}
    explicit BSpline3(const Eigen::MatrixXd &ctrl_pts) :
                                                BSplineBase<T, 3>(ctrl_pts) {
        computeBasisMatrices();
    }
    ~BSpline3();

 protected:
    void computeBasisMatrices() override;
    bool updateBasisMatrices() override;
    typename BSplineBase<T,3>::MatrixXT computeNonUniformMatrix(int i);
};

}  // namespace spline

#include "bspline3.tpp"
#endif  //BSPLINES_BSPLINE3_HPP_
