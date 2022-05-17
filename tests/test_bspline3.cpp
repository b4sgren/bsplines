#include <gtest/gtest.h>
#include "bspline3.hpp"
#include "bspline_base.hpp"

#include <eigen3/Eigen/Dense>
#include <iostream>

class Spline3_Fixture : public ::testing::Test {
 public:
    Spline3_Fixture() {
        ctrl_pts_ = Eigen::MatrixXd(8, 3);
        ctrl_pts_ << 0, 0, 0,
                    20, 10, 10,
                    15, 30, 30,
                    20, 35, 35,
                    22, 36, 40,
                    25, 37, 25,
                    40, 40, 10,
                    40, 40, 0;

        spline_ = new spline::BSpline3<double>(ctrl_pts_);
    }

    ~Spline3_Fixture() {
        // delete spline_;
    }

    Eigen::MatrixXd ctrl_pts_;
    spline::BSplineBase<double, 3> *spline_;
};

TEST_F(Spline3_Fixture, TestInitialization) {
    EXPECT_EQ(spline_->degree(), 3);
    EXPECT_EQ(spline_->getControlPoints().rows(), 8);
    EXPECT_EQ(spline_->getKnotPoints().rows(),
              ctrl_pts_.rows() + spline_->degree()+1);
}

TEST_F(Spline3_Fixture, TestSplineEvaluation) {
    double t = 2.5;
    Eigen::Vector3d x = spline_->operator()(t);

    double x_true{20.9583}, y_true{35.4167}, z_true{37.0833};

    EXPECT_NEAR(x(0), x_true, 1e-3);
    EXPECT_NEAR(x(1), y_true, 1e-3);
    EXPECT_NEAR(x(2), z_true, 1e-3);

    t = 0.0;
    x_true = 0.0;
    y_true = 0.0;
    z_true = 0.0;

    x = spline_->operator()(t);
    EXPECT_NEAR(x(0), x_true, 1e-7);
    EXPECT_NEAR(x(1), y_true, 1e-7);
    EXPECT_NEAR(x(2), z_true, 1e-7);

    t = 5.0;
    x_true = 40.0;
    y_true = 40.0;
    z_true = 0.0;

    x = spline_->operator()(t);
    EXPECT_NEAR(x(0), x_true, 1e-7);
    EXPECT_NEAR(x(1), y_true, 1e-7);
    EXPECT_NEAR(x(2), z_true, 1e-7);

}
