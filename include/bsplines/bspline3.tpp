
namespace spline {


template<typename T>
BSpline3<T>::~BSpline3() {}

template<typename T>
void BSpline3<T>::computeBasisMatrices() {
    for (int i{0}; i != this->knot_pts_.rows(); ++i)
        this->M_.push_back(computeNonUniformMatrix(i));
}

template<typename T>
bool BSpline3<T>::updateBasisMatrices() {
    int r{this->knot_pts_.rows()};
    int i = std::max(r-7, 0);
    for (; i != r; ++i)
        this->M_[i] = computeNonUniformMatrix(i);

    return true;
}

template<typename T>
typename BSplineBase<T,3>::MatrixXT BSpline3<T>::computeNonUniformMatrix(int i) {
    int s{this->knot_pts_.rows()};
    typename BSplineBase<T,3>::MatrixXT M =
                                    BSplineBase<T,3>::MatrixXT::Zero(4, 4);
    if (i+1 >= s || i+2 >= s || i+3 >= s || i-1 < 0 || i-2 < 0)
        return M;

    T ti{this->knot_pts_(i)}, tim1{this->knot_pts_(i-1)};
    T tim2{this->knot_pts_(i-2)}, tip1{this->knot_pts_(i+1)};
    T tip2{this->knot_pts_(i+2)}, tip3{this->knot_pts_(i+3)};

    T num1{tip1 - ti}, num2{ti - tim1};
    T den1{tip1 - tim1}, den2{tip1 - tim2}, den3{tip2 - tim1};
    T den5{tip3-ti}, den4{tip2-ti};
    // Watch out for integer division
    if (den1 != 0) {
        M(0, 0) = den2 != 0 ? pow(num1, 2)/(den1*den2) : 0;
        M(0, 2) = den3 != 0 ? pow(num2, 2)/(den1*den3) : 0;
        M(1, 2) = den3 != 0 ? 3*num1*num2/(den1*den3) : 0;
        M(2, 2) = den3 != 0 ? 3*pow(num1, 2)/(den1*den3) : 0;
        M(3, 2) = den3 != 0 && den5 != 0 ? -pow(num1, 2)/(den4*den3) : 0;
    }
    if (den4 != 0)
        M(3, 3) = den5 != 0 ? pow(num1, 2)/(den4*den5) : 0;

    M(0, 1) = 1 - M(0, 0) - M(0, 2);
    M(1, 0) = -3*M(0, 0);
    M(1, 1) = 3*M(0, 0) - M(1, 2);
    M(2, 0) = 3*M(0, 0);
    M(2, 1) = -3*M(0, 0) - M(2, 2);
    M(3, 0) = -M(0, 0);
    M(3, 2) -= (M(2, 2)/3.0 + M(3, 3));
    M(3, 1) = M(0, 0) - M(3, 2) - M(3, 3);

    return M;
}

}
