#ifndef POES_MATRIX_H
#define POES_MATRIX_H

#include "array.h"
#include <cmath>

namespace Tailor
{
    template<int nrow, int ncol>
        class Matrix
        {
            public:

                Matrix()
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            (*this)(r,c) = 0.;
                        }
                    }
                }

                Matrix transpose()
                {
                    auto m = *this;

                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            if (c == r) {
                                continue;
                            }

                            m(r,c) = (*this)(c,r);
                        }
                    }

                    return m;
                }

                double max() const
                {
                    double mx = TAILOR_BIG_NEG_NUM;

                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            mx = std::max(mx, (*this)(r,c));
                        }
                    }

                    return mx;
                }

                double& operator()(int row, int col)
                {
                    return data_[row * ncol + col]; 
                }

                double operator()(int row, int col) const
                {
                    return data_[row * ncol + col]; 
                }

                bool isnan() const
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            if (std::isnan((*this)(r,c)))
                            {
                                return true;
                            }
                        }
                    }

                    return false;
                }

                template<int nrowo, int ncolo> Matrix<nrow, ncolo> operator*(const Matrix<nrowo, ncolo>& other) const
                {
                    // matrix-matrix product

                    assert(nrowo == ncol);

                    Matrix<nrow, ncolo> mat;

                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncolo; ++c)
                        {
                            mat(r, c) = 0.; // take a row

                            for (int j=0; j<ncol; ++j)
                            {
                                //mat(r, c) += data_(r * ncol + j) * other(j * ncolo + c);
                                mat(r,c) += (*this)(r,j) * other(j,c);
                            }
                        }
                    }

                    return mat;
                }

                Matrix& operator-()
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            //data_[r * ncol + c] *= -1.;
                            (*this)(r,c) *= -1.;
                        }
                    }

                    return *this;
                }

                /*Matrix operator*(const Vector<ncol>& v) const
                {
                    // matrix-diag matrix product
                    // v is ncol x ncol diagonal square matrix

                    Matrix<nrow, ncol> mat;

                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            //mat(r, c) = data_(r * ncol + c) * v(c);
                            mat(r,c) = (*this)(r,c) * v(c);
                        }
                    }

                    return mat;
                }*/

                Vector<nrow> operator*(const Vector<ncol>& v) const
                {
                    // matrix-vector product

                    Vector<nrow> vec;
                    vec = 0.;

                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            //vec[r] += data_(r * ncol + c) * v(c);
                            vec[r] += (*this)(r,c) * v[c];
                        }
                    }

                    return vec;
                }

                Matrix operator*(double s) const
                {
                    Matrix<nrow, ncol> mat;

                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            //mat(r, c) = s * data_(r * ncol + c);
                            mat(r,c) = s * (*this)(r,c);
                        }
                    }

                    return mat;
                }

                Matrix& operator+=(const Matrix<nrow, ncol>& M)
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            //data_[r * ncol + c] += M(r, c);
                            (*this)(r,c) += M(r, c);
                        }
                    }

                    return *this;
                }

                Matrix& operator-=(const Matrix<nrow, ncol>& M)
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            (*this)(r,c) -= M(r, c);
                        }
                    }

                    return *this;
                }

                Matrix& add_diag(double d)
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        (*this)(r,r) += d;
                    }

                    return *this;
                }

                Matrix& operator+=(double d)
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            //data_[r * ncol + c] += d;
                            (*this)(r,c) += d;
                        }
                    }

                    return *this;
                }

                Matrix& operator*=(double d)
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            (*this)(r,c) *= d;
                        }
                    }

                    return *this;
                }

                Matrix operator+(const Matrix<nrow, ncol>& M) const
                {
                    Matrix<nrow, ncol> mat;

                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            //mat(r, c) = data_(r * ncol + c) + M(r, c);
                            mat(r,c) = (*this)(r,c) + M(r,c);
                        }
                    }

                    return mat;
                }

                Matrix operator-(const Matrix<nrow, ncol>& M) const
                {
                    Matrix<nrow, ncol> mat;

                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            //mat(r, c) = data_(r * ncol + c) - M(r, c);
                            mat(r,c) = (*this)(r,c) - M(r,c);
                        }
                    }

                    return mat;
                }

                Matrix& operator=(double s)
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            //data_[r * ncol + c] = s;
                            (*this)(r,c) = s;
                        }
                    }

                    return *this;
                }

                bool operator==(const Matrix<nrow, ncol>& M) const
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            if ((*this)(r,c) != M(r,c))
                            {
                                return false;
                            }
                        }
                    }

                    return true;
                }


            private:

                Vector<nrow * ncol> data_;
        };

    vararray operator/(const vararray& f, const Matrix<NVAR, NVAR>& A);

    template<int nrow, int ncol>
        std::ostream& operator<<(std::ostream& os, const Matrix<nrow, ncol>& mat)
        {
            for (int i=0; i<nrow; ++i)
            {
                os << "MAT(" << i << ",:): ";
                for (int j=0; j<ncol; ++j)
                {
                    os << mat(i,j) << ", "; 
                }
                os << "\n"; 
            }

            return os;
        }

    using varmat = Matrix<NVAR, NVAR>;
}

#endif
