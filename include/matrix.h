#ifndef POES_MATRIX_H
#define POES_MATRIX_H

#include "array.h"
#include <cmath>

namespace Tailor
{
    template<int nrow, int ncol, class T>
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

                Matrix(T a, T b, T c)
                {
                    assert(nrow == 3);
                    assert(ncol == 1);

                    (*this)(0) = a;
                    (*this)(1) = b;
                    (*this)(2) = c;
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

                T max() const
                {
                    T mx = TAILOR_BIG_NEG_NUM;

                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            mx = std::max(mx, (*this)(r,c));
                        }
                    }

                    return mx;
                }

                T& operator()(int row, int col)
                {
                    return data_[row * ncol + col]; 
                }

                T operator()(int row, int col) const
                {
                    return data_[row * ncol + col]; 
                }

                T& operator()(int row)
                {
                    assert(ncol == 1);
                    int col = 1;
                    return (*this)(row, col);
                }

                T operator()(int row) const
                {
                    assert(ncol == 1);
                    int col = 1;
                    return (*this)(row, col);
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

                template<int nrowo, int ncolo>
                    Matrix<nrow, ncolo, T> operator*(const Matrix<nrowo, ncolo, T>& other) const
                {
                    // matrix-matrix product

                    assert(nrowo == ncol);

                    Matrix<nrow, ncolo, T> mat;

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

                //Vector<nrow> operator*(const Vector<ncol>& v) const
                //{
                //    // matrix-vector product

                //    Vector<nrow> vec;
                //    vec = 0.;

                //    for (int r=0; r<nrow; ++r)
                //    {
                //        for (int c=0; c<ncol; ++c)
                //        {
                //            //vec[r] += data_(r * ncol + c) * v(c);
                //            vec[r] += (*this)(r,c) * v[c];
                //        }
                //    }

                //    return vec;
                //}

                Matrix operator*(T s) const
                {
                    Matrix<nrow, ncol, T> mat;

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

                Matrix operator/(T s) const
                {
                    Matrix<nrow, ncol, T> mat;

                    for (int r=0; r<nrow; ++r)
                    {
                        for (int c=0; c<ncol; ++c)
                        {
                            //mat(r, c) = s * data_(r * ncol + c);
                            mat(r,c) = (*this)(r,c) / s;
                        }
                    }

                    return mat;
                }

                Matrix& operator+=(const Matrix<nrow, ncol, T>& M)
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

                Matrix& operator-=(const Matrix<nrow, ncol, T>& M)
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

                Matrix& add_diag(T d)
                {
                    for (int r=0; r<nrow; ++r)
                    {
                        (*this)(r,r) += d;
                    }

                    return *this;
                }

                Matrix& operator+=(T d)
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

                Matrix& operator*=(T d)
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

                Matrix operator+(const Matrix<nrow, ncol, T>& M) const
                {
                    Matrix<nrow, ncol, T> mat;

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

                Matrix operator-(const Matrix<nrow, ncol, T>& M) const
                {
                    Matrix<nrow, ncol, T> mat;

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

                Matrix& operator=(T s)
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

                bool operator==(const Matrix<nrow, ncol, T>& M) const
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

                double len()
                {
                    return std::sqrt(std::pow(data_[0], 2.) + std::pow(data_[1], 2.) + std::pow(data_[2], 2.));
                }

                template<class Archive> void serialize(Archive & ar, const unsigned int version)
                {
                    ar & data_;
                }

                friend class boost::serialization::access;

            private:

                Array<T, nrow * ncol> data_;
        };

    template<int nrow, class T>
        using Vector = Matrix<nrow, 1, T>;

    using Vector5 = Vector<5, double>;
    using Vector3 = Vector<3, double>;
    using Vector3Int = Vector<3, int>;

    template<int nrow, int ncol, class T>
        std::ostream& operator<<(std::ostream& os, const Matrix<nrow, ncol, T>& mat)
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

    using Matrix5 = Matrix<5, 5, double>;

    template<int nrow, class T>
        double len(const Vector<nrow, T> a);

    template<int nrow, class T>
        Vector<nrow, T> cross(const Vector<nrow, T>& a, const Vector<nrow, T>& b);

    template<int nrow, class T>
        double dot(const Vector<nrow, T>& a, const Vector<nrow, T>& b);

    template<int nrow, class T>
        Vector<nrow, T> normals(const Vector<nrow, T>& a, const Vector<nrow, T>& b, const Vector<nrow, T>& c);
}

#endif
