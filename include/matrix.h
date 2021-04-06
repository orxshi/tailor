#ifndef POES_MATRIX_H
#define POES_MATRIX_H

//#include "array.h"
#include "core.h"
#include <boost/serialization/utility.hpp>
#include <cmath>
#include <array>
#include <cassert>
#include <iostream>

namespace Tailor
{
    template<int nrow, int ncol, class T>
        class Matrix
        {
            public:

                Matrix();
                Matrix(T a, T b, T c);
                Matrix(T s);

                Matrix transpose();
                T& operator()(int row, int col);
                T operator()(int row, int col) const;
                T& operator()(int row);
                T operator()(int row) const;
                bool isnan() const;
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
                                    mat(r,c) += (*this)(r,j) * other(j,c);
                                }
                            }
                        }

                        return mat;
                    }
                Matrix& operator-();
                Matrix operator*(T s) const;
                Matrix operator/(T s) const;
                Matrix& operator+=(const Matrix<nrow, ncol, T>& M);
                Matrix& operator-=(const Matrix<nrow, ncol, T>& M);
                Matrix& add_diag(T d);
                Matrix& operator+=(T d);
                Matrix& operator*=(T d);
                Matrix& operator/=(T s);
                Matrix operator+(const Matrix<nrow, ncol, T>& M) const;
                Matrix operator-(const Matrix<nrow, ncol, T>& M) const;
                Matrix& operator=(T s);
                bool operator==(const Matrix<nrow, ncol, T>& M) const;
                double len();
                int nelm() const;

                template<class Archive> void serialize(Archive & ar, const unsigned int version)
                {
                    ar & data_;
                }

                friend class boost::serialization::access;

            private:

                std::array<T, nrow * ncol> data_;
        };

    using Matrix5 = Matrix<5, 5, double>;

    template<int nrow, class T>
        using Vector = Matrix<nrow, 1, T>;

    using Vector5 = Vector<5, double>;
    using Vector3 = Vector<3, double>;
    using Vector3Int = Vector<3, int>;

    template<int nrow, int ncol, class T>
        std::ostream& operator<<(std::ostream& os, const Matrix<nrow, ncol, T>& mat);

    template<int nrow, class T>
        double len(const Vector<nrow, T> a);

    template<int nrow, class T>
        Vector<nrow, T> cross(const Vector<nrow, T>& a, const Vector<nrow, T>& b);

    template<int nrow, class T>
        double dot(const Vector<nrow, T>& a, const Vector<nrow, T>& b);

    template<int nrow, class T>
        Vector<nrow, T> normals(const Vector<nrow, T>& a, const Vector<nrow, T>& b, const Vector<nrow, T>& c);

    template<int nrow, class T>
        Vector<nrow, T> normalize(const Vector<nrow, T>& a);

    Vector5 operator/(const Vector5& f, const Matrix5& A);

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> operator*(T s, const Matrix<nrow, ncol, T>& m);

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> abs(const Matrix<nrow, ncol, T>& m);

    template<int nrow, int ncol, class T>
        T max(const Matrix<nrow, ncol, T>& m);
    
    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> unit_matrix();
}

#include "matrix.hpp"

#endif
