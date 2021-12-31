#include "rotationmatrix.h"

namespace Tailor
{
    Matrix3 set_r(double theta, const Vector3& u)
    {
        // https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
        // section: Rotation matrix from axis and angle

        Matrix3 m;

        double ct = std::cos(theta);
        double st = std::sin(theta);
        double mct = 1. - ct;

        m(0,0) = ct + u(0) * u(0) * mct;
        m(0,1) = u(0) * u(1) * mct - u(2) * st;
        m(0,2) = u(0) * u(2) * mct + u(1) * st;

        m(1,0) = u(1) * u(0) * mct + u(2) * st;
        m(1,1) = ct + u(1) * u(1) * mct;
        m(1,2) = u(1) * u(2) * mct - u(0) * st;

        m(2,0) = u(2) * u(0) * mct - u(1) * st;
        m(2,1) = u(2) * u(1) * mct + u(0) * st;
        m(2,2) = ct + u(2) * u(2) * mct;

        return m;
    }

    //Matrix3 RotationMatrix::set_rx(double yaw)
    //{
    //    Matrix3 m;

    //    m(0,0) = 1.;
    //    m(0,1) = 0.;
    //    m(0,2) = 0.;

    //    m(1,0) = 0.;
    //    m(1,1) = std::cos(yaw);
    //    m(1,2) = -std::sin(yaw);

    //    m(2,0) = 0.;
    //    m(2,1) = std::sin(yaw);
    //    m(2,2) = std::cos(yaw);

    //    return m;
    //}

    //Matrix3 RotationMatrix::set_ry(double pitch)
    //{
    //    Matrix3 m;

    //    m(0,0) = std::cos(pitch);
    //    m(0,1) = 0.;
    //    m(0,2) = std::sin(pitch);
    //          
    //    m(1,0) = 0.;
    //    m(1,1) = 1.;
    //    m(1,2) = 0.;
    //          
    //    m(2,0) = -std::sin(pitch);
    //    m(2,1) = 0.;
    //    m(2,2) = std::cos(pitch);

    //    return m;
    //}

    //Matrix3 RotationMatrix::set_rz(double roll)
    //{
    //    Matrix3 m;

    //    m(0,0) = std::cos(roll);
    //    m(0,1) = -std::sin(roll);
    //    m(0,2) = 0.;
    //          
    //    m(1,0) = std::sin(roll);
    //    m(1,1) = std::cos(roll);
    //    m(1,2) = 0.;
    //          
    //    m(2,0) = 0.;
    //    m(2,1) = 0.;
    //    m(2,2) = 1.;

    //    return m;
    //}

    //void RotationMatrix::rotate(double angle, int axis, const Vector3& pivot, Vector3& v)
    void RotationMatrix::rotate(double angle, const Vector3& axis, const Vector3& pivot, Vector3& v)
    {
        v = v - pivot;

        auto m = set_r(angle, axis);

        v = m * v;
        v = v + pivot;

        assert(!std::isnan(v(0)));
        assert(!std::isnan(v(1)));
        assert(!std::isnan(v(2)));
    }

    //Vector3 RotationMatrix::mul(const Matrix3& m, const Vector3& v)
    //{
    //    Vector3 c = m * v;

    //    assert(!c.isnan());

    //    return c;
    //}

    //Vector3 RotationMatrix::mul(const trimat& m, const Vector3& v)
    //{
    //    Vector3 c;

    //    c(0) = m[0] * v(0) + m[1] * v(1) + m[2] * v(2);
    //    c(1) = m[3] * v(0) + m[4] * v(1) + m[5] * v(2);
    //    c(2) = m[6] * v(0) + m[7] * v(1) + m[8] * v(2);

    //    assert(!c.isnan());

    //    return c;
    //}

    double angle(const Vector3& a, const Vector3& b)
    {
        double l = dot(a, b) / (len(a) * len(b));

        return std::acos(l);
    }
}
