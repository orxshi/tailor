#ifndef VEC2_HPP
#define	VEC2_HPP

namespace Tailor
{
    template<class T>
        vec3<T>::vec3(T x, T y, T z)
        {
            set(x, y, z);
        }

    template<class T>
        T vec3<T>::len() const
        {
            return std::sqrt(std::pow(data_[0], 2.) + std::pow(data_[1], 2.) + std::pow(data_[2], 2.));
        }

    template<class T>
        const T& vec3<T>::operator()(int i) const
        {
            assert(i >= 0);
            if (i >= 3)
            {
                std::cout << "i: " << i << std::endl;
            }
            assert(i < 3);
            return data_[i];
        }

    template<class T>
        void vec3<T>::set(T x, T y, T z)
        {
            data_[0] = x;
            data_[1] = y;
            data_[2] = z;
        }

    template<class T>
        void vec3<T>::set(T v)
        {
            data_[0] = v;
            data_[1] = v;
            data_[2] = v;
        }

    template<class T>
        void vec3<T>::set(int i, T v)
        {
            data_[i] = v;
        }

    template<class T>
        void vec3<T>::set_x(T v)
        {
            data_[0] = v;
        }

    template<class T>
        void vec3<T>::set_y(T v)
        {
            data_[1] = v;
        }

    template<class T>
        void vec3<T>::set_z(T v)
        {
            data_[2] = v;
        }

    template<class T>
        vec3<T> vec3<T>::operator-(const vec3& x) const
        {
            vec3 c;
            c.data_[0] = data_[0] - x(0);
            c.data_[1] = data_[1] - x(1);
            c.data_[2] = data_[2] - x(2);
            return c;
        }

    template<class T>
        vec3<T> vec3<T>::operator+(const vec3& x) const
        {
            vec3 c;
            c.data_[0] = data_[0] + x(0);
            c.data_[1] = data_[1] + x(1);
            c.data_[2] = data_[2] + x(2);
            return c;
        }

    template<class T>
        vec3<T>& vec3<T>::operator=(const vec3& other)
        {
            data_[0] = other.data_[0];
            data_[1] = other.data_[1];
            data_[2] = other.data_[2];

            assert(!std::isnan(other.data_[0]));
            assert(!std::isnan(other.data_[1]));
            assert(!std::isnan(other.data_[2]));

            if (data_[0] != other.data_[0])
            {
                std::cout << "data0: " << data_[0] << std::endl;
                std::cout << "odata0: " << other.data_[0] << std::endl;
            }
            assert(data_[0] == other.data_[0]);
            assert(data_[1] == other.data_[1]);
            assert(data_[2] == other.data_[2]);

            return *this;
        }

    template<class T>
        vec3<T>& vec3<T>::operator+=(const vec3& other)
        {
            data_[0] = data_[0] + other.data_[0];
            data_[1] = data_[0] + other.data_[1];
            data_[2] = data_[0] + other.data_[2];

            return *this;
        }

    template<class T>
        bool vec3<T>::operator==(const vec3& other) const
        {
            return (data_[0] == other(0) && data_[1] == other(1) && data_[2] == other(2));
        }

    template<class T>
        bool vec3<T>::operator!=(const vec3& other) const
        {
            return !(*this == other);
        }

    template<class T>
        vec3<T> vec3<T>::operator/(T x) const
        {
            vec3 c;
            c.data_[0] = data_[0] / x;
            c.data_[1] = data_[1] / x;
            c.data_[2] = data_[2] / x;
            return c;
        }

    template<class T>
        vec3<T> vec3<T>::operator*(T x) const
        {
            vec3 c;
            c.data_[0] = data_[0] * x;
            c.data_[1] = data_[1] * x;
            c.data_[2] = data_[2] * x;
            return c;
        }


    template<class T>
        vec3<T> normals(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c)
        {
            return cross((a-b), (c-b));
        }

    template<class T>
        vec3<T> cross(const vec3<T>& a, const vec3<T>& b)
    {
        vec3<T> c;

        c.set_x(a(1) * b(2) - a(2) * b(1));
        c.set_y(a(2) * b(0) - a(0) * b(2));
        c.set_z(a(0) * b(1) - a(1) * b(0));

        return c;
    }

    template<class T>
        T dotp(const vec3<T>& a, const vec3<T>& b)
    {
        T c = a(0) * b(0) + a(1) * b(1) + a(2) * b(2);

        return c;
    }

    template<class T>
        T angle(const vec3<T>& a, const vec3<T>& b)
    {
        return std::acos(dotp(a, b) / (len(a) * len(b))); // radian
    }

    template<class T>
        vec3<T> normalize(const vec3<T>& a)
    {
        return a / a.len();
    }
}

#endif
