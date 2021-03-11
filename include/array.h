#ifndef ARRAY_H
#define	ARRAY_H

#include <array>
#include <vector>
#include <iostream>
#include <cassert>
#include "core.h"
#include <boost/serialization/utility.hpp>

namespace Tailor
{
    template<typename T, size_t max_size>
        class Array
        {
            size_t size_;
            std::array<T, max_size> data_;

            public:

            Array(): size_(0) {}
            Array(bool init_to_max)
            {
                if (init_to_max) {
                    size_ = max_size;
                }
                else {
                    size_ = 0;
                }
            }

            void set_data(T v)
            {
                for (int i=0; i<size_; ++i) {
                    data_[i] = v;
                }
            }

            Array& operator=(const std::vector<T>& v)
            {
                assert(!v.empty());
                if (v.size() > max_size)
                {
                    std::cout << "vsize: " << v.size() << std::endl;
                    std::cout << "max_size: " << max_size << std::endl;
                }
                assert(v.size() <= max_size);
                size_ = v.size();
                for (int i=0; i<v.size(); ++i) {
                    data_[i] = v[i];
                }
                //std::copy(v.begin(), v.end(), data_.data());
                return *this;
            }

            Array(const std::vector<T>& v)
            {
                assert(!v.empty());
                assert(v.size() <= max_size);
                size_ = v.size();
                std::copy(v.begin(), v.end(), data_.data());
            }

            size_t mem() const
            {
                return (max_size * sizeof(T) + sizeof(size_t));
            }

            /*const T& operator()(int i) const
            {
                assert(i < size_);
                return data_[i];
            }*/

            T& operator[](int i)
            {
                if (i >= size_)
                {
                    std::cout << "i: " << i << std::endl;
                    std::cout << "size_: " << size_ << std::endl;
                }
                assert(i < size_);
                return data_[i];
            }
            
            const T& operator[](int i) const
            {
                if (i >= size_)
                {
                    std::cout << "i: " << i << std::endl;
                    std::cout << "size_: " << size_ << std::endl;
                }
                assert(i < size_);
                return data_[i];
            }

            void reserve(int i)
            {
                assert(i >= 0);
                assert(i <= max_size);
                size_ = i;
            }

            T& back()
            {
                assert(size_ != 0);
                return data_[size_ - 1];
            }

            const T& back() const
            {
                assert(size_ != 0);
                return data_[size_ - 1];
            }

            void add(const T& type)
            {
                assert(size_ != max_size);
                data_[size_] = type;
                ++size_;
            }

            void insert(const Array& other)
            {
                assert((other.size() + size_) <= max_size);
                for (int i=size_; i<size_+other.size(); ++i)
                {
                    data_[i] = other[i-size_];
                }
            }

            void push_back(const T& type)
            {
                add(type);
            }

            void clear()
            {
                size_ = 0;
            }

            bool empty() const
            {
                if (size_ == 0) {
                    return true;
                }

                return false;
            }

            const T* pbegin() const
            {
                return data_.data();
            }

            const T* pend() const
            {
                return data_.data() + size_;
            }

            typename std::array<T, max_size>::iterator begin()
            {
                return data_.begin();
            }

            typename std::array<T, max_size>::iterator end()
            {
                return (data_.begin() + size_);
            }

            typename std::array<T, max_size>::const_iterator begin() const
            {
                return data_.begin();
            }

            typename std::array<T, max_size>::const_iterator end() const
            {
                return (data_.begin() + size_);
            }

            size_t size() const
            {
                return size_;
            }

            void erase(typename std::array<T, max_size>::iterator it)
            {
                for (auto j=it; j!=end() - 1; ++j)
                {
                    *j = *std::next(j);
                }
                --size_;
            }

            void erase(const T& type)
            {
                for (int i=0; i<size_; ++i)
                {
                    if (data_[i] == type)
                    {
                        for (int j=i; j<size_-1; ++j)
                        {
                            data_[j] = data_[j+1];
                        }
                        --size_;
                        break;
                    }
                }
            }

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & size_;
                ar & data_;
            }

            friend class boost::serialization::access;
        };

    //template<size_t max_size>
    //    class Vector: public Array<double, max_size>
    //    {
    //        public:

    //            Vector(): Array<double, max_size>(true) {}

    //            Vector& operator=(double v)
    //            {
    //                this->set_data(v);
    //                return *this;
    //            }

    //        private:
    //    };

    //template<size_t max_size> using Vector = Array<double, max_size>;
    //using vararray = Vector<NVAR>;

    //template<size_t max_size> 
    //    Vector<max_size> operator+(const Vector<max_size>& a, const Vector<max_size>& b)
    //    {
    //        Vector<max_size> c;
    //        for (int i=0; i<max_size; ++i)
    //        {
    //            c[i] = a[i] + b[i];
    //        }

    //        return c;
    //    }
    //template<size_t max_size> 
    //    Vector<max_size> operator-(const Vector<max_size>& a, const Vector<max_size>& b)
    //    {
    //        Vector<max_size> c;
    //        for (int i=0; i<max_size; ++i)
    //        {
    //            c[i] = a[i] - b[i];
    //        }

    //        return c;
    //    }

    //template<size_t max_size> 
    //    Vector<max_size> operator*(const Vector<max_size>& a, double d)
    //    {
    //        Vector<max_size> v;
    //        for (int i=0; i<max_size; ++i)
    //        {
    //            v[i] = a[i] * d;
    //        }

    //        return v;
    //    }

    //template<size_t max_size> 
    //    Vector<max_size> operator*(double d, const Vector<max_size>& a)
    //    {
    //        Vector<max_size> v;
    //        for (int i=0; i<max_size; ++i)
    //        {
    //            v[i] = a[i] * d;
    //        }

    //        return v;
    //    }

    //template<size_t max_size> 
    //    Vector<max_size> operator/(const Vector<max_size>& a, double d)
    //    {
    //        Vector<max_size> v;
    //        for (int i=0; i<max_size; ++i)
    //        {
    //            v[i] = a[i] / d;
    //        }

    //        return v;
    //    }
    //template<size_t max_size> 
    //    Vector<max_size>& operator+=(Vector<max_size>& a, const Vector<max_size>& b)
    //    {
    //        for (int i=0; i<max_size; ++i)
    //        {
    //            a[i] += b[i];
    //        }

    //        return a;
    //    }
    // template<size_t max_size> 
    //    Vector<max_size>& operator-=(Vector<max_size>& a, const Vector<max_size>& b)
    //    {
    //        for (int i=0; i<max_size; ++i)
    //        {
    //            a[i] -= b[i];
    //        }

    //        return a;
    //    }

    // template<size_t max_size> 
    //    double max(Vector<max_size>& a)
    //    {
    //        double d = TAILOR_BIG_NEG_NUM;
    //        for (int i=0; i<max_size; ++i)
    //        {
    //            d = std::max(d, a[i]);
    //        }

    //        return d;
    //    }
}

#endif
