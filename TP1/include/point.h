#ifndef POINT_H
#define POINT_H

#include <concepts>
#include <type_traits>
#include <cassert>
#include <cmath>
#include <array>

template<typename T>
concept Coordinate_concept = requires(T a, T b) 
{
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
};

enum class Axis
{
    X = 0,
    Y = 1,
    Z = 2
};

template<typename T, unsigned int DIMENTION>
requires Coordinate_concept<T> && (DIMENTION > 0)
class Point
{
public:
    Point() {
        coordinates.fill(T(0));
    }
    template<typename... Args>
    requires (sizeof...(Args) == DIMENTION) && ((std::convertible_to<Args,T> && ...))
    Point(Args... args) : coordinates{T(args)...} {}

    T& operator[](Axis axis)
    {
        assert((int)axis < DIMENTION);
        return this->coordinates[axis];
    }
    const T& operator[](Axis axis) const 
    {
        assert((int)axis < DIMENTION);
        return coordinates[axis];
    }
    T& operator[](unsigned int axis)
    {
        assert(axis < DIMENTION);
        return this->coordinates[axis];
    }
    const T& operator[](unsigned int axis) const
    {
        assert(axis < DIMENTION);
        return this->coordinates[axis];
    }


    
    Point<T,DIMENTION>& operator+=(const Point<T,DIMENTION>& rhs) 
    {
        for (unsigned int i=0; i<DIMENTION; ++i) coordinates[i] += rhs.coordinates[i];
        return *this;
    }

    Point<T,DIMENTION>& operator-=(const Point<T,DIMENTION>& rhs) 
    {
        for (unsigned int i=0; i<DIMENTION; ++i) coordinates[i] -= rhs.coordinates[i];
        return *this;
    }

    Point<T,DIMENTION>& operator*=(T scalar) 
    {
        for (unsigned int i=0; i<DIMENTION; ++i) coordinates[i] *= scalar;
        return *this;
    }

    Point<T,DIMENTION>& operator/=(T scalar) 
    {
        for (unsigned int i=0; i<DIMENTION; ++i) coordinates[i] /= scalar;
        return *this;
    }

    T norm() const 
    {
        T sum = T(0);
        for(unsigned int i=0; i<DIMENTION; ++i) sum += coordinates[i]*coordinates[i];
        return std::sqrt(sum);
    }

    Point<T,DIMENTION> operator+(const Point<T,DIMENTION>& rhs) 
    {
        *this += rhs;
        return *this;
    }

    Point<T,DIMENTION> operator-(const Point<T,DIMENTION>& rhs) 
    {
        *this -= rhs;
        return *this;
    }

    Point<T,DIMENTION> operator*(T scalar) 
    {
        *this *= scalar;
        return *this;
    }

    Point<T,DIMENTION> operator/(T scalar) 
    {
        *this /= scalar;
        return *this;
    }

private:
    std::array<T, DIMENTION> coordinates;
};


#endif //POINT_H