#ifndef POINT_H
#define POINT_H

#include <concepts>
#include <type_traits>
#include <cassert>

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
    T& operator[](Axis axis)
    {
        assert((int)axis < DIMENTION);
        return this->coordinates[axis];
    }
        T& operator[](unsigned int axis)
    {
        assert(axis < DIMENTION);
        return this->coordinates[axis];
    }

private:
    T coordinates[DIMENTION];
};


#endif //POINT_H