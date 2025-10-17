#ifndef VERTEX_H
#define VERTEX_H

#include <point.h>
#include <face.h>
#include <typedef.h>

template<typename T, unsigned int DIMENTION>
class Vertex
{
public:
    Vertex() = default;
    Point<T, DIMENTION> & getPosition()
    {
        return position;
    }
    const Point<T, DIMENTION> & getPosition() const
    {
        return position;
    }
    faceIndex & faceNeighboor()
    {
        return _faceNeighboor;
    }

private:
    /**
     * @brief Position of the vertex
     * 
     */
    Point<T, DIMENTION> position;

    /**
     * @brief One of the triangles this vertex is part of.
     * 
     */
    faceIndex _faceNeighboor;
};

#endif //VERTEX_H