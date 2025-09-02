#ifndef VERTEX_H
#define VERTEX_H

#include<point.h>

template<typename T, unsigned int DIMENTION>
class Vertex
{
public:
    Point<T, DIMENTION> & getPosition()
    {
        return position;
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
    int triangleNeighboor;
};

#endif //VERTEX_H