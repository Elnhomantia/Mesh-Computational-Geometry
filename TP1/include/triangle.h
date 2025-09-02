#ifndef TRIANGLE_H
#define TRIANGLE_H

#include<vertex.h>
#include<cassert>

template<typename T, unsigned int DIMENTION, unsigned int EDGES = 3>
requires (EDGES > 2)
class Triangle
{
public:
    Vertex<T, DIMENTION> & vertexAt(unsigned int index)
    {
        assert(index < EDGES);
        return vertices[index];
    }
    int & triangleNeighboorIDAt(unsigned int index)
    {
        assert(index < EDGES);
        return triangleNeighboorID[index];
    }

private:
    /**
     * @brief EDGES vertices making a polygon.
     */
    Vertex<T, DIMENTION> vertices[EDGES];
    /**
     * @brief Neighboors polygon. It share an edge opposing the vertex of the same index.
     * For example, an ABC triangle has 3 vertices and 3 edges :
     * vertices[0] = A, vertices[1] = B and vertices[2] = C.
     * The neighboor triangle triangleNeighboorID[0] is sharing the BC edge.
     */
    int triangleNeighboorID[EDGES];
};

#endif //TRIANGLE_H