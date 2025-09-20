#ifndef FACE_H
#define FACE_H

#include <vertex.h>
#include <cassert>
#include <typedef.h>

template<typename T, unsigned int DIMENTION, unsigned int VERTICE_NUMBER = 3>
requires (VERTICE_NUMBER > 2)
class Face
{
public:
    vertexIndex & vertexAt(unsigned int index)
    {
        assert((unsigned int)index < VERTICE_NUMBER);
        return vertices[index];
    }
    faceIndex & facesNeighboorIDAt(unsigned int index)
    {
        assert(index < VERTICE_NUMBER);
        return facesNeighboorID[index];
    }

private:
    /**
     * @brief VERTICE_NUMBER vertices making a polygon.
     */
    vertexIndex vertices[VERTICE_NUMBER];
    /**
     * @brief Neighboors polygon. It share an edge opposing the vertex of the same index.
     * For example, an ABC triangle has 3 vertices and 3 edges :
     * vertices[0] = A, vertices[1] = B and vertices[2] = C.
     * The neighboor triangle facesNeighboorID[0] is sharing the BC edge.
     */
    faceIndex facesNeighboorID[VERTICE_NUMBER];
};

#endif //FACE_H