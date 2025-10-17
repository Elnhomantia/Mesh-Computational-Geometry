#ifndef CIRCULATOR_ON_VERTICES_H
#define CIRCULATOR_ON_VERTICES_H

#include <mesh.h>

#define INVALID_VERTEX_INDEX static_cast<vertexIndex>(-1)
#define INVALID_FACE_INDEX static_cast<faceIndex>(-1)

class CirculatorOnVertices
{
public:
    CirculatorOnVertices(Mesh& m, vertexIndex v, faceIndex startFace);

    Vertex<float,3>& operator*();

    CirculatorOnVertices& operator++();

    inline faceIndex getFaceIndex() { return currentFace; }
    inline vertexIndex getVertexIndex() { return currentVertex; }

    bool isValid() const;

private:
    Mesh& mesh;
    vertexIndex vIndex;
    faceIndex currentFace;
    vertexIndex currentVertex;
    faceIndex startFace;
    bool firstIteration;

    void updateCurrentVertex();

    faceIndex nextIncidentFace();
};

#endif //CIRCULATOR_ON_VERTICES_H