#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>

#include <vertex.h>
#include <face.h>


class Mesh
{
public:
    void loadOffFile(const std::string & file);
    bool check();

    inline std::vector<Vertex<float, 3>>& getVertices() { return vertices; }
    inline std::vector<Face<float, 3, 3>>& getFaces() { return faces; }

    std::vector<Face<float, 3, 3>> getFaceNeighboors(vertexIndex index);
    std::vector<Vertex<float, 3>> getVertexNeighboors(vertexIndex index);

private:
    std::vector<Vertex<float, 3>> vertices;
    std::vector<Face<float, 3, 3>> faces;
};

#endif //MESH_H