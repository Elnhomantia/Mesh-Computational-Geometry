#include <mesh.h>
#include "OFFparser.cpp"
#include "meshChecker.cpp"
#include <set>

void Mesh::loadOffFile(const std::string & file)
{
    parseOFF<float, 3, 3>(file, vertices, faces);
}

bool Mesh::check()
{
    return checkMesh(faces);
}

std::vector<faceIndex> Mesh::getFaceNeighboors(vertexIndex index)
{
    faceIndex f = vertices[index].faceNeighboor();
    
    std::set<faceIndex> s;
    return std::vector<faceIndex>(s.begin(), s.end());
}
std::vector<Vertex<float, 3>> Mesh::getVertexNeighboors(vertexIndex index)
{
    faceIndex f = vertices[index].faceNeighboor();
    std::set<vertexIndex> s;
}