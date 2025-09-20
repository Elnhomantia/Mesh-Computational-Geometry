#include <vector>
#include <string>
#include <face.h>
#include <vertex.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <map>

template <typename T, unsigned int DIMENTION, unsigned int VERTICE_NUMBER>
void parseOFF(const std::string &filename,
              std::vector<Vertex<T, DIMENTION>> &vertices,
              std::vector<Face<T, DIMENTION, VERTICE_NUMBER>> &faces)
{
    std::ifstream file(filename);
    if (!file)
        throw std::runtime_error("File won't open");

    std::string line;
    std::getline(file, line);
    if (line != "OFF")
        throw std::runtime_error("The file is not an OFF format");

    unsigned int numVertices = 0, numFaces = 0;
    std::getline(file, line);
    std::istringstream header(line);
    header >> numVertices >> numFaces;

    if (!header)
        throw std::runtime_error("Error when reading counts");

    vertices.resize(numVertices);
    for (unsigned int i = 0; i < numVertices; ++i)
    {
        std::getline(file, line);
        std::istringstream ss(line);
        for (unsigned int j = 0; j < DIMENTION; ++j)
        {
            if (!(ss >> vertices[i].getPosition()[j]))
                throw std::runtime_error("Error when reading vertices coordinates");
        }
    }

    //vertex index * 2 -> global face index + local face index
    std::map<std::pair<vertexIndex, vertexIndex>, std::pair<faceIndex, unsigned int>> vertexNeighboor;

    faces.resize(numFaces);
    for (unsigned int i = 0; i < numFaces; ++i)
    {
        std::getline(file, line);
        std::istringstream ss(line);
        int edgesInFile;
        ss >> edgesInFile;

        if (edgesInFile != VERTICE_NUMBER)
            throw std::runtime_error("Error when reading faces, number of edge is wrong");

        for (unsigned int e = 0; e < VERTICE_NUMBER; ++e)
        {
            vertexIndex vertexId;
            if (!(ss >> vertexId))
                throw std::runtime_error("Error when reading faces coordinates");
            faces[i].vertexAt(e) = vertexId;
        }
        
        for (unsigned int e = 0; e < VERTICE_NUMBER; ++e)
        {
            vertexIndex v1 = faces[i].vertexAt(e);
            vertexIndex v2 = faces[i].vertexAt((e+1)%VERTICE_NUMBER);

            faceIndex minVertexIndex = std::min(v1, v2);
            faceIndex maxVertexIndex = std::max(v1, v2);

            auto found = vertexNeighboor.find(std::make_pair(minVertexIndex, maxVertexIndex));
            if(found == vertexNeighboor.end()) //not found
            {
                vertexNeighboor.insert(
                    std::make_pair(
                    std::make_pair(minVertexIndex, maxVertexIndex),
                    std::make_pair(i, (e+2)%VERTICE_NUMBER)));
            }
            else
            {
                std::pair<faceIndex, unsigned int> foundPair = found->second;
                faces[i].facesNeighboorIDAt((e+2)%VERTICE_NUMBER) = foundPair.first;
                faces[foundPair.first].facesNeighboorIDAt(foundPair.second) = i;
            }
            vertices[faces[i].vertexAt(e)].faceNeighboor() = i;

        }
    }
}