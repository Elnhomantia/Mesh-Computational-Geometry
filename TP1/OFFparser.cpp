#include<vector>
#include<string>
#include<triangle.h>
#include<vertex.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>

template<typename T, unsigned int DIMENTION, unsigned int EDGES>
void parseOFF(const std::string& filename,
              std::vector<Vertex<T, DIMENTION>>& vertices,
              std::vector<Triangle<T, DIMENTION, EDGES>>& triangles)
{
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("File won't open");

    std::string line;
    std::getline(file, line);
    if (line != "OFF") throw std::runtime_error("The file is not an OFF format");

    unsigned int numVertices = 0, numTriangles = 0;
    std::getline(file, line);
    std::istringstream header(line);
    header >> numVertices >> numTriangles;

    if (!header) throw std::runtime_error("Error when reading counts");

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

    triangles.resize(numTriangles);
    for (unsigned int i = 0; i < numTriangles; ++i) 
    {
        std::getline(file, line);
        std::istringstream ss(line);
        int edgesInFile;
        ss >> edgesInFile;

        if (edgesInFile != EDGES)
            throw std::runtime_error("Error when reading faces, number of edge is wrong");

        for(unsigned int e = 0; e < EDGES; ++e)
        {
            if (!ss >> triangles[i].triangleNeighboorIDAt(e))
                throw std::runtime_error("Error when reading faces coordinates");
        }
    }
}