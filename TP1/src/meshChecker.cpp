#include <face.h>
#include <vertex.h>
#include <vector>
#include <iostream>

template <typename T, unsigned int DIMENTION, unsigned int VERTICE_NUMBER>
bool checkMesh(std::vector<Face<T, DIMENTION, VERTICE_NUMBER>> &faces)
{
    for(unsigned int i = 0; i < faces.size(); ++i)
    {
        if(!checkNeighboor(faces, (faceIndex)i))
        {
            return false;
        }
    }
    return true;
}

template <typename T, unsigned int DIMENTION, unsigned int VERTICE_NUMBER>
bool checkNeighboor(std::vector<Face<T, DIMENTION, VERTICE_NUMBER>> &faces,
                    const faceIndex index)
{
    Face<T, DIMENTION, VERTICE_NUMBER> & me = faces[(int)index];
    for(unsigned int i = 0; i < 3; ++i)
    {
        if(faces[(int)me.facesNeighboorIDAt(i)].facesNeighboorIDAt(0) != index
            && faces[(int)me.facesNeighboorIDAt(i)].facesNeighboorIDAt(1) != index
            && faces[(int)me.facesNeighboorIDAt(i)].facesNeighboorIDAt(2) != index)
        {
            std::cout << "Mesh error at : " << i << std::endl;
            return false;
        }
    }
    return true;
}