#include <circulatorOnVertices.h>

CirculatorOnVertices::CirculatorOnVertices(Mesh& m, vertexIndex v, faceIndex startFace)
    : mesh(m), vIndex(v), currentFace(startFace), startFace(startFace), firstIteration(true)
{
    updateCurrentVertex();
}

Vertex<float,3>& CirculatorOnVertices::operator*() 
{
    return mesh.getVertices()[currentVertex];
}

CirculatorOnVertices& CirculatorOnVertices::operator++() {
    // Move to the next face around the vertex
    currentFace = nextIncidentFace();
    // Stop if we've returned to start or reached boundary
    if((currentFace == INVALID_FACE_INDEX) || ((currentFace == startFace) && !firstIteration))
    {
        currentVertex = INVALID_VERTEX_INDEX;
    } 
    else
    {
        updateCurrentVertex();
    }

    firstIteration = false;
    return *this;
}

bool CirculatorOnVertices::isValid() const 
{
    return currentVertex != INVALID_VERTEX_INDEX;
}

void CirculatorOnVertices::updateCurrentVertex() 
{
    auto& f = mesh.getFaces()[currentFace];

    // Find the index of vIndex in the current face
    for (unsigned int i = 0; i < 3; ++i) 
    {
        if (f.vertexAt(i) == vIndex)
        {
            currentVertex = f.vertexAt((i + 1) % 3);
            break;
        }
    }
}

faceIndex CirculatorOnVertices::nextIncidentFace() 
{
    auto& f = mesh.getFaces()[currentFace];

    // Find vIndex in the face and return the CCW neighbor
    for (unsigned int i = 0; i < 3; ++i) {
        if (f.vertexAt(i) == vIndex) {
            return f.facesNeighboorIDAt((i + 2) % 3);
        }
    }

    // No neighbor found
    return INVALID_FACE_INDEX;
}