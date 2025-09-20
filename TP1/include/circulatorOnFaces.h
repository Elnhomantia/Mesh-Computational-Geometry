#ifndef CIRCULATOR_ON_FACES_H
#define CIRCULATOR_ON_FACES_H

#include <mesh.h>

#define INVALID_FACE_INDEX static_cast<faceIndex>(-1)

class CirculatorOnFaces
{
    public:
    CirculatorOnFaces(Mesh& m, vertexIndex v, faceIndex startFace);

    // Accès à la face courante
    Face<float, 3, 3>& operator*();

    // Passer à la face suivante autour du sommet
    CirculatorOnFaces& operator++();

    // Vérifier validité
    bool operator!=(const CirculatorOnFaces& other) const;

    bool isValid() const;

private:
    Mesh& mesh;
    vertexIndex vIndex;
    faceIndex currentFace;
    faceIndex startFaceIndex;
    bool firstIteration;

    faceIndex nextIncidentFace();
};

#endif //CIRCULATOR_ON_FACES_H