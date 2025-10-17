#ifndef CIRCULATOR_ON_FACES_H
#define CIRCULATOR_ON_FACES_H

#include <mesh.h>

#define INVALID_FACE_INDEX static_cast<faceIndex>(-1)

class CirculatorOnFaces
{
    public:
    CirculatorOnFaces(Mesh& m, vertexIndex v, faceIndex startFace);

    Face<float, 3, 3>& operator*();

    CirculatorOnFaces& operator++();

    bool operator!=(const CirculatorOnFaces& other) const;

    bool isValid() const;

    inline faceIndex getCurrentFaceIndex() { return currentFace; }

private:
    Mesh& mesh;
    vertexIndex vIndex;
    faceIndex currentFace;
    faceIndex startFaceIndex;

    faceIndex nextIncidentFace();
};

#endif //CIRCULATOR_ON_FACES_H