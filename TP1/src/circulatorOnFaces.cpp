#include <circulatorOnFaces.h>

CirculatorOnFaces::CirculatorOnFaces(Mesh& m, vertexIndex v, faceIndex startFace)
        : mesh(m), vIndex(v), currentFace(startFace), startFaceIndex(startFace), firstIteration(true) {}

    // Accès à la face courante
    Face<float, 3, 3>& CirculatorOnFaces::operator*() 
    {
        return mesh.getFaces()[currentFace];
    }

    // Passer à la face suivante autour du sommet
    CirculatorOnFaces& CirculatorOnFaces::operator++() 
    {
        // Trouver le prochain voisin
        currentFace = nextIncidentFace();
        // Détection de boucle (retour au départ)
        if (currentFace == startFaceIndex) {
            if (!firstIteration) {
                currentFace = INVALID_FACE_INDEX; // fin de circulation
            }
        }
        firstIteration = false;
        return *this;
    }

    // Vérifier validité
    bool CirculatorOnFaces::operator!=(const CirculatorOnFaces& other) const 
    {
        return currentFace != other.currentFace;
    }

    bool CirculatorOnFaces::isValid() const { return currentFace != INVALID_FACE_INDEX; }



    faceIndex CirculatorOnFaces::nextIncidentFace() 
    {
        auto& f = mesh.getFaces()[currentFace];

        // On cherche le sommet vIndex dans la face courante
        for (unsigned int i = 0; i < 3; i++) {
            if (f.vertexAt(i) == vIndex) {
                // Le voisin qui partage l’arête opposée
                return f.facesNeighboorIDAt((i + 2) % 3);
            }
        }
        return INVALID_FACE_INDEX;
    }