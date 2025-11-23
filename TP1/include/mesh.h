#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <functional>

#include <vertex.h>
#include <face.h>

class CirculatorOnFaces;
class CirculatorOnVertices;

class Mesh
{
public:
    void loadOffFile(const std::string &file);
    void writeOffFile(const std::string &file);
    bool check();

    inline std::vector<Vertex<float, 3>> &getVertices() { return vertices; }
    inline std::vector<Face<float, 3, 3>> &getFaces() { return faces; }

    std::vector<Face<float, 3, 3>> getFaceNeighboors(vertexIndex index);
    std::vector<Vertex<float, 3>> getVertexNeighboors(vertexIndex index);

    CirculatorOnFaces incidentFaces(vertexIndex vi);
    CirculatorOnVertices adjacentVertices(vertexIndex vi);

    float area(faceIndex f);

    float computeCotangentLaplacian(vertexIndex v, const std::function<float(vertexIndex)> &projection);

    void writeCotangentCurvatureOFF(const std::string &filename);
    void writeTemperatureOFF(const std::string &filename, unsigned int iteration);
    void writeCurvatureOFF(const std::string &filename);

    void faceSplit(faceIndex f, vertexIndex d);
    void edgeSplit(faceIndex f, vertexIndex d, vertexIndex a, vertexIndex b);
    void flip(faceIndex fl, faceIndex fr);

    float orientationTest(Point<float, 3> &a, Point<float, 3> &b, Point<float, 3> &c);
    float inTriangle(const Point<float, 3> &p,
                     const Point<float, 3> &a,
                     const Point<float, 3> &b, 
                     const Point<float, 3> &c);
    float inTriangle(Point<float, 3> &p, Face<float, 3, 3> &f);

    float inCircle(Point<float, 3> &p,
                   Point<float, 3> &a,
                   Point<float, 3> &b,
                   Point<float, 3> &c);
    float inSphere(Point<float, 3> &p,
                   Point<float, 3> &a,
                   Point<float, 3> &b,
                   Point<float, 3> &c,
                   Point<float, 3> &d);
    void makeDelaunay2D();
    void insertPoint(Point<float,3> p);

    bool isFictitious(Face<float, 3, 3> &f);
    static bool isFictitious(Point<float, 3> &p);

    float localSurfaceAround(vertexIndex v);
    void mergeCloseVertices(float threshold, float surfaceThreshold = 0.05f);

private:
    std::vector<Vertex<float, 3>> vertices;
    std::vector<Face<float, 3, 3>> faces;
};

#endif // MESH_H