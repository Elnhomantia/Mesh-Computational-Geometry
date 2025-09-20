#include <mesh.h>
#include "OFFparser.cpp"
#include "meshChecker.cpp"

#include <set>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iostream>

#include <circulatorOnFaces.h>
#include <circulatorOnVertices.h>

void Mesh::loadOffFile(const std::string & file)
{
    parseOFF<float, 3, 3>(file, vertices, faces);
}

bool Mesh::check()
{
    return checkMesh(faces);
}

CirculatorOnFaces Mesh::incidentFaces(vertexIndex vi) 
{
    faceIndex f = vertices[vi].faceNeighboor();
    return CirculatorOnFaces(*this, vi, f);
}
CirculatorOnVertices Mesh::adjacent_vertices(vertexIndex vi) 
{
    faceIndex f = vertices[vi].faceNeighboor();
    return CirculatorOnVertices(*this, vi, f);
}

std::vector<Point<float,3>> Mesh::computeUniformLaplacian()
{
    std::vector<Point<float,3>> laplacian(vertices.size());

    for (vertexIndex vi = 0; vi < vertices.size(); ++vi)
    {
        const Vertex<float,3>& v = vertices[vi];
        Point<float,3> sum; // defaults to (0,0,0)
        int degree = 0;

        for (auto circ = adjacent_vertices(vi); circ.isValid(); ++circ) {
            const Vertex<float,3>& neighbor = *circ;
            sum += neighbor.getPosition();
            ++degree;
        }

        if (degree > 0)
            laplacian[vi] = (sum / float(degree)) - v.getPosition();
        else
            laplacian[vi] = Point<float,3>(0.0f,0.0f,0.0f); // isolated vertex
    }

    return laplacian;
}

void Mesh::writeCurvatureOFF(const std::string& filename)
{
    // Step 1: compute Laplacian vectors
    std::vector<Point<float,3>> lap = computeUniformLaplacian();

    // Step 2: compute curvature magnitude
    std::vector<float> curvature(vertices.size());
    float minC = 1e30f, maxC = -1e30f;
    for (size_t i = 0; i < vertices.size(); ++i) {
        curvature[i] = lap[i].norm();
        minC = std::min(minC, curvature[i]);
        maxC = std::max(maxC, curvature[i]);
    }

    // Step 3: open OFF file with color
    std::ofstream ofs(filename);
    ofs << "COFF\n";
    ofs << vertices.size() << " " << faces.size() << " 0\n";

    // Step 4: write vertices with RGB based on curvature
    for (size_t i = 0; i < vertices.size(); ++i) {
        const auto& p = vertices[i].getPosition();

        // normalize curvature to [0,1]
        float c = (curvature[i] - minC) / (maxC - minC);

        // map to HSV color (blue low curvature, red high)
        float H = (1.0f - c) * 0.66f; // hue [0.66=blue -> 0=red]
        float S = 1.0f;
        float V = 1.0f;

        // HSV -> RGB conversion
        float r,g,b;
        int iH = int(H*6);
        float f = H*6 - iH;
        float pV = V*(1-S);
        float q = V*(1-f*S);
        float t = V*(1-(1-f)*S);

        switch(iH % 6){
            case 0: r=V; g=t; b=pV; break;
            case 1: r=q; g=V; b=pV; break;
            case 2: r=pV; g=V; b=t; break;
            case 3: r=pV; g=q; b=V; break;
            case 4: r=t; g=pV; b=V; break;
            case 5: r=V; g=pV; b=q; break;
        }

        ofs << p[0] << " " << p[1] << " " << p[2] << " "
            << int(r*255) << " " << int(g*255) << " " << int(b*255) << " 255\n";
    }

    // Step 5: write faces
    for (size_t i = 0; i < faces.size(); ++i) {
        auto& f = faces[i];
        ofs << 3 << " " << f.vertexAt(0) << " " << f.vertexAt(1) << " " << f.vertexAt(2) << "\n";
    }

    ofs.close();
    std::cout << "OFF with curvature coloring written to " << filename << std::endl;
}

void Mesh::split(faceIndex f, vertexIndex v)
{
    //v must be positionned in f
    //ABC -> ABD + ADC + BCD
    /*
           C
          / \
         / | \
        /  D  \
       / /   \ \
      A ------- B
    */
    Face<float, 3, 3> tr1;
    Face<float, 3, 3> tr2;
    Face<float, 3, 3> tr3;

    tr1.vertexAt(0) = faces[f].vertexAt(0); //A
    tr1.vertexAt(1) = faces[f].vertexAt(1); //B
    tr1.vertexAt(2) = v; //D

    tr2.vertexAt(0) = faces[f].vertexAt(0); //A
    tr2.vertexAt(1) = v; //D
    tr2.vertexAt(2) = faces[f].vertexAt(2); //C

    tr3.vertexAt(0) = faces[f].vertexAt(1); //B
    tr3.vertexAt(1) = faces[f].vertexAt(2); //C
    tr3.vertexAt(2) = v; //D

    //Add new triangles to the pool.
    faces.push_back(tr2);
    faceIndex tr2Index = faces.size() -1;
    faces.push_back(tr3);
    faceIndex tr3Index = faces.size() -1;

    tr1.facesNeighboorIDAt(0) = tr3Index;//BCD
    tr1.facesNeighboorIDAt(1) = tr2Index;//ADC
    tr1.facesNeighboorIDAt(2) = faces[f].facesNeighboorIDAt(2);//ABC neighbbor opposed to C

    tr2.facesNeighboorIDAt(0) = tr3Index;//BCD
    tr2.facesNeighboorIDAt(1) = faces[f].facesNeighboorIDAt(1);//ABC neighbbor opposed to B
    tr2.facesNeighboorIDAt(2) = f;//ABD

    tr3.facesNeighboorIDAt(0) = tr2Index;//ADC
    tr3.facesNeighboorIDAt(1) = f;//ABD
    tr3.facesNeighboorIDAt(2) = faces[f].facesNeighboorIDAt(0);//ABC neighbbor opposed to 1

    //replace source triangle (Must be at the end, we still use neighboor from it)
    //can't pop an element from vector, it will change indices making the whole mesh wrong.
    faces[f] = tr1;
}