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
#include <vector>
#include <ctime>

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

static inline double cotangent(const Point<double,3>& u, const Point<double,3>& v, const Point<double,3>& w)
{
    Point<double,3> wu = u - w;
    Point<double,3> wv = v - w;

    double dot = Point<double,3>::dotProduct(wu, wv);
    Point<double,3> cross = Point<double,3>::crossProduct(wu, wv);
    double crossNorm = cross.norm(); //orientation?

    return dot / crossNorm;
}

float Mesh::area(faceIndex f)
{
    auto & face = faces[f];

    Point<float,3> p0 = vertices[face.vertexAt(0)].getPosition();
    Point<float,3> p1 = vertices[face.vertexAt(1)].getPosition();
    Point<float,3> p2 = vertices[face.vertexAt(2)].getPosition();

    Point<float,3> u = p1 - p0;
    Point<float,3> v = p2 - p0;

    return 0.5f * Point<float,3>::crossProduct(u, v).norm();
}

float Mesh::computeCotangentLaplacian(vertexIndex v, const std::function<float(vertexIndex)>& projection)
{
    double acc = 0.0f;
    double area = 0.0f;

    for(auto cf = incidentFaces(v); cf.isValid(); ++cf)
    {
        area += this->area(cf.getCurrentFaceIndex()) / 3.0f; //cumulate a thrid of area

        //find local index of v in a
        vertexIndex lav = (*cf).vertexAt(0) == v ? 0 : ((*cf).vertexAt(1) == v ? 1 : 2);

        Point<double, 3> i(
            static_cast<double>(vertices[v].getPosition()[0]),
            static_cast<double>(vertices[v].getPosition()[1]),
            static_cast<double>(vertices[v].getPosition()[2])
        );
        Point<double, 3> j(
            static_cast<double>(vertices[(*cf).vertexAt((lav + 2) % 3)].getPosition()[0]),
            static_cast<double>(vertices[(*cf).vertexAt((lav + 2) % 3)].getPosition()[1]),
            static_cast<double>(vertices[(*cf).vertexAt((lav + 2) % 3)].getPosition()[2])
        );
        Point<double, 3> aij(
            static_cast<double>(vertices[(*cf).vertexAt((lav + 1) % 3)].getPosition()[0]),
            static_cast<double>(vertices[(*cf).vertexAt((lav + 1) % 3)].getPosition()[1]),
            static_cast<double>(vertices[(*cf).vertexAt((lav + 1) % 3)].getPosition()[2])
        );

        double alpha = cotangent(i, j, aij);

        auto & neighboorFace = faces[(*cf).facesNeighboorIDAt((lav +1) %3)];
        vertexIndex lbv = neighboorFace.vertexAt(0) == v ? 0 : (neighboorFace.vertexAt(1) == v ? 1 : 2);

        Point<double, 3> bij(
            static_cast<double>(vertices[neighboorFace.vertexAt((lbv + 2) % 3)].getPosition()[0]),
            static_cast<double>(vertices[neighboorFace.vertexAt((lbv + 2) % 3)].getPosition()[1]),
            static_cast<double>(vertices[neighboorFace.vertexAt((lbv + 2) % 3)].getPosition()[2])
        );

        double beta = cotangent(i, j, bij);

        //std::cout << alpha << " " << beta << std::endl;

        acc += (alpha + beta) * (projection(v) - projection((*cf).vertexAt((lav + 2) % 3)));
    }
    //std::cout << area <<  std::endl;
    return 1.0f / (2.0f * static_cast<float>(area)) * static_cast<float>(acc);
}

void Mesh::writeCotangentCurvatureOFF(const std::string& filename)
{
    
}

void Mesh::writeTemperatureOFF(const std::string& filename, unsigned int iteration)
{
    std::vector<float> temperature(vertices.size(), 0.0f);
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    size_t seedVertex = std::rand() % vertices.size();
    temperature[seedVertex] = 1.0f;
    std::vector<float> temp(vertices.size());
    for(unsigned int i = 0; i < iteration; ++i)
    {
        for(vertexIndex vertex = 0; vertex < vertices.size(); ++vertex)
        {
            float dt = 0.1f;
            temp[vertex] = temperature[vertex] + dt * computeCotangentLaplacian(vertex, [&temperature](vertexIndex vertex)
            { return temperature[vertex]; });
            // something degenerate the temperature...
            // assert(!std::isnan(temperature[vertex]));
        }
        //write COFF
        //normalisation
        float tmin = 0.0f;
        float tmax = 1.0f;
        float range = tmax - tmin;

        //COFF
        std::ostringstream fname;
        fname << filename << "_" << i << ".off";
        std::ofstream ofs(fname.str());
        if (!ofs) throw std::runtime_error("Cannot open file " + fname.str());

        ofs << "COFF\n";
        ofs << vertices.size() << " " << faces.size() << " 0\n";

        for (vertexIndex vertex = 0; vertex < vertices.size(); ++vertex)
        {
            const auto& p = vertices[vertex].getPosition();
            float t = temperature[vertex] / range;
            //std::cout << t << tmin << range << std::endl;

            int r = static_cast<int>(255 * t);
            int g = 0;
            int b = static_cast<int>(255 * (1.0f - t));

            ofs << p[0] << " " << p[1] << " " << p[2] << " "
                << r << " " << g << " " << b << " 255\n";
        }

        for (auto& f : faces)
        {
            ofs << 3 << " " << f.vertexAt(0) << " "
                       << f.vertexAt(1) << " "
                       << f.vertexAt(2) << "\n";
        }

        ofs.close();
        std::cout << "Iteration " << i << " -> " << fname.str() << " written." << std::endl;


        temperature = temp;
    }
}

static void HSVtoRGB(float h, float s, float v, int& r, int& g, int& b)
{
    float c = v * s;
    float hh = h / 60.0f;
    float x = c * (1 - fabs(fmod(hh, 2.0f) - 1));
    float m = v - c;

    float r_, g_, b_;
    if(hh >= 0 && hh < 1) { r_=c; g_=x; b_=0; }
    else if(hh < 2) { r_=x; g_=c; b_=0; }
    else if(hh < 3) { r_=0; g_=c; b_=x; }
    else if(hh < 4) { r_=0; g_=x; b_=c; }
    else if(hh < 5) { r_=x; g_=0; b_=c; }
    else { r_=c; g_=0; b_=x; }

    r = static_cast<int>((r_ + m) * 255);
    g = static_cast<int>((g_ + m) * 255);
    b = static_cast<int>((b_ + m) * 255);
}

void Mesh::writeCurvatureOFF(const std::string& filename)
{
    std::ofstream out(filename);
    if(!out) {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return;
    }
    out << "COFF\n";
    out << vertices.size() << " " << faces.size() << " 0\n";

    // Step 1: compute curvature for all vertices
    std::vector<float> curvature(vertices.size());
    float minH = std::numeric_limits<float>::max();
    float maxH = std::numeric_limits<float>::lowest();

    for(vertexIndex v = 0; v < vertices.size(); ++v)
    {
        float X = computeCotangentLaplacian(v, [&](vertexIndex vi)
        { return vertices[vi].getPosition()[0]; });
        float Y = computeCotangentLaplacian(v, [&](vertexIndex vi)
        { return vertices[vi].getPosition()[1]; });
        float Z = computeCotangentLaplacian(v, [&](vertexIndex vi)
        { return vertices[vi].getPosition()[2]; });

        Point<float, 3> Ds(X, Y, Z);
        float H = Ds.norm() / 2.0f;   // mean curvature
        curvature[v] = H;

        minH = std::min(minH, H);
        maxH = std::max(maxH, H);
    }

    for(vertexIndex v = 0; v < vertices.size(); ++v)
    {
        float H = curvature[v];

        float t = (H - minH) / (maxH - minH + 1e-8f);

        float hue = (1.0f - t) * 240.0f;
        float sat = 1.0f;
        float val = 1.0f;
        int r, g, b;
        HSVtoRGB(hue, sat, val, r, g, b);

        const auto& pos = vertices[v].getPosition();
        out << pos[0] << " " << pos[1] << " " << pos[2] << " "
            << r << " " << g << " " << b << " 255\n";
    }

    for (auto& f : faces)
    {
        out << 3 << " " << f.vertexAt(0) << " "
                   << f.vertexAt(1) << " "
                   << f.vertexAt(2) << "\n";
    }

    out.close();
    std::cout << "Curvature " << filename << " written." << std::endl;
}


void Mesh::faceSplit(faceIndex f, vertexIndex v)
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

void Mesh::edgeSplit(faceIndex f, vertexIndex v)
{

}

void Mesh::flip(faceIndex fl, faceIndex fr)
{
/*
    A             A 
   / \           /|\
  /   \         / | \
 B-----C   ->  B  |  C
  \   /         \ | /
   \ /           \|/
    D             D
*/
}