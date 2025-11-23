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
#include <queue>

void Mesh::loadOffFile(const std::string &file)
{
    parseOFF<float, 3, 3>(file, vertices, faces);
}
void Mesh::writeOffFile(const std::string &file)
{
    writeOFF(file, vertices, faces); 
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
CirculatorOnVertices Mesh::adjacentVertices(vertexIndex vi)
{
    faceIndex f = vertices[vi].faceNeighboor();
    return CirculatorOnVertices(*this, vi, f);
}

static inline double cotangent(const Point<double, 3> &u, const Point<double, 3> &v, const Point<double, 3> &w)
{
    Point<double, 3> wu = u - w;
    Point<double, 3> wv = v - w;

    double dot = Point<double, 3>::dotProduct(wu, wv);
    Point<double, 3> cross = Point<double, 3>::crossProduct(wu, wv);
    double crossNorm = cross.norm(); // orientation?

    return dot / crossNorm;
}

float Mesh::area(faceIndex f)
{
    auto &face = faces[f];

    Point<float, 3> p0 = vertices[face.vertexAt(0)].getPosition();
    Point<float, 3> p1 = vertices[face.vertexAt(1)].getPosition();
    Point<float, 3> p2 = vertices[face.vertexAt(2)].getPosition();

    Point<float, 3> u = p1 - p0;
    Point<float, 3> v = p2 - p0;

    return 0.5f * Point<float, 3>::crossProduct(u, v).norm();
}

float Mesh::computeCotangentLaplacian(vertexIndex v, const std::function<float(vertexIndex)> &projection)
{
    double acc = 0.0f;
    double area = 0.0f;

    for (auto cf = incidentFaces(v); cf.isValid(); ++cf)
    {
        area += this->area(cf.getCurrentFaceIndex()) / 3.0f; // cumulate a thrid of area

        // find local index of v in a
        vertexIndex lav = (*cf).vertexAt(0) == v ? 0 : ((*cf).vertexAt(1) == v ? 1 : 2);

        Point<double, 3> i(
            static_cast<double>(vertices[v].getPosition()[0]),
            static_cast<double>(vertices[v].getPosition()[1]),
            static_cast<double>(vertices[v].getPosition()[2]));
        Point<double, 3> j(
            static_cast<double>(vertices[(*cf).vertexAt((lav + 2) % 3)].getPosition()[0]),
            static_cast<double>(vertices[(*cf).vertexAt((lav + 2) % 3)].getPosition()[1]),
            static_cast<double>(vertices[(*cf).vertexAt((lav + 2) % 3)].getPosition()[2]));
        Point<double, 3> aij(
            static_cast<double>(vertices[(*cf).vertexAt((lav + 1) % 3)].getPosition()[0]),
            static_cast<double>(vertices[(*cf).vertexAt((lav + 1) % 3)].getPosition()[1]),
            static_cast<double>(vertices[(*cf).vertexAt((lav + 1) % 3)].getPosition()[2]));

        double alpha = cotangent(i, j, aij);

        auto &neighboorFace = faces[(*cf).facesNeighboorIDAt((lav + 1) % 3)];
        vertexIndex lbv = neighboorFace.vertexAt(0) == v ? 0 : (neighboorFace.vertexAt(1) == v ? 1 : 2);

        Point<double, 3> bij(
            static_cast<double>(vertices[neighboorFace.vertexAt((lbv + 2) % 3)].getPosition()[0]),
            static_cast<double>(vertices[neighboorFace.vertexAt((lbv + 2) % 3)].getPosition()[1]),
            static_cast<double>(vertices[neighboorFace.vertexAt((lbv + 2) % 3)].getPosition()[2]));

        double beta = cotangent(i, j, bij);

        // std::cout << alpha << " " << beta << std::endl;

        acc += (alpha + beta) * (projection(v) - projection((*cf).vertexAt((lav + 2) % 3)));
    }
    // std::cout << area <<  std::endl;
    return 1.0f / (2.0f * static_cast<float>(area)) * static_cast<float>(acc);
}

void Mesh::writeCotangentCurvatureOFF(const std::string &filename)
{
}

void Mesh::writeTemperatureOFF(const std::string &filename, unsigned int iteration)
{
    std::vector<float> temperature(vertices.size(), 0.0f);
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    size_t seedVertex = std::rand() % vertices.size();
    temperature[seedVertex] = 1.0f;
    std::vector<float> temp(vertices.size());
    for (unsigned int i = 0; i < iteration; ++i)
    {
        for (vertexIndex vertex = 0; vertex < vertices.size(); ++vertex)
        {
            float dt = 0.1f;
            temp[vertex] = temperature[vertex] + dt * computeCotangentLaplacian(vertex, [&temperature](vertexIndex vertex)
                                                                                { return temperature[vertex]; });
            // something degenerate the temperature...
            // assert(!std::isnan(temperature[vertex]));
        }
        // write COFF
        // normalisation
        float tmin = 0.0f;
        float tmax = 1.0f;
        float range = tmax - tmin;

        // COFF
        std::ostringstream fname;
        fname << filename << "_" << i << ".off";
        std::ofstream ofs(fname.str());
        if (!ofs)
            throw std::runtime_error("Cannot open file " + fname.str());

        ofs << "COFF\n";
        ofs << vertices.size() << " " << faces.size() << " 0\n";

        for (vertexIndex vertex = 0; vertex < vertices.size(); ++vertex)
        {
            const auto &p = vertices[vertex].getPosition();
            float t = temperature[vertex] / range;
            // std::cout << t << tmin << range << std::endl;

            int r = static_cast<int>(255 * t);
            int g = 0;
            int b = static_cast<int>(255 * (1.0f - t));

            ofs << p[0] << " " << p[1] << " " << p[2] << " "
                << r << " " << g << " " << b << " 255\n";
        }

        for (auto &f : faces)
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

static void HSVtoRGB(float h, float s, float v, int &r, int &g, int &b)
{
    float c = v * s;
    float hh = h / 60.0f;
    float x = c * (1 - fabs(fmod(hh, 2.0f) - 1));
    float m = v - c;

    float r_, g_, b_;
    if (hh >= 0 && hh < 1)
    {
        r_ = c;
        g_ = x;
        b_ = 0;
    }
    else if (hh < 2)
    {
        r_ = x;
        g_ = c;
        b_ = 0;
    }
    else if (hh < 3)
    {
        r_ = 0;
        g_ = c;
        b_ = x;
    }
    else if (hh < 4)
    {
        r_ = 0;
        g_ = x;
        b_ = c;
    }
    else if (hh < 5)
    {
        r_ = x;
        g_ = 0;
        b_ = c;
    }
    else
    {
        r_ = c;
        g_ = 0;
        b_ = x;
    }

    r = static_cast<int>((r_ + m) * 255);
    g = static_cast<int>((g_ + m) * 255);
    b = static_cast<int>((b_ + m) * 255);
}

void Mesh::writeCurvatureOFF(const std::string &filename)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return;
    }
    out << "COFF\n";
    out << vertices.size() << " " << faces.size() << " 0\n";

    // Step 1: compute curvature for all vertices
    std::vector<float> curvature(vertices.size());
    float minH = std::numeric_limits<float>::max();
    float maxH = std::numeric_limits<float>::lowest();

    for (vertexIndex v = 0; v < vertices.size(); ++v)
    {
        float X = computeCotangentLaplacian(v, [&](vertexIndex vi)
                                            { return vertices[vi].getPosition()[0]; });
        float Y = computeCotangentLaplacian(v, [&](vertexIndex vi)
                                            { return vertices[vi].getPosition()[1]; });
        float Z = computeCotangentLaplacian(v, [&](vertexIndex vi)
                                            { return vertices[vi].getPosition()[2]; });

        Point<float, 3> Ds(X, Y, Z);
        float H = Ds.norm() / 2.0f; // mean curvature
        curvature[v] = H;

        minH = std::min(minH, H);
        maxH = std::max(maxH, H);
    }

    for (vertexIndex v = 0; v < vertices.size(); ++v)
    {
        float H = curvature[v];

        float t = (H - minH) / (maxH - minH + 1e-8f);

        float hue = (1.0f - t) * 240.0f;
        float sat = 1.0f;
        float val = 1.0f;
        int r, g, b;
        HSVtoRGB(hue, sat, val, r, g, b);

        const auto &pos = vertices[v].getPosition();
        out << pos[0] << " " << pos[1] << " " << pos[2] << " "
            << r << " " << g << " " << b << " 255\n";
    }

    for (auto &f : faces)
    {
        out << 3 << " " << f.vertexAt(0) << " "
            << f.vertexAt(1) << " "
            << f.vertexAt(2) << "\n";
    }

    out.close();
    std::cout << "Curvature " << filename << " written." << std::endl;
}

void Mesh::faceSplit(faceIndex f, vertexIndex d)
{
    // vertex d must be positionned in face f
    // ABC -> ABD + ADC + BCD
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

    tr1.vertexAt(0) = faces[f].vertexAt(0); // A
    tr1.vertexAt(1) = faces[f].vertexAt(1); // B
    tr1.vertexAt(2) = d;                    // D

    tr2.vertexAt(0) = faces[f].vertexAt(0); // A
    tr2.vertexAt(1) = d;                    // D
    tr2.vertexAt(2) = faces[f].vertexAt(2); // C

    tr3.vertexAt(0) = faces[f].vertexAt(1); // B
    tr3.vertexAt(1) = faces[f].vertexAt(2); // C
    tr3.vertexAt(2) = d;                    // D

    // Add new triangles to the pool.
    faces.push_back(tr2);
    faceIndex tr2Index = faces.size() - 1;
    faces.push_back(tr3);
    faceIndex tr3Index = faces.size() - 1;

    tr1.facesNeighboorIDAt(0) = tr3Index;                       // BCD
    tr1.facesNeighboorIDAt(1) = tr2Index;                       // ADC
    tr1.facesNeighboorIDAt(2) = faces[f].facesNeighboorIDAt(2); // ABC neighbbor opposed to C

    tr2.facesNeighboorIDAt(0) = tr3Index;                       // BCD
    tr2.facesNeighboorIDAt(1) = faces[f].facesNeighboorIDAt(1); // ABC neighbbor opposed to B
    tr2.facesNeighboorIDAt(2) = f;                              // ABD

    tr3.facesNeighboorIDAt(0) = tr2Index;                       // ADC
    tr3.facesNeighboorIDAt(1) = f;                              // ABD
    tr3.facesNeighboorIDAt(2) = faces[f].facesNeighboorIDAt(0); // ABC neighbbor opposed to 1

    // replace source triangle (Must be at the end, we still use neighboor from it)
    // can't pop an element from vector, it will change indices making the whole mesh wrong.
    faces[f] = tr1;
}

void Mesh::edgeSplit(faceIndex f, vertexIndex d, vertexIndex a, vertexIndex b)
{
    // vertex d must be positioned on edge (a,b)
    // ABC -> ADC + BDC
    // If AB is shared by another face (ABX), also split it into ADX + BDX
    /*
           C
          /|\
         / | \
        /  |  \
       A---D---B
        \  |  /
         \ | /
          \|/
           X
    */

    // find C
    vertexIndex c;
    int idxA = -1, idxB = -1, idxC = -1;
    for (int i = 0; i < 3; i++)
    {
        vertexIndex v = faces[f].vertexAt(i);
        if (v == a)
            idxA = i;
        else if (v == b)
            idxB = i;
        else
        {
            c = v;
            idxC = i;
        }
    }

    // neighbor across AB (opposite to C)
    faceIndex f2 = faces[f].facesNeighboorIDAt(idxC);
    bool boundary = (f2 == -1);

    // create triangles inside f
    Face<float, 3, 3> tr1; // ADC
    Face<float, 3, 3> tr2; // BDC

    tr1.vertexAt(0) = a;
    tr1.vertexAt(1) = d;
    tr1.vertexAt(2) = c;

    tr2.vertexAt(0) = b;
    tr2.vertexAt(1) = c;
    tr2.vertexAt(2) = d;

    // keep old neighbors
    faceIndex neighOppA = faces[f].facesNeighboorIDAt(idxA);
    faceIndex neighOppB = faces[f].facesNeighboorIDAt(idxB);
    faceIndex neighOppC = faces[f].facesNeighboorIDAt(idxC);

    // push tr2 (BDC)
    faces.push_back(tr2);
    faceIndex tr2Index = faces.size() - 1;

    // replace f with tr1 (ADC)
    faces[f] = tr1;

    // if AB has a neighbor, split it too
    faceIndex tr3Index = -1;
    faceIndex tr4Index = -1;
    if (!boundary)
    {
        // find X in f2 (third vertex, not a or b)
        vertexIndex x;
        int idxA2 = -1, idxB2 = -1, idxX = -1;
        for (int i = 0; i < 3; i++)
        {
            vertexIndex v = faces[f2].vertexAt(i);
            if (v == a)
                idxA2 = i;
            else if (v == b)
                idxB2 = i;
            else
            {
                x = v;
                idxX = i;
            }
        }

        faceIndex neighOppA2 = faces[f2].facesNeighboorIDAt(idxA2);
        faceIndex neighOppB2 = faces[f2].facesNeighboorIDAt(idxB2);
        faceIndex neighOppX = faces[f2].facesNeighboorIDAt(idxX);

        // ADX
        Face<float, 3, 3> tr3;
        tr3.vertexAt(0) = a;
        tr3.vertexAt(1) = d;
        tr3.vertexAt(2) = x;

        // BDX
        Face<float, 3, 3> tr4;
        tr4.vertexAt(0) = b;
        tr4.vertexAt(1) = x;
        tr4.vertexAt(2) = d;

        // push tr4
        faces.push_back(tr4);
        tr4Index = faces.size() - 1;

        // replace f2 with tr3
        faces[f2] = tr3;
        tr3Index = f2;

        //--- Step 5b: set neighbors for tr3 (ADX)
        faces[tr3Index].facesNeighboorIDAt(0) = tr4Index;  // opposite A (edge DX) → BDX
        faces[tr3Index].facesNeighboorIDAt(1) = neighOppX; // opposite D (edge AX) → old neighbor
        faces[tr3Index].facesNeighboorIDAt(2) = f;         // opposite X (edge AD) → ADC

        // For tr4 (BDX)
        faces[tr4Index].facesNeighboorIDAt(0) = tr3Index;  // opposite B (edge DX) → ADX
        faces[tr4Index].facesNeighboorIDAt(1) = tr2Index;  // opposite X (edge BD) → BDC
        faces[tr4Index].facesNeighboorIDAt(2) = neighOppX; // opposite D (edge BX) → old neighbor
    }

    // set neighbors for tr1 (ADC) and tr2 (BDC)
    faces[f].facesNeighboorIDAt(0) = tr2Index;                        // opposite A (edge DC) → BDC
    faces[f].facesNeighboorIDAt(1) = neighOppC;                       // opposite D (edge AC) → old neighbor
    faces[f].facesNeighboorIDAt(2) = boundary ? neighOppA : tr3Index; // opposite C (edge AD) → ADX if exists

    faces[tr2Index].facesNeighboorIDAt(0) = f;                               // opposite B (edge DC) → ADC
    faces[tr2Index].facesNeighboorIDAt(1) = boundary ? neighOppB : tr4Index; // opposite C (edge BD) → BDX if exists
    faces[tr2Index].facesNeighboorIDAt(2) = neighOppC;                       // opposite D (edge BC) → old neighbor
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
    //--- find common vertices
    std::vector<vertexIndex> common;
    vertexIndex onlyFl = -1, onlyFr = -1;
    for (int i = 0; i < 3; ++i)
    {
        vertexIndex v = faces[fl].vertexAt(i);
        bool inFr = false;
        for (int j = 0; j < 3; ++j)
            if (v == faces[fr].vertexAt(j))
                inFr = true;
        if (inFr)
            common.push_back(v);
        else
            onlyFl = v;
    }
    for (int i = 0; i < 3; ++i)
    {
        vertexIndex v = faces[fr].vertexAt(i);
        bool inFl = false;
        for (int j = 0; j < 3; ++j)
            if (v == faces[fl].vertexAt(j))
                inFl = true;
        if (!inFl)
            onlyFr = v;
    }

    // not adjacent or degenerate
    if (common.size() != 2 || onlyFl == -1 || onlyFr == -1)
    {
        // faces are not sharing a single edge -> nothing to do
        return;
    }

    // name vertices: A = onlyFl, D = onlyFr
    vertexIndex A = onlyFl;
    vertexIndex D = onlyFr;

    //--- find indices inside fl so we can determine u and v in CCW order:
    int iA = -1;
    for (int i = 0; i < 3; ++i)
        if (faces[fl].vertexAt(i) == A)
        {
            iA = i;
            break;
        }
    int iU = (iA + 1) % 3;
    int iV = (iA + 2) % 3;
    vertexIndex u = faces[fl].vertexAt(iU);
    vertexIndex v = faces[fl].vertexAt(iV);
    // now fl is (A, u, v) in CCW order

    // find indices of D, u, v inside fr
    int jD = -1, jU = -1, jV = -1;
    for (int i = 0; i < 3; ++i)
    {
        vertexIndex vv = faces[fr].vertexAt(i);
        if (vv == D)
            jD = i;
        else if (vv == u)
            jU = i;
        else if (vv == v)
            jV = i;
    }
    if (jD == -1 || jU == -1 || jV == -1)
    {
        // inconsistent topology
        return;
    }

    //--- save old external neighbors (may be -1 for boundary)
    // fl neighbors:
    faceIndex neigh_fl_iU = faces[fl].facesNeighboorIDAt(iU); // across (v,A)
    faceIndex neigh_fl_iV = faces[fl].facesNeighboorIDAt(iV); // across (A,u)
    // fr neighbors:
    faceIndex neigh_fr_jU = faces[fr].facesNeighboorIDAt(jU); // across (D,v)
    faceIndex neigh_fr_jV = faces[fr].facesNeighboorIDAt(jV); // across (u,D)

    //--- build new triangles (keep CCW)
    Face<float, 3, 3> fNew; // will replace faces[fl] : (A, u, D)
    fNew.vertexAt(0) = A;
    fNew.vertexAt(1) = u;
    fNew.vertexAt(2) = D;

    Face<float, 3, 3> gNew; // will replace faces[fr] : (D, v, A)
    gNew.vertexAt(0) = D;
    gNew.vertexAt(1) = v;
    gNew.vertexAt(2) = A;

    //--- set neighbors for new faces using saved old neighbors
    // fNew (A, u, D)
    //  opposite A -> face across (u,D)  == old neigh_fr_jV
    //  opposite u -> the other new face (fr)
    //  opposite D -> face across (A,u)  == old neigh_fl_iV
    fNew.facesNeighboorIDAt(0) = neigh_fr_jV;
    fNew.facesNeighboorIDAt(1) = fr;
    fNew.facesNeighboorIDAt(2) = neigh_fl_iV;

    // gNew (D, v, A)
    //  opposite D -> face across (v,A)  == old neigh_fl_iU
    //  opposite v -> the other new face (fl)
    //  opposite A -> face across (D,v)  == old neigh_fr_jU
    gNew.facesNeighboorIDAt(0) = neigh_fl_iU;
    gNew.facesNeighboorIDAt(1) = fl;
    gNew.facesNeighboorIDAt(2) = neigh_fr_jU;

    // replace in-place to keep indices
    faces[fl] = fNew;
    faces[fr] = gNew;

    // neigh across (v,A) used to point to fl -> must now point to fr
    if (neigh_fl_iU != -1)
    {
        for (int k = 0; k < 3; ++k)
        {
            if (faces[neigh_fl_iU].facesNeighboorIDAt(k) == fl)
            {
                faces[neigh_fl_iU].facesNeighboorIDAt(k) = fr;
                // do not break: handle possible degenerate repeats
            }
        }
    }

    // neigh across (u,D) used to point to fr -> must now point to fl
    if (neigh_fr_jV != -1)
    {
        for (int k = 0; k < 3; ++k)
        {
            if (faces[neigh_fr_jV].facesNeighboorIDAt(k) == fr)
            {
                faces[neigh_fr_jV].facesNeighboorIDAt(k) = fl;
            }
        }
    }

    // Note: neigh_fl_iV (across A,u) keeps pointing to fl (no change needed),
    //       neigh_fr_jU (across D,v) keeps pointing to fr (no change needed).
}

float Mesh::orientationTest(Point<float, 3> &a, Point<float, 3> &b, Point<float, 3> &c)
{
    // Only 2D, consider z = 0;
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

float Mesh::inTriangle(const Point<float, 3> &p,
                       const Point<float, 3> &a,
                       const Point<float, 3> &b,
                       const Point<float, 3> &c)
{
    float denom = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1]);

    // Compute barycentric coordinates
    float alpha = ((b[1] - c[1]) * (p[0] - c[0]) + (c[0] - b[0]) * (p[1] - c[1])) / denom;
    float beta = ((c[1] - a[1]) * (p[0] - c[0]) + (a[0] - c[0]) * (p[1] - c[1])) / denom;
    float gamma = 1.0f - alpha - beta;

    // On boundary?
    if (alpha == 0.0f || beta == 0.0f || gamma == 0.0f)
        return 0.0f;

    // Inside if all positive, outside otherwise
    if (alpha > 0.0f && beta > 0.0f && gamma > 0.0f)
        return 1.0f;
    else
        return -1.0f;
}
float Mesh::inTriangle(Point<float, 3> &p, Face<float, 3, 3> &f)
{
    return inTriangle(p,
                      vertices[f.vertexAt(0)].getPosition(),
                      vertices[f.vertexAt(1)].getPosition(),
                      vertices[f.vertexAt(2)].getPosition());
}

float Mesh::inCircle(Point<float, 3> &p,
                     Point<float, 3> &a,
                     Point<float, 3> &b,
                     Point<float, 3> &c)
{
    // Geometric predicate : test if point p lies inside the circumcircle of triangle (a,b,c)
    // > 0 : inside
    // < 0 : outside
    // = 0 : on the circle
    // Works in 2D (z ignored)

    float ax = a[0] - p[0];
    float ay = a[1] - p[1];
    float bx = b[0] - p[0];
    float by = b[1] - p[1];
    float cx = c[0] - p[0];
    float cy = c[1] - p[1];

    float det =
        (ax * ax + ay * ay) * (bx * cy - cx * by) - (bx * bx + by * by) * (ax * cy - cx * ay) + (cx * cx + cy * cy) * (ax * by - bx * ay);

    return det;
}

float Mesh::inSphere(Point<float, 3> &p,
                     Point<float, 3> &a,
                     Point<float, 3> &b,
                     Point<float, 3> &c,
                     Point<float, 3> &d)
{
    // Geometric predicate : test if point p lies inside the circumsphere of tetrahedron (a,b,c,d)
    // > 0 : inside
    // < 0 : outside
    // = 0 : on the sphere

    float ax = a[0] - p[0];
    float ay = a[1] - p[1];
    float az = a[2] - p[2];
    float bx = b[0] - p[0];
    float by = b[1] - p[1];
    float bz = b[2] - p[2];
    float cx = c[0] - p[0];
    float cy = c[1] - p[1];
    float cz = c[2] - p[2];
    float dx = d[0] - p[0];
    float dy = d[1] - p[1];
    float dz = d[2] - p[2];

    float det =
        (ax * ax + ay * ay + az * az) * (bx * (cy * dz - cz * dy) - by * (cx * dz - cz * dx) + bz * (cx * dy - cy * dx)) - (bx * bx + by * by + bz * bz) * (ax * (cy * dz - cz * dy) - ay * (cx * dz - cz * dx) + az * (cx * dy - cy * dx)) + (cx * cx + cy * cy + cz * cz) * (ax * (by * dz - bz * dy) - ay * (bx * dz - bz * dx) + az * (bx * dy - by * dx)) - (dx * dx + dy * dy + dz * dz) * (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx));

    return det;
}

struct Edge
{
    faceIndex f1, f2;
    int localEdge; // index of edge in f1 (0..2)
    Edge(faceIndex _f1, faceIndex _f2, int _e)
        : f1(_f1), f2(_f2), localEdge(_e) {}
};

bool Mesh::isFictitious(Point<float, 3> &p)
{
    return p[2] == std::numeric_limits<float>::infinity();
}

bool Mesh::isFictitious(Face<float, 3, 3> &f)
{
    return isFictitious(vertices[f.vertexAt(0)].getPosition()) || isFictitious(vertices[f.vertexAt(1)].getPosition()) || isFictitious(vertices[f.vertexAt(2)].getPosition());
}

void Mesh::makeDelaunay2D()
{
    // Transform a naive 2D triangulation into a Delaunay triangulation
    // using Lawson's edge flipping algorithm.

    std::queue<Edge> edgesToCheck;

    // Initialize a queue with all valid edges (non-fictitious, internal)
    for (faceIndex f = 0; f < faces.size(); ++f)
    {
        for (int i = 0; i < 3; ++i)
        {
            faceIndex neighbor = faces[f].facesNeighboorIDAt(i);
            if (neighbor == INVALID_FACE)
                continue;

            vertexIndex v0 = faces[f].vertexAt(i);
            vertexIndex v1 = faces[f].vertexAt((i + 1) % 3);

            Point<float, 3> &p0 = vertices[v0].getPosition();
            Point<float, 3> &p1 = vertices[v1].getPosition();

            // skip fictitious edge
            if (isFictitious(p0) || isFictitious(p1))
                continue;

            edgesToCheck.push(Edge(f, neighbor, i));
        }
    }

    while (!edgesToCheck.empty())
    {
        Edge e = edgesToCheck.front();
        edgesToCheck.pop();

        faceIndex f1 = e.f1;
        faceIndex f2 = e.f2;

        if (f1 == INVALID_FACE || f2 == INVALID_FACE)
            continue;

        // vertices of f1
        vertexIndex a = faces[f1].vertexAt((e.localEdge + 2) % 3);
        vertexIndex b = faces[f1].vertexAt(e.localEdge);
        vertexIndex c = faces[f1].vertexAt((e.localEdge + 1) % 3);

        // opposite vertex in f2
        vertexIndex u0 = faces[f2].vertexAt(0);
        vertexIndex u1 = faces[f2].vertexAt(1);
        vertexIndex u2 = faces[f2].vertexAt(2);

        vertexIndex edge0 = faces[f1].vertexAt((e.localEdge + 2) % 3);
        vertexIndex edge1 = faces[f1].vertexAt((e.localEdge + 1) % 3);

        vertexIndex d;
        if (u0 != edge0 && u0 != edge1)
            d = u0;
        else if (u1 != edge0 && u1 != edge1)
            d = u1;
        else
            d = u2;

        Point<float, 3> &pa = vertices[a].getPosition();
        Point<float, 3> &pb = vertices[b].getPosition();
        Point<float, 3> &pc = vertices[c].getPosition();
        Point<float, 3> &pd = vertices[d].getPosition();

        // skip if any vertex is fictitious
        if (isFictitious(pa) || isFictitious(pb) || isFictitious(pc) || isFictitious(pd))
            continue;

        // local Delaunay test
        float test = inCircle(pd, pa, pb, pc);

        if (test > 0.0f) // not locally Delaunay
        {
            flip(f1, f2);

            // enqueue affected edges for re-check
            for (int i = 0; i < 3; ++i)
            {
                faceIndex nf1 = faces[f1].facesNeighboorIDAt(i);
                faceIndex nf2 = faces[f2].facesNeighboorIDAt(i);

                vertexIndex va = faces[f1].vertexAt(i);
                vertexIndex vb = faces[f1].vertexAt((i + 1) % 3);
                Point<float, 3> &pva = vertices[va].getPosition();
                Point<float, 3> &pvb = vertices[vb].getPosition();

                // skip fictitious edges
                if (nf1 != INVALID_FACE && !isFictitious(pva) && !isFictitious(pvb))
                    edgesToCheck.push(Edge(f1, nf1, i));
                if (nf2 != INVALID_FACE && !isFictitious(pva) && !isFictitious(pvb))
                    edgesToCheck.push(Edge(f2, nf2, i));
            }
        }
    }
}

void Mesh::insertPoint(Point<float, 3> p)
{
    // create new vertex
    vertexIndex pIndex = static_cast<vertexIndex>(vertices.size());
    Vertex<float, 3> v;
    v.getPosition() = p;
    vertices.push_back(v);

    vertexIndex fictitiousVertex = INVALID_VERTEX;

    // find a triangle containing the point
    for (faceIndex fIndex = 0; fIndex < faces.size(); ++fIndex)
    {
        if (inTriangle(p, faces[fIndex]))
        {
            if (!isFictitious(faces[fIndex]))
            {
                // point is inside a real triangle -> split and done
                faceSplit(fIndex, pIndex);
                return;
            }
        }
        for(unsigned int i = 0; i < 3; ++i)
        {
            if(isFictitious(vertices[faces[fIndex].vertexAt(i)].getPosition()))
            {
                fictitiousVertex = faces[fIndex].vertexAt(i);
            }
        }
    }

    //circulator
    for(CirculatorOnFaces c = incidentFaces(fictitiousVertex); c.isValid(); ++c)
    {
        if (inTriangle(p, *c)) //probably not good, super vertex is in 3D
        {
            faceSplit(c.getCurrentFaceIndex(), pIndex);
            break;
        }
    }
    //flip other super faces 
    ///@todo
}

float Mesh::localSurfaceAround(vertexIndex v)
{
    float total = 0.0f;
    for (size_t fi = 0; fi < faces.size(); ++fi)
    {
        for (int k = 0; k < 3; ++k)
        {
            if (faces[fi].vertexAt(k) == v)
            {
                total += area(static_cast<faceIndex>(fi));
                break;
            }
        }
    }
    return total;
}

void Mesh::mergeCloseVertices(float threshold, float surfaceThreshold)
{
    float thresholdSq = threshold * threshold;
    const size_t n = vertices.size();
    std::vector<int> mapping(n);
    for (size_t i = 0; i < n; ++i)
        mapping[i] = static_cast<int>(i);

    for (size_t i = 0; i < n; ++i)
    {
        if (mapping[i] != static_cast<int>(i)) continue;
        Point<float, 3> &pi = vertices[i].getPosition();

        for (size_t j = i + 1; j < n; ++j)
        {
            if (mapping[j] != static_cast<int>(j)) continue;
            Point<float, 3> &pj = vertices[j].getPosition();

            float dx = pi[0] - pj[0];
            float dy = pi[1] - pj[1];
            float dz = pi[2] - pj[2];
            float distSq = dx*dx + dy*dy + dz*dz;

            if (distSq < thresholdSq)
            {
                float sBefore = 0.0f;
                for (size_t fi = 0; fi < faces.size(); ++fi)
                {
                    for (int k = 0; k < 3; ++k)
                        if (faces[fi].vertexAt(k) == i || faces[fi].vertexAt(k) == j)
                        {
                            sBefore += area(static_cast<faceIndex>(fi));
                            break;
                        }
                }

                mapping[j] = static_cast<int>(i);

                float sAfter = 0.0f;
                for (size_t fi = 0; fi < faces.size(); ++fi)
                {
                    for (int k = 0; k < 3; ++k)
                        if (faces[fi].vertexAt(k) == i || faces[fi].vertexAt(k) == j)
                        {
                            sAfter += area(static_cast<faceIndex>(fi));
                            break;
                        }
                }

                float diff = 0.0f;
                if (sBefore > 1e-12f)
                    diff = std::abs(sAfter - sBefore) / sBefore;

                if (diff > surfaceThreshold)
                {
                    mapping[j] = static_cast<int>(j); // Nope don't merge, lose too much surface
                    std::cout << "[mergeCloseVertices][WARN] (" << i << "," << j
                              << ") Local surface change = "
                              << diff * 100.0f << "%\n"
                              << "Did not merge" << std::endl;
                }
            }
        }
    }

    for (auto &f : faces)
        for (int k = 0; k < 3; ++k)
            f.vertexAt(k) = mapping[f.vertexAt(k)];

    std::vector<Vertex<float, 3>> newVertices;
    newVertices.reserve(vertices.size());
    std::unordered_map<int, int> newIndex;

    for (size_t i = 0; i < n; ++i)
        if (mapping[i] == static_cast<int>(i))
        {
            int ni = static_cast<int>(newVertices.size());
            newVertices.push_back(vertices[i]);
            newIndex[i] = ni;
        }

    for (auto &f : faces)
        for (int k = 0; k < 3; ++k)
            f.vertexAt(k) = newIndex[mapping[f.vertexAt(k)]];

    vertices = std::move(newVertices);
    std::cout << "[mergeCloseVertices] Merge done. Now at : "
              << vertices.size() << " vertices." << std::endl;
}


