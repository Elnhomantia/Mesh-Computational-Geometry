#include <iostream>
#include <filesystem>

#include <mesh.h>

int main()
{
    std::filesystem::path cwd = std::filesystem::current_path();
    std::cout << "Current working directory: " << cwd << std::endl;

    /*
    Mesh m;
    m.loadOffFile("queen.off");
    m.check();
    // temperature directory need to exist in build directory
    m.writeTemperatureOFF("temperature/iteration", 20);
    m.writeCurvatureOFF("curvature.off");
    */

    /*
    Mesh t;
    t.loadOffFile("simpleTestDelauney.off");
    //t.makeDelaunay2D();
    t.insertPoint(Point<float, 3>(0.6f, 0.7f, 0.0f));
    //t.insertPoint(Point<float, 3>(-1.0f, -1.0f, 0.0f));
    t.writeOffFile("simpleTestDelauneyOut.off");
    */

    Mesh merge;
    merge.loadOffFile("mergeTest.off");
    merge.mergeCloseVertices(0.4);
    merge.writeOffFile("mergeTestResult.off");

    Mesh mergeQueen;
    mergeQueen.loadOffFile("queen.off");
    mergeQueen.mergeCloseVertices(0.004);
    mergeQueen.writeOffFile("queenMerged.off");

    return EXIT_SUCCESS;
}