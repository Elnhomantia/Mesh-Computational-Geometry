#include <iostream>
#include <filesystem>

#include <mesh.h>

int main()
{
    std::filesystem::path cwd = std::filesystem::current_path();
    std::cout << "Current working directory: " << cwd << std::endl;

    Mesh m;
    m.loadOffFile("queen.off");
    m.check();
    // temperature directory need to exist in build directory
    m.writeTemperatureOFF("temperature/iteration", 20);
    m.writeCurvatureOFF("curvature.off");

    return EXIT_SUCCESS;
}