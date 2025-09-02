#include <iostream>
#include <filesystem>

#include"OFFparser.cpp"

int main()
{
    std::filesystem::path cwd = std::filesystem::current_path();
    std::cout << "Current working directory: " << cwd << std::endl;

    std::vector<Vertex<float, 3>> vertices;
    std::vector<Triangle<float, 3, 3>> faces;

    try
    {
        parseOFF<float, 3, 3>("../queen.off", vertices, faces);
    }
    catch(const std::runtime_error & e)
    {
        std::cout << "Error when parsing OFF file : " << e.what() << std::endl;
    }
    return 0;
}