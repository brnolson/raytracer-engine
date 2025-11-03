#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "scene.h"

// parser for scene files
void parseSceneFile(const std::string& filename, Scene& scene) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open " << filename << std::endl;
        return;
    }

    Material currentMaterial;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        size_t colon = line.find(':');
        if (colon == std::string::npos) continue;

        std::string cmd = line.substr(0, colon);
        std::string params = line.substr(colon + 1);
        std::istringstream iss(params);

        if (cmd == "camera_pos")
            iss >> scene.cameraPos.x >> scene.cameraPos.y >> scene.cameraPos.z;
        else if (cmd == "camera_fwd") {
            iss >> scene.cameraFwd.x >> scene.cameraFwd.y >> scene.cameraFwd.z;
            scene.cameraFwdWasSet = true;
        }
        else if (cmd == "camera_up")
            iss >> scene.cameraUp.x >> scene.cameraUp.y >> scene.cameraUp.z;
        else if (cmd == "camera_fov_ha")
            iss >> scene.fovY;
        else if (cmd == "film_resolution")
            iss >> scene.width >> scene.height;
        else if (cmd == "output_image")
            iss >> scene.outputFile;
        else if (cmd == "background")
            iss >> scene.backgroundColor.x >> scene.backgroundColor.y >> scene.backgroundColor.z;
        else if (cmd == "ambient_light")
            iss >> scene.ambientLight.x >> scene.ambientLight.y >> scene.ambientLight.z;
        else if (cmd == "max_depth")
            iss >> scene.maxDepth;
        else if (cmd == "material") {
            iss >> currentMaterial.ambient.x >> currentMaterial.ambient.y >> currentMaterial.ambient.z
                >> currentMaterial.diffuse.x >> currentMaterial.diffuse.y >> currentMaterial.diffuse.z
                >> currentMaterial.specular.x >> currentMaterial.specular.y >> currentMaterial.specular.z
                >> currentMaterial.shininess
                >> currentMaterial.transmissive.x >> currentMaterial.transmissive.y >> currentMaterial.transmissive.z
                >> currentMaterial.ior;
        }
        else if (cmd == "sphere") {
            float x, y, z, r;
            iss >> x >> y >> z >> r;
            scene.spheres.emplace_back(Vec3(x, y, z), r, currentMaterial);
        }
        else if (cmd == "point_light") {
            float r, g, b, x, y, z;
            iss >> r >> g >> b >> x >> y >> z;
            scene.pointLights.emplace_back(Color(r, g, b), Vec3(x, y, z));
        }
        else if (cmd == "spot_light") {
            float r, g, b, px, py, pz, dx, dy, dz, a1, a2;
            iss >> r >> g >> b >> px >> py >> pz >> dx >> dy >> dz >> a1 >> a2;
            scene.spotLights.emplace_back(Color(r, g, b), Vec3(px, py, pz), Vec3(dx, dy, dz), a1, a2);
        }
    }

    scene.setupCamera();
    file.close();
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <scene.txt>\n";
        return 1;
    }

    Scene scene;
    parseSceneFile(argv[1], scene);

    std::cout << "Rendering " << scene.width << "x" << scene.height
        << " (" << scene.spheres.size() << " spheres, "
        << scene.pointLights.size() << " point lights, "
        << scene.spotLights.size() << " spot lights)\n";

    std::vector<unsigned char> pixels(scene.width * scene.height * 3);

    // render pixels
    for (int j = 0; j < scene.height; ++j) {
        if (j % 50 == 0)
            std::cout << (100.0f * j / scene.height) << "%\r" << std::flush;

        for (int i = 0; i < scene.width; ++i) {
            Ray ray = scene.generateRay(i, j);
            Color color = clampColor(scene.trace(ray));

            int idx = (j * scene.width + i) * 3;
            pixels[idx + 0] = (unsigned char)(color.x * 255);
            pixels[idx + 1] = (unsigned char)(color.y * 255);
            pixels[idx + 2] = (unsigned char)(color.z * 255);
        }
    }

    std::cout << "100%\n";

    // flip image horizontally
    for (int y = 0; y < scene.height; ++y) {
        for (int x = 0; x < scene.width / 2; ++x) {
            int left = (y * scene.width + x) * 3;
            int right = (y * scene.width + (scene.width - 1 - x)) * 3;

            std::swap(pixels[left + 0], pixels[right + 0]);
            std::swap(pixels[left + 1], pixels[right + 1]);
            std::swap(pixels[left + 2], pixels[right + 2]);
        }
    }

    std::string filename = scene.outputFile;

    // file to .bmp
    if (filename.size() < 4 || filename.substr(filename.size() - 4) != ".bmp") {
        filename = filename.substr(0, filename.find_last_of('.')) + ".bmp";
    }

    if (stbi_write_bmp(filename.c_str(), scene.width, scene.height, 3, pixels.data()))
        std::cout << "Saved: " << filename << "\n";
    else
        std::cerr << "Error saving image!\n";

    return 0;
}