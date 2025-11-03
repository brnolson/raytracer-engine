#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <chrono>
#include "scene.h"

float randomFloat() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dis(0.0f, 1.0f);
    return dis(gen);
}

void parseSceneFile(const std::string& filename, Scene& scene) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open " << filename << std::endl;
        return;
    }

    Material currentMat;
    std::vector<Vec3> verts;
    std::vector<Vec3> normals;
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
        else if (cmd == "samples_per_pixel")
            iss >> scene.samplesPerPixel;
        else if (cmd == "material") {
            iss >> currentMat.ambient.x >> currentMat.ambient.y >> currentMat.ambient.z
                >> currentMat.diffuse.x >> currentMat.diffuse.y >> currentMat.diffuse.z
                >> currentMat.specular.x >> currentMat.specular.y >> currentMat.specular.z
                >> currentMat.shininess
                >> currentMat.transmissive.x >> currentMat.transmissive.y >> currentMat.transmissive.z
                >> currentMat.ior;
        }
        else if (cmd == "max_vertices") {
            int maxV;
            iss >> maxV;
            verts.reserve(maxV);
        }
        else if (cmd == "max_normals") {
            int maxN;
            iss >> maxN;
            normals.reserve(maxN);
        }
        else if (cmd == "vertex") {
            float x, y, z;
            iss >> x >> y >> z;
            verts.push_back(Vec3(x, y, z));
        }
        else if (cmd == "normal") {
            float x, y, z;
            iss >> x >> y >> z;
            normals.push_back(Vec3(x, y, z));
        }
        else if (cmd == "sphere") {
            float x, y, z, r;
            iss >> x >> y >> z >> r;
            scene.spheres.push_back(Sphere(Vec3(x, y, z), r, currentMat));
        }
        else if (cmd == "triangle") {
            int v0, v1, v2;
            iss >> v0 >> v1 >> v2;
            
            // bounds check
            if (v0 < verts.size() && v1 < verts.size() && v2 < verts.size()) {
                scene.triangles.push_back(Triangle(verts[v0], verts[v1], verts[v2], currentMat));
            }
        }
        else if (cmd == "normal_triangle") {
            int v0, v1, v2, n0, n1, n2;
            iss >> v0 >> v1 >> v2 >> n0 >> n1 >> n2;
            
            if (v0 < verts.size() && v1 < verts.size() && v2 < verts.size() &&
                n0 < normals.size() && n1 < normals.size() && n2 < normals.size()) {
                scene.triangles.push_back(Triangle(
                    verts[v0], verts[v1], verts[v2],
                    normals[n0], normals[n1], normals[n2],
                    currentMat
                ));
            }
        }
        else if (cmd == "plane") {
            float px, py, pz, nx, ny, nz;
            iss >> px >> py >> pz >> nx >> ny >> nz;
            scene.planes.push_back(Plane(Vec3(px, py, pz), Vec3(nx, ny, nz), currentMat));
        }
        else if (cmd == "box") {
            float x0, y0, z0, x1, y1, z1;
            iss >> x0 >> y0 >> z0 >> x1 >> y1 >> z1;
            scene.boxes.push_back(Box(Vec3(x0, y0, z0), Vec3(x1, y1, z1), currentMat));
        }
        else if (cmd == "point_light") {
            float r, g, b, x, y, z;
            iss >> r >> g >> b >> x >> y >> z;
            scene.pointLights.push_back(PointLight(Color(r, g, b), Vec3(x, y, z)));
        }
        else if (cmd == "directional_light") {
            float r, g, b, dx, dy, dz;
            iss >> r >> g >> b >> dx >> dy >> dz;
            scene.directionalLights.push_back(DirectionalLight(Color(r, g, b), Vec3(dx, dy, dz)));
        }
        else if (cmd == "spot_light") {
            float r, g, b, px, py, pz, dx, dy, dz, a1, a2;
            iss >> r >> g >> b >> px >> py >> pz >> dx >> dy >> dz >> a1 >> a2;
            scene.spotLights.push_back(SpotLight(Color(r, g, b), Vec3(px, py, pz), Vec3(dx, dy, dz), a1, a2));
        }
    }

    scene.setupCamera();
    file.close();
    
    if (!scene.triangles.empty()) {
        scene.buildBVH();
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <scene.txt>\n";
        return 1;
    }

    Scene scene;
    parseSceneFile(argv[1], scene);

    std::cout << "Rendering " << scene.width << "x" << scene.height;
    if (scene.samplesPerPixel > 1) {
        std::cout << " (" << scene.samplesPerPixel << "spp)";  // forgot space
    }
    std::cout << " [" << scene.spheres.size() << " spheres, "
        << scene.triangles.size() << " tris, "
        << scene.planes.size() << " planes, "
        << scene.boxes.size() << " boxes]";
    if (scene.useBVH) {
        std::cout << " +BVH";
    }
    std::cout << "\nLights: " << scene.pointLights.size() << " pt, "
        << scene.directionalLights.size() << " dir, "
        << scene.spotLights.size() << " spot\n";

    std::vector<unsigned char> pixels(scene.width * scene.height * 3);
    auto t_start = std::chrono::high_resolution_clock::now();

    for (int j = 0; j < scene.height; ++j) {
        if (j % 50 == 0)
            std::cout << (int)(100.0f * j / scene.height) << "%\r" << std::flush;

        for (int i = 0; i < scene.width; ++i) {
            Color col(0, 0, 0);
            
            if (scene.samplesPerPixel > 1) {
                for (int s = 0; s < scene.samplesPerPixel; ++s) {
                    Ray r = scene.generateJitteredRay(i, j, randomFloat(), randomFloat());
                    col += scene.trace(r);
                }
                col = col / (float)scene.samplesPerPixel;
            } else {
                Ray r = scene.generateRay(i, j);
                col = scene.trace(r);
            }
            
            // clamp and convert
            col.x = std::min(1.0f, std::max(0.0f, col.x));
            col.y = std::min(1.0f, std::max(0.0f, col.y));
            col.z = std::min(1.0f, std::max(0.0f, col.z));

            int idx = (j * scene.width + i) * 3;
            pixels[idx] = (unsigned char)(col.x * 255);
            pixels[idx + 1] = (unsigned char)(col.y * 255);
            pixels[idx + 2] = (unsigned char)(col.z * 255);
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(t_end - t_start);
    
    std::cout << "100%\n";
    std::cout << "Took " << elapsed.count() << " seconds\n";

    // BMP stores pixels flipped horizontally for some reason
    for (int y = 0; y < scene.height; ++y) {
        for (int x = 0; x < scene.width / 2; ++x) {
            int left_idx = (y * scene.width + x) * 3;
            int right_idx = (y * scene.width + (scene.width - 1 - x)) * 3;

            std::swap(pixels[left_idx], pixels[right_idx]);
            std::swap(pixels[left_idx + 1], pixels[right_idx + 1]);
            std::swap(pixels[left_idx + 2], pixels[right_idx + 2]);
        }
    }

    std::string fname = scene.outputFile;

    if (fname.size() < 4 || fname.substr(fname.size() - 4) != ".bmp") {
        fname = fname.substr(0, fname.find_last_of('.')) + ".bmp";
    }

    if (stbi_write_bmp(fname.c_str(), scene.width, scene.height, 3, pixels.data()))
        std::cout << "Wrote " << fname << "\n";
    else
        std::cerr << "Failed to write image\n";

    return 0;
}