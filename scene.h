#pragma once

#include <vector>
#include <string>
#include <limits>
#include <cstdio>
#include "geometry.h"
#include "material.h"
#include "lighting.h"
#include "ray.h"


const float EPSILON = 0.001f;
const float INF = std::numeric_limits<float>::infinity();

// Stores intersection details for shading
struct HitRecord {
    float t = INF;
    Vec3 point;
    Vec3 normal;
    Material material;
};

// sphere geometry with intersection logic
struct Sphere {
    Vec3 center;
    float radius;
    Material material;

    Sphere(Vec3 c, float r, Material m) : center(c), radius(r), material(m) {}

    bool intersect(const Ray& ray, float tMin, float tMax, HitRecord& hit) const {
        Vec3 oc = ray.origin - center;
        float a = ray.direction.dot(ray.direction);
        float b = 2.0f * oc.dot(ray.direction);
        float c = oc.dot(oc) - radius * radius;
        float disc = b * b - 4 * a * c;

        if (disc < 0) return false;
        float sqrtd = std::sqrt(disc);

        // closer root first
        float t = (-b - sqrtd) / (2 * a);
        if (t < tMin || t > tMax) {
            t = (-b + sqrtd) / (2 * a);
            if (t < tMin || t > tMax) return false;
        }

        hit.t = t;
        hit.point = ray.at(t);
        hit.normal = (hit.point - center).normalized();
        hit.material = material;
        return true;
    }
};

class Scene {
public:
    // camera setup
    Vec3 cameraPos = { 0,0,0 };
    Vec3 cameraFwd = { 0,0,1 };
    Vec3 cameraUp = { 0,1,0 };
    Vec3 cameraRight = { -1,0,0 };
    float fovY = 45;
    int width = 640, height = 480, maxDepth = 5;
    bool cameraFwdWasSet = false;

    // scene settings
    Color ambientLight = { 0,0,0 };
    Color backgroundColor = { 0,0,0 };
    std::string outputFile = "raytraced.bmp";

    // scene contents
    std::vector<Sphere> spheres;
    std::vector<PointLight> pointLights;
    std::vector<SpotLight> spotLights;

    // rendering methods
    void setupCamera();
    Ray generateRay(int i, int j) const;
    bool intersectScene(const Ray& ray, float tMin, float tMax, HitRecord& hit) const;
    bool inShadow(const Vec3& point, const Vec3& lightPos) const;

    Color shade(const HitRecord& hit, const Ray& ray, int depth) const;
    Color trace(const Ray& ray, int depth = 0) const;
};
