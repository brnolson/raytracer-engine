#pragma once

#include <vector>
#include <string>
#include <limits>
#include <cstdio>
#include <algorithm>
#include <memory>
#include "geometry.h"
#include "material.h"
#include "lighting.h"
#include "ray.h"

const float EPSILON = 0.001f;
const float INF = std::numeric_limits<float>::infinity();

struct HitRecord {
    float t = INF;
    Vec3 point;
    Vec3 normal;
    Material material;
};

// bounding box for BVH
struct AABB {
    Vec3 min;
    Vec3 max;

    AABB() : min(Vec3(INF, INF, INF)), max(Vec3(-INF, -INF, -INF)) {}
    AABB(const Vec3& min, const Vec3& max) : min(min), max(max) {}

    bool intersect(const Ray& ray, float tMin, float tMax) const {
        for (int a = 0; a < 3; ++a) {
            float invD = 1.0f / (a == 0 ? ray.direction.x : (a == 1 ? ray.direction.y : ray.direction.z));
            float orig = a == 0 ? ray.origin.x : (a == 1 ? ray.origin.y : ray.origin.z);
            float minVal = a == 0 ? min.x : (a == 1 ? min.y : min.z);
            float maxVal = a == 0 ? max.x : (a == 1 ? max.y : max.z);
            
            float t0 = (minVal - orig) * invD;
            float t1 = (maxVal - orig) * invD;
            
            if (invD < 0.0f) std::swap(t0, t1);
            
            tMin = t0 > tMin ? t0 : tMin;
            tMax = t1 < tMax ? t1 : tMax;
            
            if (tMax <= tMin) return false;
        }
        return true;
    }

    static AABB merge(const AABB& a, const AABB& b) {
        Vec3 minPt(
            std::min(a.min.x, b.min.x),
            std::min(a.min.y, b.min.y),
            std::min(a.min.z, b.min.z)
        );
        Vec3 maxPt(
            std::max(a.max.x, b.max.x),
            std::max(a.max.y, b.max.y),
            std::max(a.max.z, b.max.z)
        );
        return AABB(minPt, maxPt);
    }

    Vec3 center() const {
        return (min + max) * 0.5f;
    }
};

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

    AABB getBoundingBox() const {
        Vec3 rad(radius, radius, radius);
        return AABB(center - rad, center + rad);
    }
};

struct Triangle {
    Vec3 v0, v1, v2;
    Vec3 n0, n1, n2;
    bool hasVertexNormals;
    Material material;

    Triangle(Vec3 a, Vec3 b, Vec3 c, Material m)
        : v0(a), v1(b), v2(c), material(m), hasVertexNormals(false) {}

    Triangle(Vec3 a, Vec3 b, Vec3 c, Vec3 na, Vec3 nb, Vec3 nc, Material m)
        : v0(a), v1(b), v2(c), n0(na), n1(nb), n2(nc), material(m), hasVertexNormals(true) {
        // sanity check normals
        if (na.length() < 0.001f || nb.length() < 0.001f || nc.length() < 0.001f) {
            hasVertexNormals = false;
        }
    }

    bool intersect(const Ray& ray, float tMin, float tMax, HitRecord& hit) const {
        Vec3 e1 = v1 - v0;
        Vec3 e2 = v2 - v0;
        Vec3 h = ray.direction.cross(e2);
        float a = e1.dot(h);

        if (std::abs(a) < EPSILON) return false;

        float f = 1.0f / a;
        Vec3 s = ray.origin - v0;
        float u = f * s.dot(h);

        if (u < 0.0f || u > 1.0f) return false;

        Vec3 q = s.cross(e1);
        float v = f * ray.direction.dot(q);

        if (v < 0.0f || u + v > 1.0f) return false;

        float t = f * e2.dot(q);

        if (t < tMin || t > tMax) return false;

        hit.t = t;
        hit.point = ray.at(t);
        
        if (hasVertexNormals) {
            float w = 1.0f - u - v;
            hit.normal = (n0 * w + n1 * u + n2 * v).normalized();
        } else {
            Vec3 faceNorm = e1.cross(e2).normalized();
            // flip if backfacing
            if (faceNorm.dot(ray.direction) > 0) {
                hit.normal = faceNorm * -1.0f;
            } else {
                hit.normal = faceNorm;
            }
        }

        hit.material = material;
        return true;
    }

    AABB getBoundingBox() const {
        Vec3 minPt(
            std::min(std::min(v0.x, v1.x), v2.x),
            std::min(std::min(v0.y, v1.y), v2.y),
            std::min(std::min(v0.z, v1.z), v2.z)
        );
        Vec3 maxPt(
            std::max(std::max(v0.x, v1.x), v2.x),
            std::max(std::max(v0.y, v1.y), v2.y),
            std::max(std::max(v0.z, v1.z), v2.z)
        );
        // pad a bit to avoid degenerate cases
        Vec3 pad(0.0001f, 0.0001f, 0.0001f);
        return AABB(minPt - pad, maxPt + pad);
    }
};

struct Plane {
    Vec3 point;
    Vec3 normal;
    Material material;

    Plane(Vec3 p, Vec3 n, Material m) : point(p), normal(n.normalized()), material(m) {}

    bool intersect(const Ray& ray, float tMin, float tMax, HitRecord& hit) const {
        float denom = normal.dot(ray.direction);
        
        if (std::abs(denom) < EPSILON) return false;

        float t = (point - ray.origin).dot(normal) / denom;

        if (t < tMin || t > tMax) return false;

        hit.t = t;
        hit.point = ray.at(t);
        hit.normal = denom < 0 ? normal : normal * -1.0f;
        hit.material = material;
        return true;
    }
};

struct Box {
    Vec3 minCorner;
    Vec3 maxCorner;
    Material material;

    Box(Vec3 minC, Vec3 maxC, Material m) : minCorner(minC), maxCorner(maxC), material(m) {}

    bool intersect(const Ray& ray, float tMin, float tMax, HitRecord& hit) const {
        float t1 = (minCorner.x - ray.origin.x) / ray.direction.x;
        float t2 = (maxCorner.x - ray.origin.x) / ray.direction.x;
        float t3 = (minCorner.y - ray.origin.y) / ray.direction.y;
        float t4 = (maxCorner.y - ray.origin.y) / ray.direction.y;
        float t5 = (minCorner.z - ray.origin.z) / ray.direction.z;
        float t6 = (maxCorner.z - ray.origin.z) / ray.direction.z;

        float tmin = std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
        float tmax = std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));

        if (tmax < 0 || tmin > tmax) return false;

        float t = tmin;
        if (t < tMin || t > tMax) return false;

        hit.t = t;
        hit.point = ray.at(t);
        hit.material = material;

        // figure out which face we hit
        Vec3 center = (minCorner + maxCorner) * 0.5f;
        Vec3 local = hit.point - center;
        Vec3 size = (maxCorner - minCorner) * 0.5f;
        
        float bias = 1.0001f;
        if (std::abs(local.x) > size.x * bias - EPSILON) {
            hit.normal = Vec3(local.x > 0 ? 1 : -1, 0, 0);
        } else if (std::abs(local.y) > size.y * bias - EPSILON) {
            hit.normal = Vec3(0, local.y > 0 ? 1 : -1, 0);
        } else {
            hit.normal = Vec3(0, 0, local.z > 0 ? 1 : -1);
        }

        return true;
    }

    AABB getBoundingBox() const {
        return AABB(minCorner, maxCorner);
    }
};

struct BVHNode {
    AABB box;
    std::unique_ptr<BVHNode> left;
    std::unique_ptr<BVHNode> right;
    int primStart;
    int primCount;

    BVHNode() : primStart(-1), primCount(0) {}

    bool isLeaf() const { return primCount > 0; }
};

class Scene {
public:
    Vec3 cameraPos = { 0,0,0 };
    Vec3 cameraFwd = { 0,0,1 };
    Vec3 cameraUp = { 0,1,0 };
    Vec3 cameraRight = { -1,0,0 };
    float fovY = 45;
    int width = 640, height = 480, maxDepth = 5;
    bool cameraFwdWasSet = false;
    bool useBVH = false;
    int samplesPerPixel = 1;

    Color ambientLight = { 0,0,0 };
    Color backgroundColor = { 0,0,0 };
    std::string outputFile = "raytraced.bmp";

    std::vector<Sphere> spheres;
    std::vector<Triangle> triangles;
    std::vector<Plane> planes;
    std::vector<Box> boxes;
    std::vector<PointLight> pointLights;
    std::vector<DirectionalLight> directionalLights;
    std::vector<SpotLight> spotLights;

    std::unique_ptr<BVHNode> bvhRoot;
    std::vector<int> triangleIndices;

    void setupCamera();
    void buildBVH();
    Ray generateRay(int i, int j) const;
    Ray generateJitteredRay(int i, int j, float jx, float jy) const;
    bool intersectScene(const Ray& ray, float tMin, float tMax, HitRecord& hit) const;
    bool intersectBVH(const Ray& ray, float tMin, float tMax, HitRecord& hit) const;
    bool intersectBVHNode(const BVHNode* node, const Ray& ray, float tMin, float& tMax, HitRecord& hit) const;
    bool inShadow(const Vec3& point, const Vec3& lightPos) const;
    bool inShadowDirectional(const Vec3& point, const Vec3& lightDir) const;

    Color shade(const HitRecord& hit, const Ray& ray, int depth) const;
    Color trace(const Ray& ray, int depth = 0) const;

private:
    std::unique_ptr<BVHNode> buildBVHRecursive(int start, int end, int depth);
    AABB computeBounds(int start, int end) const;
};