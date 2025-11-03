#include "scene.h"
#include <chrono>

void Scene::setupCamera() {
    if (cameraFwdWasSet) cameraFwd = (-cameraFwd).normalized();
    else cameraFwd = cameraFwd.normalized();
    
    cameraUp = cameraUp.normalized();
    cameraRight = cameraFwd.cross(cameraUp).normalized();
    cameraUp = cameraRight.cross(cameraFwd).normalized();

    printf("\nCamera: pos(%.2f,%.2f,%.2f) fov=%.1f\n",
        cameraPos.x, cameraPos.y, cameraPos.z, fovY);
}

Ray Scene::generateRay(int i, int j) const {
    float aspect = float(width) / height;
    float h = tan(fovY * M_PI / 180.0f);
    float w = aspect * h;
    float u = (2.0f * (i + 0.5f) / width - 1.0f) * w;
    float v = (1.0f - 2.0f * (j + 0.5f) / height) * h;
    Vec3 dir = (cameraFwd + cameraRight * u + cameraUp * v).normalized();
    return Ray(cameraPos, dir);
}

Ray Scene::generateJitteredRay(int i, int j, float jx, float jy) const {
    float aspect = float(width) / float(height);
    float halfHeight = tan(fovY * M_PI / 180.0f);
    float halfWidth = aspect * halfHeight;
    
    float u = (2.0f * (i + jx) / width - 1.0f) * halfWidth;
    float v = (1.0f - 2.0f * (j + jy) / height) * halfHeight;
    
    Vec3 dir = (cameraFwd + cameraRight * u + cameraUp * v).normalized();
    return Ray(cameraPos, dir);
}

AABB Scene::computeBounds(int start, int end) const {
    AABB result;
    for (int i = start; i < end; ++i) {
        int idx = triangleIndices[i];
        AABB box = triangles[idx].getBoundingBox();
        result = AABB::merge(result, box);
    }
    return result;
}

std::unique_ptr<BVHNode> Scene::buildBVHRecursive(int start, int end, int depth) {
    auto node = std::make_unique<BVHNode>();
    node->box = computeBounds(start, end);
    
    int numPrims = end - start;
    
    if (numPrims <= 4) {
        node->primStart = start;
        node->primCount = numPrims;
        return node;
    }
    
    // split on longest axis
    Vec3 extents = node->box.max - node->box.min;
    int splitAxis = 0;
    if (extents.y > extents.x) splitAxis = 1;
    if (extents.z > (splitAxis == 0 ? extents.x : extents.y)) splitAxis = 2;
    
    // sort by centroid along axis
    std::sort(triangleIndices.begin() + start, triangleIndices.begin() + end,
        [this, splitAxis](int a, int b) {
            Vec3 centerA = triangles[a].getBoundingBox().center();
            Vec3 centerB = triangles[b].getBoundingBox().center();
            float valA = (splitAxis == 0 ? centerA.x : (splitAxis == 1 ? centerA.y : centerA.z));
            float valB = (splitAxis == 0 ? centerB.x : (splitAxis == 1 ? centerB.y : centerB.z));
            return valA < valB;
        });
    
    int mid = start + numPrims / 2;
    node->left = buildBVHRecursive(start, mid, depth + 1);
    node->right = buildBVHRecursive(mid, end, depth + 1);
    
    return node;
}

void Scene::buildBVH() {
    if (triangles.empty()) {
        useBVH = false;
        printf("No triangles, skipping BVH build\n");
        return;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // initialize index array
    triangleIndices.resize(triangles.size());
    for (size_t i = 0; i < triangles.size(); ++i) {
        triangleIndices[i] = i;
    }
    
    bvhRoot = buildBVHRecursive(0, triangles.size(), 0);
    useBVH = true;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    printf("Built BVH for %zu triangles in %lld ms\n", triangles.size(), (long long)ms.count());
}

bool Scene::intersectBVHNode(const BVHNode* node, const Ray& ray, float tMin, float& tMax, HitRecord& hit) const {
    if (!node->box.intersect(ray, tMin, tMax)) {
        return false;
    }
    
    if (node->isLeaf()) {
        bool hitSomething = false;
        for (int i = 0; i < node->primCount; ++i) {
            int triIdx = triangleIndices[node->primStart + i];
            HitRecord temp;
            if (triangles[triIdx].intersect(ray, tMin, tMax, temp)) {
                hitSomething = true;
                tMax = temp.t;  // update for closer hits
                hit = temp;
            }
        }
        return hitSomething;
    }
    
    // recurse into children
    bool hitL = intersectBVHNode(node->left.get(), ray, tMin, tMax, hit);
    bool hitR = intersectBVHNode(node->right.get(), ray, tMin, tMax, hit);
    
    return hitL || hitR;
}

bool Scene::intersectBVH(const Ray& ray, float tMin, float tMax, HitRecord& hit) const {
    if (!bvhRoot) return false;
    return intersectBVHNode(bvhRoot.get(), ray, tMin, tMax, hit);
}

bool Scene::intersectScene(const Ray& ray, float tMin, float tMax, HitRecord& hit) const {
    bool foundHit = false;
    float closestT = tMax;
    
    // check spheres
    for (const auto& sphere : spheres) {
        HitRecord temp;
        if (sphere.intersect(ray, tMin, closestT, temp)) {
            foundHit = true;
            closestT = temp.t;
            hit = temp;
        }
    }

    // check triangles (via BVH if available)
    if (useBVH && bvhRoot) {
        HitRecord temp;
        if (intersectBVH(ray, tMin, closestT, temp)) {
            foundHit = true;
            closestT = temp.t;
            hit = temp;
        }
    } else {
        for (const auto& tri : triangles) {
            HitRecord temp;
            if (tri.intersect(ray, tMin, closestT, temp)) {
                foundHit = true;
                closestT = temp.t;
                hit = temp;
            }
        }
    }

    // check planes
    for (const auto& plane : planes) {
        HitRecord temp;
        if (plane.intersect(ray, tMin, closestT, temp)) {
            foundHit = true;
            closestT = temp.t;
            hit = temp;
        }
    }

    // check boxes
    for (const auto& box : boxes) {
        HitRecord temp;
        if (box.intersect(ray, tMin, closestT, temp)) {
            foundHit = true;
            closestT = temp.t;
            hit = temp;
        }
    }

    return foundHit;
}

bool Scene::inShadow(const Vec3& pt, const Vec3& lightPos) const {
    Vec3 toLight = lightPos - pt;
    float lightDist = toLight.length();
    Ray shadowRay(pt, toLight);
    HitRecord temp;
    return intersectScene(shadowRay, EPSILON, lightDist - EPSILON, temp);
}

bool Scene::inShadowDirectional(const Vec3& pt, const Vec3& lightDir) const {
    Ray shadowRay(pt, lightDir * -1.0f);
    HitRecord temp;
    return intersectScene(shadowRay, EPSILON, 10000.0f, temp);
}

Color Scene::shade(const HitRecord& hit, const Ray& ray, int depth) const {
    Color result = ambientLight * hit.material.ambient;
    Vec3 V = (cameraPos - hit.point).normalized();

    // point lights
    for (const auto& light : pointLights) {
        if (inShadow(hit.point, light.position)) continue;
        
        Vec3 L = (light.position - hit.point).normalized();
        float NdotL = std::max(0.0f, hit.normal.dot(L));
        
        Vec3 H = (L + V).normalized();
        float NdotH = std::max(0.0f, hit.normal.dot(H));
        float specular = pow(NdotH, hit.material.shininess);
        
        float d = (light.position - hit.point).length();
        float attenuation = 1.0f / (d * d + 0.1f);
        
        result += (hit.material.diffuse * NdotL + hit.material.specular * specular) 
                  * light.intensity * attenuation;
    }

    // directional lights
    for (const auto& light : directionalLights) {
        if (inShadowDirectional(hit.point, light.direction)) continue;
        
        Vec3 L = (light.direction * -1.0f).normalized();
        float NdotL = std::max(0.0f, hit.normal.dot(L));
        
        Vec3 H = (L + V).normalized();
        float NdotH = std::max(0.0f, hit.normal.dot(H));
        float specular = pow(NdotH, hit.material.shininess);
        
        result += (hit.material.diffuse * NdotL + hit.material.specular * specular) 
                  * light.intensity;
    }

    // spot lights
    for (const auto& light : spotLights) {
        if (inShadow(hit.point, light.position)) continue;
        
        float spotFalloff = light.getFalloff(hit.point);
        if (spotFalloff <= 0.0f) continue;
        
        Vec3 L = (light.position - hit.point).normalized();
        float NdotL = std::max(0.0f, hit.normal.dot(L));
        
        Vec3 H = (L + V).normalized();
        float NdotH = std::max(0.0f, hit.normal.dot(H));
        float specular = pow(NdotH, hit.material.shininess);
        
        float d = (light.position - hit.point).length();
        float attenuation = 1.0f / (d * d + 0.1f);
        
        result += (hit.material.diffuse * NdotL + hit.material.specular * specular) 
                  * light.intensity * attenuation * spotFalloff;
    }

    // reflections & refractions
    if (depth < maxDepth) {
        // reflection
        float avgSpec = (hit.material.specular.x + hit.material.specular.y + hit.material.specular.z) / 3.0f;
        if (avgSpec > 0.01f) {
            Vec3 reflected = ray.direction - hit.normal * (2.0f * ray.direction.dot(hit.normal));
            reflected = reflected.normalized();
            Ray reflRay(hit.point + reflected * EPSILON, reflected);
            Color reflColor = trace(reflRay, depth + 1);
            result += reflColor * hit.material.specular;
        }

        // refraction
        float avgTrans = (hit.material.transmissive.x + hit.material.transmissive.y + hit.material.transmissive.z) / 3.0f;
        if (avgTrans > 0.01f) {
            float eta = 1.0f / hit.material.ior;
            float cosTheta = -hit.normal.dot(ray.direction);
            float discriminant = 1.0f - eta * eta * (1.0f - cosTheta * cosTheta);
            
            if (discriminant >= 0.0f) {
                Vec3 refracted = ray.direction * eta + hit.normal * (eta * cosTheta - sqrt(discriminant));
                refracted = refracted.normalized();
                Ray refrRay(hit.point + refracted * EPSILON, refracted);
                Color refrColor = trace(refrRay, depth + 1);
                result += refrColor * hit.material.transmissive;
            }
        }
    }

    return result;
}

Color Scene::trace(const Ray& ray, int depth) const {
    if (depth >= maxDepth) return Color(0, 0, 0);
    
    HitRecord hit;
    if (intersectScene(ray, EPSILON, INF, hit)) {
        return shade(hit, ray, depth);
    }
    
    return backgroundColor;
}