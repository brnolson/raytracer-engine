#include "scene.h"

void Scene::setupCamera() {
    // flip forward vector if it was explicitly set
    if (cameraFwdWasSet) cameraFwd = (-cameraFwd).normalized();
    else cameraFwd = cameraFwd.normalized();

    // recompute orthonormal basis
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

bool Scene::intersectScene(const Ray& ray, float tMin, float tMax, HitRecord& hit) const {
    bool hitAnything = false;
    float closest = tMax;
    for (const auto& sphere : spheres) {
        HitRecord temp;
        if (sphere.intersect(ray, tMin, closest, temp)) {
            hitAnything = true;
            closest = temp.t;
            hit = temp;
        }
    }
    return hitAnything;
}

bool Scene::inShadow(const Vec3& point, const Vec3& lightPos) const {
    Vec3 toLight = lightPos - point;
    float dist = toLight.length();
    Ray shadowRay(point, toLight);
    HitRecord temp;
    return intersectScene(shadowRay, EPSILON, dist - EPSILON, temp);
}

Color Scene::shade(const HitRecord& hit, const Ray& ray, int depth) const {
    Color colorAccum = ambientLight * hit.material.ambient;
    Vec3 viewDir = (cameraPos - hit.point).normalized();

    // point lights
    for (const auto& light : pointLights) {
        if (inShadow(hit.point, light.position)) continue;
        Vec3 L = (light.position - hit.point).normalized();
        float diff = std::max(0.0f, hit.normal.dot(L));
        Vec3 H = (L + viewDir).normalized();
        float spec = pow(std::max(0.0f, hit.normal.dot(H)), hit.material.shininess);
        float dist = (light.position - hit.point).length();
        float attenuation = 1.0f / (dist * dist + 0.1f);
        colorAccum += (hit.material.diffuse * diff + hit.material.specular * spec)
            * light.intensity * attenuation;
    }

	// spot lights
    for (const auto& light : spotLights) {
        if (inShadow(hit.point, light.position)) continue;
        float falloff = light.getFalloff(hit.point);
        if (falloff <= 0) continue;
        Vec3 L = (light.position - hit.point).normalized();
        float diff = std::max(0.0f, hit.normal.dot(L));
        Vec3 H = (L + viewDir).normalized();
        float spec = pow(std::max(0.0f, hit.normal.dot(H)), hit.material.shininess);
        float dist = (light.position - hit.point).length();
        float attenuation = 1.0f / (dist * dist + 0.1f);
        colorAccum += (hit.material.diffuse * diff + hit.material.specular * spec)
            * light.intensity * attenuation * falloff;
    }

    // reflection and refraction
    if (depth < maxDepth) {
        // reflection
        float refl = (hit.material.specular.x + hit.material.specular.y + hit.material.specular.z) / 3.0f;
        if (refl > 0.01f) {
            Vec3 R = ray.direction - hit.normal * (2.0f * ray.direction.dot(hit.normal));
            R = R.normalized();
            Ray reflRay(hit.point + R * EPSILON, R);
            Color reflected = trace(reflRay, depth + 1);
            colorAccum += reflected * hit.material.specular;
        }

        // refraction
        float transp = (hit.material.transmissive.x + hit.material.transmissive.y + hit.material.transmissive.z) / 3.0f;
        if (transp > 0.01f) {
            float eta = 1.0f / hit.material.ior;
            float cosI = -hit.normal.dot(ray.direction);
            float k = 1.0f - eta * eta * (1.0f - cosI * cosI);
            if (k >= 0) {
                Vec3 T = ray.direction * eta + hit.normal * (eta * cosI - sqrt(k));
                T = T.normalized();
                Ray refractRay(hit.point + T * EPSILON, T);
                Color transmitted = trace(refractRay, depth + 1);
                colorAccum += transmitted * hit.material.transmissive;
            }
        }
    }


    return colorAccum;
}

Color Scene::trace(const Ray& ray, int depth) const {
    if (depth >= maxDepth) return { 0, 0, 0 };
    HitRecord hit;
    if (intersectScene(ray, EPSILON, INF, hit))
        return shade(hit, ray, depth);
    return backgroundColor;
}
