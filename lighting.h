#pragma once
#include "geometry.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct PointLight {
    Color intensity;
    Vec3 position;

    PointLight() = default;
    PointLight(const Color& intensity, const Vec3& position)
        : intensity(intensity), position(position) {}
};

struct DirectionalLight {
    Color intensity;
    Vec3 direction;

    DirectionalLight() = default;
    DirectionalLight(const Color& intensity, const Vec3& direction)
        : intensity(intensity), direction(direction.normalized()) {}
};

struct SpotLight {
    Color intensity;
    Vec3 position;
    Vec3 direction;
    float angle1;
    float angle2;

    SpotLight() = default;
    SpotLight(const Color& intensity, const Vec3& position, const Vec3& direction,
        float angle1, float angle2)
        : intensity(intensity), position(position), direction(direction.normalized()),
        angle1(angle1), angle2(angle2) {}

    float getFalloff(const Vec3& p) const {
        Vec3 toPoint = (p - position).normalized();
        float angle = std::acos(toPoint.dot(direction)) * 180.0f / M_PI;

        if (angle < angle1) return 1.0f;
        if (angle > angle2) return 0.0f;

        float t = (angle - angle1) / (angle2 - angle1);
        return 1.0f - t;
    }
};