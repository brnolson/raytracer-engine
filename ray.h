#pragma once
#include "geometry.h"

// Represents a ray with origin and normalized direction
struct Ray {
    Vec3 origin;
    Vec3 direction;

    Ray() = default;

    // Constructs a ray and normalizes the direction
    Ray(const Vec3& o, const Vec3& d) : origin(o), direction(d.normalized()) {}

    // Computes a point along the ray at distance t
    Vec3 at(float t) const {
        return origin + direction * t;
    }
};