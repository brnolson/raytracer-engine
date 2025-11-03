#pragma once
#include <cmath>
#include <algorithm>

struct Vec3 {
    float x, y, z;

    Vec3(float v = 0) : x(v), y(v), z(v) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3& v) const { return { x + v.x, y + v.y, z + v.z }; }
    Vec3 operator-(const Vec3& v) const { return { x - v.x, y - v.y, z - v.z }; }
    Vec3 operator*(float s) const { return { x * s, y * s, z * s }; }
    Vec3 operator/(float s) const { return { x / s, y / s, z / s }; }

    // for multiplying colors/vectors component-wise
    Vec3 operator*(const Vec3& v) const { return { x * v.x, y * v.y, z * v.z }; }

    Vec3& operator+=(const Vec3& v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    Vec3 operator-() const {
        return { -x, -y, -z };
    }



    float dot(const Vec3& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    Vec3 cross(const Vec3& v) const {
        return {
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        };
    }

    float length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vec3 normalized() const {
        float len = length();
        return len > 0 ? (*this) / len : *this;
    }
};

// keeps RGB values in valid range
inline Vec3 clampColor(const Vec3& c) {
    return {
        std::clamp(c.x, 0.0f, 1.0f),
        std::clamp(c.y, 0.0f, 1.0f),
        std::clamp(c.z, 0.0f, 1.0f)
    };
}

using Color = Vec3;