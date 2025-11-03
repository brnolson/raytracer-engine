#pragma once
#include "geometry.h"

struct Material {
    Color ambient = { 0, 0, 0 };
    Color diffuse = { 1, 1, 1 };
    Color specular = { 0, 0, 0 };
    Color transmissive = { 0, 0, 0 };
    float shininess = 5.0f;
    float ior = 1.0f;  // index of refraction
};