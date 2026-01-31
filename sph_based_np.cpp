#include <iostream>
#include <vector>
#include <raylib.h>
#include <cmath>


constexpr float EPS = 1e-6f;


float W_poly6(float radius, float h);
float gradW_spiky(float radius, float h);
float laplacianW_viscosity(float radius, float h);

float W_poly6(float r, float h) {
    if (r < 0.0f || r > h) return 0.0f;

    float h2 = h * h;
    float diff = h2 - r * r;
    float diff3 = diff * diff * diff;

    float coeff = 4.0f / (PI * pow(h, 8));
    return coeff * diff3;
}

float gradW_spiky(float r, float h) {
    if (r <= EPS || r > h) return 0.0f;

    float diff = h - r;
    float coeff = -30.0f / (PI * pow(h, 5));

    return coeff * diff * diff;
}

float laplacianW_viscosity(float r, float h) {
    if (r < 0.0f || r > h) return 0.0f;

    float coeff = 20.0f / (PI * pow(h, 5));
    return coeff * (h - r);
}


