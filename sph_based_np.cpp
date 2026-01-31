#include <iostream>
#include <vector>
#include <raylib.h>
#include <cmath>


constexpr float H = 16.0f;
constexpr float MASS = 1.0f;

constexpr float REST_DENSITY = 1000.0f;
constexpr float GAS_CONST   = 2000.0f;
constexpr float VISCOSITY   = 0.1f;

constexpr float DT = 0.001f;
constexpr float EPS = 1e-6f;

struct SPHsystem {
    int count;

    std::vector<Vector2> pos;
    std::vector<Vector2> vel;
    std::vector<Vector2> force;

    std::vector<float> density;
    std::vector<float> pressure;
};

/**
 * SECTION 1
 * Smoothing kernels used by the official SPH paper
 * implementation
 */
float W_poly6(float r, float h);
float gradW_spiky(float r, float h);
float laplacianW_viscosity(float r, float h);

/**
 * SECTION 2
 * Vector2 Helper functions
 */
float Vector2Length(Vector2 rij);
Vector2 Vector2Scale(Vector2 rij, float factor);
Vector2 Vector2Add(Vector2 v1, Vector2 v2);

/**
 * SECTION 3
 * Physics calculations
 */
void computeDensity(SPHsystem& sph);
void computePressure(SPHsystem& sph);
void computePressureFoce(SPHsystem& sph);
void computeViscosityFoce(SPHsystem& sph);
void integrate(SPHsystem& sph);


/* SECITON 1 */
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

/* SECTION 2 */
inline float Vector2Length(Vector2 rij) {
    return sqrtf(rij.x * rij.x + rij.y * rij.y);
}

inline Vector2 Vector2Scale(Vector2 rij, float fac) {
    Vector2 result;
    result.x = rij.x * fac;
    result.y = rij.y * fac;
    return result;
}

inline Vector2 Vector2Add(Vector2 v1, Vector2 v2) {
    Vector2 result;
    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    return result;
}

/* SECTION 3 */
void computeDensity(SPHsystem& sph) {
    for (int i = 0; i < sph.count; i++) {
        float rho = 0.0f;

        for (int j : neighbours(i)) {
            Vector2 rij = sph.pos[i] - sph.pos[j];
            float r = Vector2Length(rij);
            rho += MASS * W_poly6(r, H);
        }
        sph.density[i] = rho;
    }
}

void computePressure(SPHsystem& sph) {
    for (int i = 0; i < sph.count; i++) {
        sph.pressure[i] = 
            GAS_CONST * (sph.density[i] - REST_DENSITY);
    }
}
void computePressureForce(SPHsystem& sph) {
    for (int i = 0; i < sph.count; i++) {
        Vector2 f_pressure = {0, 0};

        for (int j : neighbours(i)) {
            if (i == j) continue;

            Vector2 rij = sph.pos[i] - sph.pos[j];
            float r = Vector2Length(rij);
            if (r <= EPS || r > H) continue;

            Vector2 dir = Vector2Scale(rij, 1.0f / r);
            float grad = gradW_spiky(r, H);

            float scalar = 
                -MASS * (sph.pressure[i] + sph.pressure[j]) / 
                (2.0f * sph.density[j]);
            
            f_pressure = Vector2Add(
                f_pressure,
                Vector2Scale(dir, scalar * grad)
            );
        }
        sph.force[i] = f_pressure;
    }
}

void integrate(SPHsystem& sph){
    for (int i = 0; i < sph.count; i++) {
        Vector2 accel = Vector2Scale(
            sph.force[i], 1.0f / sph.density[i]
        );

        sph.vel[i] = Vector2Add(
            sph.vel[i], Vector2Scale(accel, DT)
        );

        sph.pos[i] = Vector2Add(
            sph.pos[i], Vector2Scale(sph.vel[i], DT)
        );

        sph.force[i] = {0, 0};
    }
}
