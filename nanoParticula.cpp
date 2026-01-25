#include "raylib.h"
#include <cstdint>
#include <random>
#include <cmath>

#define WIDTH 800
#define HEIGHT 600
#define NUM_PARTICLES 10

typedef struct {
    float x_pos, y_pos;
    float v_x, v_y;
    float radius;
} Particle;

Particle particles[NUM_PARTICLES];

float dot(float x1, float y1, float x2, float y2) {
    return x1 * x2 + y1 * y2;
}

void UpdateParticle(Particle* particle) {
    particle->x_pos += particle->v_x;
    particle->y_pos += particle->v_y;

    float x_curr = particle->x_pos;
    float y_curr = particle->y_pos;
    float rad = particle->radius;

    if (x_curr - rad < 0) {
        particle->x_pos = rad;
        particle->v_x = -particle->v_x;
    }
    if (x_curr + rad > WIDTH) {
        particle->x_pos = WIDTH - rad;
        particle->v_x = -particle->v_x;
    }
    if (y_curr + rad > HEIGHT) {
        particle->y_pos = HEIGHT - rad;
        particle->v_y = -particle->v_y;
    }
    if (y_curr - rad < 0) {
        particle->y_pos = rad;
        particle->v_y = -particle->v_y;
    }
}

void DrawParticle(Particle* particle) {
    DrawCircle(particle->x_pos, particle->y_pos, particle->radius, WHITE);
}

void DrawParticles() {
    for (int i =0; i < NUM_PARTICLES; i++){
        DrawParticle(particles + i);
    }
}

void UpdateParticles() {
    for (int i =0; i < NUM_PARTICLES; i++){
        UpdateParticle(particles + i);
    }
}

void InitParticles() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> rad(5, 15);
    std::uniform_int_distribution<> v_x(5, 10);
    std::uniform_int_distribution<> v_y(2,10);
    std::uniform_int_distribution<> x_pos(0, WIDTH);
    std::uniform_int_distribution<> y_pos(0, HEIGHT);

    for (int i =0; i < NUM_PARTICLES; i++){
        particles[i].radius = rad(gen);
        particles[i].v_x = v_x(gen);
        particles[i].v_y = v_y(gen);
        particles[i].x_pos = x_pos(gen);
        particles[i].y_pos = y_pos(gen);
    }
}

void CheckParticleCollision() {
    Particle* curr;
    Particle* other;

    for (int i = 0; i < NUM_PARTICLES; i++){
        curr = particles + i;
        for (int j = 0; j < NUM_PARTICLES; j++){
            if (j == i)
                continue;
            other = particles + j;
            float dx = curr->x_pos - other->x_pos;
            float dy = curr->y_pos - other->y_pos;
            float dist_sq = dx * dx + dy * dy;
            float rad_sum = curr->radius + other->radius;
            
            if (dist_sq <= rad_sum * rad_sum) {

                // Unit normal Vector at collision
                float norm_x = curr->x_pos - other->x_pos;
                float norm_y = curr->y_pos - other->y_pos;
                norm_x /= sqrt(norm_x * norm_x + norm_y * norm_y);
                norm_y /= sqrt(norm_x * norm_x + norm_y * norm_y);

                // Unit tangent Vector
                float tang_x = -norm_y;
                float tang_y = norm_x;

                // Velocity component projection along tangent and normal
                float v1_norm = dot(curr->v_x, curr->v_y, norm_x, norm_y);
                float v2_norm = dot(other->v_x, other->v_y, norm_x, norm_y);
                float v1_tang = dot(curr->v_x, curr->v_y, tang_x, tang_y);
                float v2_tang = dot(other->v_x, other->v_y, tang_x, tang_y);

                // Since mass is considered the same, velocities exchange, elastic collision
                float temp = v1_norm;
                v1_norm = v2_norm;
                v2_norm = temp;

                // Reconstruct Velocities
                curr->v_x = v1_norm * norm_x + v1_tang * tang_x;
                curr->v_y = v1_norm * norm_y + v1_tang * tang_y;
                other->v_x = v2_norm * norm_x + v2_tang * tang_x;
                other->v_y = v2_norm * norm_y + v2_tang * tang_y;

                // Seperate overlapping particles
                float overlap = 0.5f * (rad_sum - sqrt(dist_sq));
                curr->x_pos += overlap * norm_x;
                curr->y_pos += overlap * norm_y;
                other->x_pos -= overlap * norm_x;
                other->y_pos -= overlap * norm_y;

            }
        }
    }
}


int main() {
    InitWindow(800, 600, "Raylib + CMake on WSL");

    // To avoid running at full capacity, 
    SetTargetFPS(60);

    InitParticles();

    while (!WindowShouldClose()) {
        BeginDrawing();
            DrawFPS(5, 5);
            ClearBackground(BLACK);
            UpdateParticles();
            CheckParticleCollision();
            DrawParticles();
        EndDrawing();
    }

    CloseWindow();
    return 0;
}