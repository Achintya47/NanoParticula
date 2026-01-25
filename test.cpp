#include "raylib.h"
#include <cstdint>
#include <random>

#define WIDTH 800
#define HEIGHT 600
#define NUM_PARTICLES 100

typedef struct {
    float x_pos, y_pos;
    float v_x, v_y;
    float radius;
} Particle;

Particle particles[NUM_PARTICLES];

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