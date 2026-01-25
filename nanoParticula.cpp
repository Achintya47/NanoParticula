#include "raylib.h"
#include <cstdint>
#include <random>
#include <cmath>
#include <iostream>
#include <vector>


#define WIDTH 800
#define HEIGHT 600
#define NUM_PARTICLES 100
#define SCALING_CONSTANT 1.5f
#define DAMPING_COEFF 0.99f
#define GRAVITY 0.3f
#define RESTITUTION 0.8f

constexpr int CELL_SIZE = 50;
constexpr int GRID_WIDTH = WIDTH / CELL_SIZE;
constexpr int GRID_HEIGHT = HEIGHT / CELL_SIZE;
std::vector<std::vector<int>> grid(GRID_WIDTH * GRID_HEIGHT);

typedef struct {
    float x_pos, y_pos;
    float v_x, v_y;
    float radius;
    Color color;
    float mass;

} Particle;

// AoS structure
Particle particles[NUM_PARTICLES];

inline float dot(float x1, float y1, float x2, float y2) {
    return x1 * x2 + y1 * y2;
}

void ResetGrid(){
    for (auto& cell : grid) cell.clear();

    for (int i = 0; i < NUM_PARTICLES; i++) {
        int cellX = particles[i].x_pos / CELL_SIZE;
        int cellY = particles[i].y_pos / CELL_SIZE;

        if (cellX < 0) cellX = 0;
        if (cellX >= GRID_WIDTH) cellX = GRID_WIDTH - 1;
        if (cellY < 0) cellY = 0;
        if (cellY >= GRID_HEIGHT) cellY = GRID_HEIGHT - 1;

        int cellIndex = cellY * GRID_WIDTH + cellX;
        grid[cellIndex].push_back(i);
    }
}

void CheckParticleCollisionGrid() {
    for (int cellY = 0; cellY < GRID_HEIGHT; cellY++) {
        for (int cellX = 0; cellX < GRID_WIDTH; cellX++) {
            int cellIndex = cellY * GRID_WIDTH + cellX;

            // Check this cell + 8 neighbors
            for (int ny = -1; ny <= 1; ny++) {
                for (int nx = -1; nx <= 1; nx++) {
                    int nX = cellX + nx;
                    int nY = cellY + ny;

                    if (nX < 0 || nX >= GRID_WIDTH || nY < 0 || nY >= GRID_HEIGHT)
                        continue;

                    int neighborIndex = nY * GRID_WIDTH + nX;

                    // Compare particles in cell vs neighbor
                    for (int i : grid[cellIndex]) {
                        for (int j : grid[neighborIndex]) {
                            if (i >= j) continue; // avoid double checks

                            Particle* curr = &particles[i];
                            Particle* other = &particles[j];

                            float dx = curr->x_pos - other->x_pos;
                            float dy = curr->y_pos - other->y_pos;
                            float dist_sq = dx * dx + dy * dy;
                            float rad_sum = curr->radius + other->radius;

                            if (dist_sq <= rad_sum * rad_sum) {
                                // Handle collision (same as your CheckParticleCollision)
                                // Normalization, velocity exchange, restitution, overlap correction...
                            }
                        }
                    }
                }
            }
        }
    }
}


void UpdateParticle(Particle* particle) {
    particle->x_pos += particle->v_x;
    particle->y_pos += particle->v_y;

    particle->v_x *= DAMPING_COEFF;
    particle->v_y *= DAMPING_COEFF;

    particle->v_y += GRAVITY;

    if (fabs(particle->v_x) < 0.01f) particle->v_x = 0;
    if (fabs(particle->v_y) < 0.01f) particle->v_y = 0;

    float x_curr = particle->x_pos;
    float y_curr = particle->y_pos;
    float rad = particle->radius;

    if (x_curr - rad < 0) {
        particle->x_pos = rad;
        particle->v_x = -particle->v_x * RESTITUTION;
    }
    if (x_curr + rad > WIDTH) {
        particle->x_pos = WIDTH - rad;
        particle->v_x = -particle->v_x * RESTITUTION;
    }
    if (y_curr + rad > HEIGHT) {
        particle->y_pos = HEIGHT - rad;
        particle->v_y = -particle->v_y * RESTITUTION;
    }
    if (y_curr - rad < 0) {
        particle->y_pos = rad;
        particle->v_y = -particle->v_y * RESTITUTION;
    }
}

void DrawParticle(Particle* particle) {
    DrawCircle(particle->x_pos, particle->y_pos, particle->radius, particle->color);
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

        particles[i].color = (Color){
            (unsigned char)GetRandomValue(50, 255),   // R
            (unsigned char)GetRandomValue(50, 255),   // G
            (unsigned char)GetRandomValue(50, 255),   // B
            255                                       // Alpha
        };
        particles[i].mass = particles[i].radius * particles[i].radius * SCALING_CONSTANT;
    }
}


void HandleCollision(Particle* curr, Particle* other, float rad_sum, float dist_sq) {

    // Unit normal Vector at collision
    float norm_x = curr->x_pos - other->x_pos;
    float norm_y = curr->y_pos - other->y_pos;
    float magnitude = sqrt(norm_x * norm_x + norm_y * norm_y);
    norm_x /= magnitude;
    norm_y /= magnitude;

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
    v1_norm = ((v1_norm * (curr->mass - other->mass) + 2 * other->mass * v2_norm) / (curr->mass + other->mass)) * RESTITUTION;
    v2_norm = ((v2_norm * (other->mass - curr->mass) + 2 * curr->mass * temp) / (curr->mass + other->mass)) * RESTITUTION;


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


void CheckParticleCollision() {
    Particle* curr;
    Particle* other;

    for (int i = 0; i < NUM_PARTICLES; i++){
        curr = particles + i;
        for (int j = i + 1; j < NUM_PARTICLES; j++){
            if (j == i)
                continue;
            other = particles + j;
            float dx = curr->x_pos - other->x_pos;
            float dy = curr->y_pos - other->y_pos;
            float dist_sq = dx * dx + dy * dy;
            float rad_sum = curr->radius + other->radius;
            
            if (dist_sq <= rad_sum * rad_sum) {
                HandleCollision(curr, other, rad_sum, dist_sq);
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