#include "raylib.h"
#include <cstdint>
#include <random>
#include <cmath>
#include <iostream>
#include <vector>


#define WIDTH 800
#define HEIGHT 600
#define NUM_PARTICLES 1000

#define SCALING_CONSTANT 1.5f
#define DAMPING_COEFF 0.99999f
#define GRAVITY 0.9f
#define RESTITUTION 0.3f
#define VELOCITY_EPSILON 0.01f

constexpr int CELL_SIZE = 50;
constexpr int GRID_WIDTH = WIDTH / CELL_SIZE;
constexpr int GRID_HEIGHT = HEIGHT / CELL_SIZE;
std::vector<std::vector<int>> grid(GRID_WIDTH * GRID_HEIGHT);

long long collisionChecks = 0; // Global Counter
long long totalChecks = 0; // Accumulate across frames
long long frameCount = 0; // Number of frames simulated
long long actualCollisions = 0;

constexpr int GRAPH_WIDTH = 200;
constexpr int GRAPH_HEIGHT = 100;
constexpr int HISTORY_SIZE = 200;

std::vector<int> collisionHistory(HISTORY_SIZE, 0);
int historyIndex = 0;


typedef struct {
    float x_pos, y_pos;
    float v_x, v_y;
    float radius;
    Color color;
    float mass;

} Particle;



struct Particles2 {
    std::vector<float> x_pos, y_pos;
    std::vector<float> v_x, v_y;
    std::vector<float> mass;
    std::vector<float> radius;

    Particles2(size_t n) {
        x_pos.resize(n);
        y_pos.resize(n);
        radius.resize(n);
        v_x.resize(n);
        v_y.resize(n);
        mass.resize(n);
    }
};


// AoS structure
Particle particles[NUM_PARTICLES];

Particle mouseParticle;

inline float dot(float x1, float y1, float x2, float y2) {
    return x1 * x2 + y1 * y2;
}

void ResetGrid(Particles2& particles){
    for (auto& cell : grid) cell.clear();

    for (int i = 0; i < NUM_PARTICLES; i++) {
        int cellX = particles.x_pos[i] / CELL_SIZE;
        int cellY =  particles.y_pos[i] / CELL_SIZE;

        if (cellX < 0) cellX = 0;
        if (cellX >= GRID_WIDTH) cellX = GRID_WIDTH - 1;
        if (cellY < 0) cellY = 0;
        if (cellY >= GRID_HEIGHT) cellY = GRID_HEIGHT - 1;

        int cellIndex = cellY * GRID_WIDTH + cellX;
        grid[cellIndex].push_back(i);
    }
}

void ReportCollisionStats() {
    frameCount++;
    totalChecks += collisionChecks;

    long long avgChecks = totalChecks / frameCount;

    static char buffer[128];
    snprintf(buffer, sizeof(buffer), "Avg collision checks/frame: %lld", avgChecks);
    // Draw on screen (top-left corner, below FPS)
    DrawText(buffer, 5, 25, 20, GREEN);

    // Also show current frame’s checks
    static char buffer2[128];
    snprintf(buffer2, sizeof(buffer2), "Collision checks this frame: %lld", collisionChecks);
    DrawText(buffer2, 5, 50, 20, YELLOW);

    // Also show current frame’s checks
    static char buffer3[128];
    snprintf(buffer2, sizeof(buffer3), "Actual Collisions this frame: %lld", actualCollisions);
    DrawText(buffer2, 5, 75, 20, BLUE);

    collisionChecks = 0; // reset for next frame
    actualCollisions = 0; // reset for next frame
}

// void HandleCollision(Particle* curr, Particle* other, float rad_sum, float dist_sq) {

//     // Unit normal Vector at collision
//     float norm_x = curr->x_pos - other->x_pos;
//     float norm_y = curr->y_pos - other->y_pos;
//     float magnitude = sqrt(norm_x * norm_x + norm_y * norm_y);
//     norm_x /= magnitude;
//     norm_y /= magnitude;

//     // Unit tangent Vector
//     float tang_x = -norm_y;
//     float tang_y = norm_x;

//     // Velocity component projection along tangent and normal
//     float v1_norm = dot(curr->v_x, curr->v_y, norm_x, norm_y);
//     float v2_norm = dot(other->v_x, other->v_y, norm_x, norm_y);
//     float v1_tang = dot(curr->v_x, curr->v_y, tang_x, tang_y);
//     float v2_tang = dot(other->v_x, other->v_y, tang_x, tang_y);

//     // Since mass is considered the same, velocities exchange, elastic collision
//     float temp = v1_norm;
//     v1_norm = ((v1_norm * (curr->mass - other->mass) + 2 * other->mass * v2_norm) / (curr->mass + other->mass)) * RESTITUTION;
//     v2_norm = ((v2_norm * (other->mass - curr->mass) + 2 * curr->mass * temp) / (curr->mass + other->mass)) * RESTITUTION;


//     // Reconstruct Velocities
//     curr->v_x = v1_norm * norm_x + v1_tang * tang_x;
//     curr->v_y = v1_norm * norm_y + v1_tang * tang_y;
//     other->v_x = v2_norm * norm_x + v2_tang * tang_x;
//     other->v_y = v2_norm * norm_y + v2_tang * tang_y;

//     // Seperate overlapping particles
//     float overlap = 0.5f * (rad_sum - sqrt(dist_sq));
//     curr->x_pos += overlap * norm_x;
//     curr->y_pos += overlap * norm_y;
//     other->x_pos -= overlap * norm_x;
//     other->y_pos -= overlap * norm_y;
// }


void HandleCollision(Particles2& particles, size_t curr, size_t other, float rad_sum, float dist_sq) {

    // Unit normal Vector at collision
    float norm_x = particles.x_pos[curr] - particles.x_pos[other];
    float norm_y = particles.y_pos[curr] - particles.y_pos[other];
    float magnitude = sqrt(norm_x * norm_x + norm_y * norm_y);

    if (magnitude == 0.0f) return;

    norm_x /= magnitude;
    norm_y /= magnitude;

    // Unit tangent Vector
    float tang_x = -norm_y;
    float tang_y = norm_x;

    // Velocity component projection along tangent and normal
    float v1_norm = dot(particles.v_x[curr], particles.v_y[curr], norm_x, norm_y);
    float v2_norm = dot(particles.v_x[other], particles.v_y[other], norm_x, norm_y);
    float v1_tang = dot(particles.v_x[curr], particles.v_y[curr], tang_x, tang_y);
    float v2_tang = dot(particles.v_x[other], particles.v_y[other], tang_x, tang_y);

    // Since mass is considered the same, velocities exchange, elastic collision
    float temp = v1_norm;
    v1_norm = ((v1_norm * (particles.mass[curr] - particles.mass[other]) + 2 * particles.mass[other] * v2_norm) / (particles.mass[curr] + particles.mass[other]));
    v2_norm = ((v2_norm * (particles.mass[other] - particles.mass[curr]) + 2 * particles.mass[curr] * temp) / (particles.mass[curr] + particles.mass[other]));


    // Reconstruct Velocities
    particles.v_x[curr] = v1_norm * norm_x + v1_tang * tang_x;
    particles.v_y[curr] = v1_norm * norm_y + v1_tang * tang_y;
    particles.v_x[other] = v2_norm * norm_x + v2_tang * tang_x;
    particles.v_y[other] = v2_norm * norm_y + v2_tang * tang_y;

    float overlap = (rad_sum - sqrt(dist_sq));
    float totalMass = particles.mass[curr] + particles.mass[other];
    float ratioCurr = particles.mass[other] / totalMass;
    float ratioOther = particles.mass[curr] / totalMass;

    // Apply mass-weighted correction
    particles.x_pos[curr] += overlap * ratioCurr * norm_x;
    particles.y_pos[curr] += overlap * ratioCurr * norm_y;
    particles.x_pos[other] -= overlap * ratioOther * norm_x;
    particles.y_pos[other] -= overlap * ratioOther * norm_y;

}


// void CheckParticleCollisionGrid() {
//     for (int cellY = 0; cellY < GRID_HEIGHT; cellY++) {
//         for (int cellX = 0; cellX < GRID_WIDTH; cellX++) {
//             int cellIndex = cellY * GRID_WIDTH + cellX;

//             // Check this cell + 8 neighbors
//             for (int ny = -1; ny <= 1; ny++) {
//                 for (int nx = -1; nx <= 1; nx++) {
//                     int nX = cellX + nx;
//                     int nY = cellY + ny;

//                     if (nX < 0 || nX >= GRID_WIDTH || nY < 0 || nY >= GRID_HEIGHT)
//                         continue;

                    

//                     int neighborIndex = nY * GRID_WIDTH + nX;

//                     // Compare particles in cell vs neighbor
//                     for (int i : grid[cellIndex]) {
//                         for (int j : grid[neighborIndex]) {
//                             if (i >= j) continue; // avoid double checks
//                             collisionChecks++;

//                             Particle* curr = &particles[i];
//                             Particle* other = &particles[j];

//                             float dx = curr->x_pos - other->x_pos;
//                             float dy = curr->y_pos - other->y_pos;
//                             float dist_sq = dx * dx + dy * dy;
//                             float rad_sum = curr->radius + other->radius;

//                             if (dist_sq <= rad_sum * rad_sum) {
//                                 actualCollisions++;
//                                 HandleCollision(curr, other, rad_sum, dist_sq);
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }

//  PASS GRID AS A PARAMETER LATER in Main
void CheckParticleCollisionGrid(Particles2& particles) {
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
                            collisionChecks++;

                            int curr = i;
                            int other = j;

                            float dx = particles.x_pos[curr] - particles.x_pos[other];
                            float dy = particles.y_pos[curr] - particles.y_pos[other];
                            float dist_sq = dx * dx + dy * dy;
                            float rad_sum = particles.radius[curr] + particles.radius[other];

                            if (dist_sq <= rad_sum * rad_sum) {
                                actualCollisions++;
                                HandleCollision(particles, curr, other, rad_sum, dist_sq);
                            }
                        }
                    }
                }
            }
        }
    }
}


// void UpdateParticle(Particle* particle) {
//     particle->x_pos += particle->v_x;
//     particle->y_pos += particle->v_y;

//     particle->v_x *= DAMPING_COEFF;
//     particle->v_y *= DAMPING_COEFF;

//     particle->v_y += GRAVITY;

//     if (fabs(particle->v_x) < 0.01f) particle->v_x = 0;
//     if (fabs(particle->v_y) < 0.01f) particle->v_y = 0;

//     float x_curr = particle->x_pos;
//     float y_curr = particle->y_pos;
//     float rad = particle->radius;

//     if (x_curr - rad < 0) {
//         particle->x_pos = rad;
//         particle->v_x = -particle->v_x * RESTITUTION;
//     }
//     if (x_curr + rad > WIDTH) {
//         particle->x_pos = WIDTH - rad;
//         particle->v_x = -particle->v_x * RESTITUTION;
//     }
//     if (y_curr + rad > HEIGHT) {
//         particle->y_pos = HEIGHT - rad;
//         particle->v_y = -particle->v_y * RESTITUTION;
//     }
//     if (y_curr - rad < 0) {
//         particle->y_pos = rad;
//         particle->v_y = -particle->v_y * RESTITUTION;
//     }
// }

void UpdateParticle(Particles2& particles, int curr) {
    particles.x_pos[curr] += particles.v_x[curr];
    particles.y_pos[curr] += particles.v_y[curr];

    particles.v_y[curr] += GRAVITY;
    particles.v_x[curr] *= DAMPING_COEFF;
    particles.v_y[curr] *= DAMPING_COEFF;    

    if (std::fabs(particles.v_x[curr]) < VELOCITY_EPSILON) particles.v_x[curr] = 0;
    if (std::fabs(particles.v_y[curr]) < VELOCITY_EPSILON) particles.v_y[curr] = 0;

    float rad = particles.radius[curr];

    if (particles.x_pos[curr] - rad < 0) {
        particles.x_pos[curr] = rad;
        particles.v_x[curr] = -particles.v_x[curr] * RESTITUTION;
    }
    if (particles.x_pos[curr] + rad > WIDTH) {
        particles.x_pos[curr] = WIDTH - rad;
        particles.v_x[curr] = -particles.v_x[curr] * RESTITUTION;
    }
    if (particles.x_pos[curr] + rad > HEIGHT) {
        particles.y_pos[curr] = HEIGHT - rad;
        particles.v_y[curr] = -particles.v_y[curr] * RESTITUTION;
    }
    if (particles.x_pos[curr] - rad < 0) {
        particles.y_pos[curr] = rad;
        particles.v_y[curr] = -particles.v_y[curr] * RESTITUTION;
    }
}

void DrawParticle(Particles2& particles, int curr) {
    DrawCircle(particles.x_pos[curr], particles.y_pos[curr], particles.radius[curr], BLUE);
}

void DrawParticles(Particles2& particles) {
    for (int i =0; i < NUM_PARTICLES; i++){
        DrawParticle(particles, i);
    }
}

void UpdateParticles(Particles2& particles) {
    for (int i =0; i < NUM_PARTICLES; i++){
        UpdateParticle(particles, i);
    }
}

void InitParticles(Particles2& particles) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> rad(1,2);
    std::uniform_int_distribution<> v_x(5, 10);
    std::uniform_int_distribution<> v_y(2,10);
    std::uniform_int_distribution<> x_pos(0, WIDTH);
    std::uniform_int_distribution<> y_pos(0, HEIGHT);

    // Mouse Particle
    mouseParticle.radius = 20; // bigger so it's visible
    mouseParticle.mass = mouseParticle.radius * mouseParticle.radius * SCALING_CONSTANT;
    mouseParticle.color = RED; // distinct color
    

    for (int i =0; i < NUM_PARTICLES; i++){
        particles.radius[i] = rad(gen);
        particles.v_x[i]= v_x(gen);
        particles.v_y[i] = v_y(gen);
        particles.x_pos[i] = x_pos(gen);
        particles.y_pos[i] = y_pos(gen);

        // particles[i].color = (Color){
        //     (unsigned char)GetRandomValue(50, 255),   // R
        //     (unsigned char)GetRandomValue(50, 255),   // G
        //     (unsigned char)GetRandomValue(50, 255),   // B
        //     255                                       // Alpha
        // };
        particles.mass[i] = particles.radius[i] * particles.radius[i] * SCALING_CONSTANT;
    }
}



int main() {
    InitWindow(800, 600, "Raylib + CMake on WSL");

    // To avoid running at full capacity, 
    SetTargetFPS(60);
    Particles2 particles(NUM_PARTICLES);

    InitParticles(particles);

    while (!WindowShouldClose()) {
        BeginDrawing();
            DrawFPS(5, 5);
            ClearBackground(BLACK);

            ReportCollisionStats();
            // DrawCollisionGraph(WIDTH - GRAPH_WIDTH - 10, 10);

            ResetGrid(particles);
            CheckParticleCollisionGrid(particles);
            // DrawCircle(mouseParticle.x_pos, mouseParticle.y_pos, mouseParticle.radius, mouseParticle.color);

            UpdateParticles(particles);
            // CheckParticleCollision();
            DrawParticles(particles);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}