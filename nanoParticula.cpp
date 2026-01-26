#include "raylib.h"
#include <cstdint>
#include <random>
#include <cmath>
#include <iostream>
#include <vector>


#define WIDTH 800
#define HEIGHT 600
#define NUM_PARTICLES 2000

#define SCALING_CONSTANT 1.5f
#define DAMPING_COEFF 0.99999f
#define GRAVITY 0.9f
#define RESTITUTION 0.3f
#define VELOCITY_EPSILON 0.01f
#define MU 0.4f

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




void HandleCollision(Particles2& particles, size_t curr, size_t other, float rad_sum, float dist) {

    // --- Contact normal
    float norm_x = particles.x_pos[curr] - particles.x_pos[other];
    float norm_y = particles.y_pos[curr] - particles.y_pos[other];

    if (dist < 1e-6f) {
        norm_x = 1.0f;
        norm_y = 0.0f;
        dist = 1.0f;
    }

    norm_x /= dist;
    norm_y /= dist;

    // --- Tangent
    float tang_x = -norm_y;
    float tang_y = norm_x;

    // --- Relative velocity
    float rel_vx = particles.v_x[curr] - particles.v_x[other];
    float rel_vy = particles.v_y[curr] - particles.v_y[other];

    float vn = rel_vx * norm_x + rel_vy * norm_y;

    // Separating contact
    if (vn > 0.0f)
        return;

    float invMass1 = 1.0f / particles.mass[curr];
    float invMass2 = 1.0f / particles.mass[other];

    // --- Normal impulse (restitution)
    float jn = -(1.0f + RESTITUTION) * vn;
    jn /= (invMass1 + invMass2);

    // --- Viscoelastic damping (normal direction)
    float kd = 0.05f;
    jn += -kd * vn;

    // Apply normal impulse
    particles.v_x[curr] += jn * norm_x * invMass1;
    particles.v_y[curr] += jn * norm_y * invMass1;
    particles.v_x[other] -= jn * norm_x * invMass2;
    particles.v_y[other] -= jn * norm_y * invMass2;

    // --- Recompute relative velocity AFTER normal impulse
    rel_vx = particles.v_x[curr] - particles.v_x[other];
    rel_vy = particles.v_y[curr] - particles.v_y[other];

    float vt = rel_vx * tang_x + rel_vy * tang_y;

    // --- Tangential (friction) impulse
    float jt = -vt / (invMass1 + invMass2);

    float maxFriction = MU * std::fabs(jn);

    // Coulomb friction clamp
    if (std::fabs(jt) > maxFriction)
        jt = (jt > 0.0f ? maxFriction : -maxFriction);

    // Static friction (sticking)
    if (std::fabs(vt) < 0.01f)
        jt = -vt / (invMass1 + invMass2);

    // Apply friction impulse
    particles.v_x[curr] += jt * tang_x * invMass1;
    particles.v_y[curr] += jt * tang_y * invMass1;
    particles.v_x[other] -= jt * tang_x * invMass2;
    particles.v_y[other] -= jt * tang_y * invMass2;

    float penetration = rad_sum - dist;
    if (penetration > 0.0f) {
        float percent = 0.8f;
        float slop = 0.01f;
        float corr = std::max(penetration - slop, 0.0f)
                    / (invMass1 + invMass2) * percent;

        particles.x_pos[curr] += corr * norm_x * invMass1;
        particles.y_pos[curr] += corr * norm_y * invMass1;
        particles.x_pos[other] -= corr * norm_x * invMass2;
        particles.y_pos[other] -= corr * norm_y * invMass2;
    }

}



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
                                HandleCollision(particles, curr, other, rad_sum, std::sqrt(dist_sq));
                            }
                        }
                    }
                }
            }
        }
    }
}

inline void ResolveWallContact(float& vx, float& vy,float norm_x, float norm_y) 
{
    // Tangent
    float tang_x = -norm_y;
    float tang_y =  norm_x;

    // Decompose velocity
    float vn = vx * norm_x + vy * norm_y;
    float vt = vx * tang_x + vy * tang_y;

    // Normal impulse
    float jn = -(1.0f + RESTITUTION) * vn;

    vx += jn * norm_x;
    vy += jn * norm_y;

    // Friction impulse
    float jt = -vt;
    float maxFriction = MU * std::fabs(jn);

    if (std::fabs(jt) > maxFriction)
        jt = (jt > 0.0f ? maxFriction : -maxFriction);

    vx += jt * tang_x;
    vy += jt * tang_y;
}


void UpdateParticle(Particles2& particles, int curr) {
    particles.x_pos[curr] += particles.v_x[curr];
    particles.y_pos[curr] += particles.v_y[curr];

    particles.v_y[curr] += GRAVITY;
    particles.v_x[curr] *= DAMPING_COEFF;
    particles.v_y[curr] *= DAMPING_COEFF;    

    float rad = particles.radius[curr];

    if (particles.x_pos[curr] - rad < 0.0f) {
        particles.x_pos[curr] = rad;
        ResolveWallContact(particles.v_x[curr], particles.v_y[curr], 1.0f, 0.0f);
    }
    if (particles.x_pos[curr] + rad > WIDTH) {
        particles.x_pos[curr] = WIDTH - rad;
        ResolveWallContact(particles.v_x[curr], particles.v_y[curr], -1.0f, 0.0f);
    }
    if (particles.y_pos[curr] + rad > HEIGHT) {
        particles.y_pos[curr] = HEIGHT - rad;
        ResolveWallContact(particles.v_x[curr], particles.v_y[curr], 0.0f, -1.0f);
    }
    if (particles.y_pos[curr] - rad < 0) {
        particles.y_pos[curr] = rad;
        ResolveWallContact(particles.v_x[curr], particles.v_y[curr], 0.0f, 1.0f);
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

            for (int it = 0; it < 5; it++) {
                ResetGrid(particles);
                CheckParticleCollisionGrid(particles);
            }

            UpdateParticles(particles);

            DrawParticles(particles);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}