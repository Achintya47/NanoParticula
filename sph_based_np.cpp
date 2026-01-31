#include <iostream>
#include <vector>
#include <raylib.h>
#include <cmath>


constexpr float H = 16.0f;
constexpr float MASS = 64000.0f;

constexpr float REST_DENSITY = 1000.0f;
constexpr float GAS_CONST   = 1200.0f;
constexpr float VISCOSITY   = 0.3f;

constexpr float DT = 0.0008f;
constexpr float EPS = 1e-6f;

constexpr float MAX_SPEED = 500.0f;
constexpr float WALL_DAMPING = 0.5f;
constexpr float PARTICLE_RADIUS = 3.0f;
constexpr float VELOCITY_DAMPING = 0.99f;

constexpr Vector2 GRAVITY = {0.0f, 1500.0f};


#define SIM_WIDTH 800
#define SIM_HEIGHT 600

struct SPHsystem {
    int count;

    std::vector<Vector2> pos;
    std::vector<Vector2> vel;
    std::vector<Vector2> force;

    std::vector<float> density;
    std::vector<float> pressure;
};

struct NeighbourGrid {
    float cellSize;
    int gridWidth;
    int gridHeight;

    std::vector<std::vector<int>> cells;
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
inline float Vector2LengthSqr(Vector2 rij);
Vector2 Vector2Scale(Vector2 rij, float factor);
Vector2 Vector2Add(Vector2 v1, Vector2 v2);
Vector2 Vector2Subtract(Vector2 v1, Vector2 v2);

/**
 * SECTION 3
 * Physics calculations
 */
void computeDensity(SPHsystem& sph, const NeighbourGrid& grid);
void computePressure(SPHsystem& sph);
void computePressureForce(SPHsystem& sph, const NeighbourGrid& grid);
void computeViscosityFoce(SPHsystem& sph, const NeighbourGrid& grid);
void applyGravity(SPHsystem& sph);
void integrate(SPHsystem& sph);


/**
 * SECTION 4
 * Grid Helpers
 */
void initNeighbourGrid(NeighbourGrid& grid);
inline int cellIndex(int x, int y, int gridWidth);
void buildNeighbourGrid(NeighbourGrid& grid, 
    const std::vector<Vector2>& pos, int count);

template <typename Func>
void ForEachNeighbour(
    int i,
    const NeighbourGrid& grid, 
    const std::vector<Vector2>& pos,
    float h,
    Func&& func);

/**
 * SECTION 5
 * Raylib and Visualization part
 */
void drawParticles(const SPHsystem& sph);

/**
 * SECTION 6
 * Other Helpers
 */

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

inline float Vector2LengthSqr(Vector2 rij) {
    return (rij.x * rij.x + rij.y * rij.y);
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

inline Vector2 Vector2Subtract(Vector2 v1, Vector2 v2) {
    Vector2 result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    return result;
}

/* SECTION 3 */
void computeDensity(SPHsystem& sph, const NeighbourGrid& grid) {
    for (int i = 0; i < sph.count; i++) {
        float rho = 0.0f;

        ForEachNeighbour(i, grid, sph.pos, H, [&](int j) {
            Vector2 rij = Vector2Subtract(sph.pos[i], sph.pos[j]);
            float r = Vector2Length(rij);
            rho += MASS * W_poly6(r, H);
        });

        sph.density[i] = rho;
    }
}


void computePressure(SPHsystem& sph) {
    for (int i = 0; i < sph.count; i++) {
        float p = GAS_CONST * (sph.density[i] - REST_DENSITY);
        sph.pressure[i] = (p > 0.0f) ? p : 0.0f;
    }
}

void computePressureForce(SPHsystem& sph, const NeighbourGrid& grid) {
    for (int i = 0; i < sph.count; i++) {
        Vector2 f_pressure = {0, 0};

        ForEachNeighbour(i, grid, sph.pos, H, [&](int j) {
            if (i == j) return;
            Vector2 rij = Vector2Subtract(sph.pos[i], sph.pos[j]);
            float r2 = Vector2LengthSqr(rij);
            if (r2 <= EPS*EPS) return;

            float r = sqrtf(r2);
            Vector2 dir = Vector2Scale(rij, 1.0f / r);
            float grad = gradW_spiky(r, H);

            float scalar =
                -MASS * (sph.pressure[i] + sph.pressure[j]) /
                (2.0f * sph.density[j]);

            f_pressure = Vector2Add(
                f_pressure,
                Vector2Scale(dir, scalar * grad)
            );
        });

        sph.force[i] = f_pressure;
    }
}

void computeViscosityForce(SPHsystem& sph, const NeighbourGrid& grid) {
    for (int i = 0; i < sph.count; i++) {
        Vector2 f_visc = {0, 0};

        ForEachNeighbour(i, grid, sph.pos, H, [&](int j) {
            if (i == j) return;

            Vector2 rij = Vector2Subtract(sph.pos[i], sph.pos[j]);
            float r2 = Vector2LengthSqr(rij);
            if (r2 > H * H) return;

            float r = sqrtf(r2);
            float lap = laplacianW_viscosity(r, H);

            Vector2 velDiff =
                Vector2Subtract(sph.vel[j], sph.vel[i]);

            Vector2 term = Vector2Scale(
                velDiff,
                MASS * lap / sph.density[j]
            );

            f_visc = Vector2Add(f_visc, term);
        });

        // viscosity coefficient applied once
        sph.force[i] = Vector2Add(
            sph.force[i],
            Vector2Scale(f_visc, VISCOSITY)
        );
    }
}

void applyGravity(SPHsystem& sph) {
    for (int i = 0; i < sph.count; i++) {
        Vector2 g = Vector2Scale(GRAVITY, sph.density[i]);
        sph.force[i] = Vector2Add(sph.force[i], g);
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

        // Global velocity damping, the particles don't go WILD
        sph.vel[i] = Vector2Scale(
            sph.vel[i], VELOCITY_DAMPING
        );

        sph.pos[i] = Vector2Add(
            sph.pos[i], Vector2Scale(sph.vel[i], DT)
        );

        sph.force[i] = {0, 0};
    }
}

/**
 * SECTION 4
 */

void initNeighbourGrid(NeighbourGrid& grid) {
    grid.cellSize = H;

    grid.gridWidth  = (int)ceil(SIM_WIDTH  / grid.cellSize);
    grid.gridHeight = (int)ceil(SIM_HEIGHT / grid.cellSize);

    int totalCells = grid.gridWidth * grid.gridHeight;
    grid.cells.resize(totalCells);
}

inline int cellIndex(int x, int y, int gridWidth) {
    return y * gridWidth + x;
}

void buildNeighbourGrid(NeighbourGrid& grid,
    const std::vector<Vector2>& pos,
    int count) {

    // Clear old data
    for (auto& cell : grid.cells)
        cell.clear();

    // Insert particles
    for (int i = 0; i < count; i++) {
        int cx = (int)(pos[i].x / grid.cellSize);
        int cy = (int)(pos[i].y / grid.cellSize);

        // Clamp to grid
        cx = std::max(0, std::min(cx, grid.gridWidth  - 1));
        cy = std::max(0, std::min(cy, grid.gridHeight - 1));

        int idx = cellIndex(cx, cy, grid.gridWidth);
        grid.cells[idx].push_back(i);
    }
}

template <typename Func>
void ForEachNeighbour(
    int i,
    const NeighbourGrid& grid, 
    const std::vector<Vector2>& pos,
    float h,
    Func&& func) {

        // Find current particle i's cell index
        int cx = (int)(pos[i].x / grid.cellSize);
        int cy = (int)(pos[i].y / grid.cellSize);

        // Loop over nearby cells, typically a 3x3 loop
        // can be reduced from 9 iterations to 5-6 iterations
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {

                // Calculate cell indexes
                int nx = cx + dx;
                int ny = cy + dy;

                // Clamp to grid bounds
                if (nx < 0 || ny < 0 ||
                    nx >= grid.gridWidth ||
                    ny >= grid.gridHeight) continue;
                
                int cellIdx = cellIndex(nx, ny, grid.gridWidth);
                
                // Iterate over particles in that cell
                for (int j : grid.cells[cellIdx]) {
                    Vector2 rij = Vector2Subtract(pos[i], pos[j]);
                    // Only particles inside the kernel radius count
                    if (Vector2LengthSqr(rij) <= h * h)
                        func(j);
                }
            }
        }
}


/* SECTION 5 */
void drawParticles(const SPHsystem& sph) {
    constexpr float RADIUS = 3.0f;
    constexpr float MAX_SPEED = 500.0f;

    for (int i = 0; i < sph.count; i++) {
        float speed = Vector2Length(sph.vel[i]);

        float t = speed / MAX_SPEED;
        if (t > 1.0f) t = 1.0f;

        float hue = 240.f * (1.0f - t);
        Color col = ColorFromHSV(hue, 1.0f, 1.0f);

        DrawCircleV(sph.pos[i], RADIUS, col);
    }
}

void handleWallCollisions(SPHsystem& sph) {
    for (int i = 0; i < sph.count; i++) {

        // ---- LEFT WALL ----
        if (sph.pos[i].x < PARTICLE_RADIUS) {
            sph.pos[i].x = PARTICLE_RADIUS;
            sph.vel[i].x *= -WALL_DAMPING;
        }

        // ---- RIGHT WALL ----
        if (sph.pos[i].x > SIM_WIDTH - PARTICLE_RADIUS) {
            sph.pos[i].x = SIM_WIDTH - PARTICLE_RADIUS;
            sph.vel[i].x *= -WALL_DAMPING;
        }

        // ---- TOP WALL ----
        if (sph.pos[i].y < PARTICLE_RADIUS) {
            sph.pos[i].y = PARTICLE_RADIUS;
            sph.vel[i].y *= -WALL_DAMPING;
        }

        // ---- BOTTOM WALL ----
        if (sph.pos[i].y > SIM_HEIGHT - PARTICLE_RADIUS) {
            sph.pos[i].y = SIM_HEIGHT - PARTICLE_RADIUS;
            sph.vel[i].y *= -WALL_DAMPING;
        }
    }
}


/* SECTION 6 */
void initSPH(SPHsystem& sph, int nx, int ny) {
    sph.count = nx * ny;

    sph.pos.resize(sph.count);
    sph.vel.resize(sph.count);
    sph.force.resize(sph.count);
    sph.density.resize(sph.count);
    sph.pressure.resize(sph.count);

    float spacing = H * 0.5f;
    int index = 0;

    Vector2 start = {200.0f, 100.0f};

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            sph.pos[index] = {
                start.x + x * spacing,
                start.y + y * spacing
            };

            sph.vel[index] = {0.0f, 0.0f};
            sph.force[index] = {0.0f, 0.0f};
            sph.density[index] = REST_DENSITY;
            sph.pressure[index] = 0.0f;

            index++;
        }
    }
}


int main() {
    InitWindow(SIM_WIDTH, SIM_HEIGHT, "SPH Visualizer");
    SetTargetFPS(60);

    SPHsystem sph;
    NeighbourGrid grid;

    initSPH(sph, 40, 30);       // 1200 particles
    initNeighbourGrid(grid);

    while (!WindowShouldClose()) {

        // ---- PHYSICS ----
        buildNeighbourGrid(grid, sph.pos, sph.count);

        computeDensity(sph, grid);
        computePressure(sph);

        computePressureForce(sph, grid);
        computeViscosityForce(sph, grid);

        applyGravity(sph);

        integrate(sph);
        handleWallCollisions(sph);

        // ---- RENDER ----
        BeginDrawing();
        ClearBackground(BLACK);

        drawParticles(sph);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
