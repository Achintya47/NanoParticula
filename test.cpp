#include "raylib.h"
#include <cstdint>

#define WIDTH 800
#define HEIGHT 600
int main() {
    InitWindow(800, 600, "Raylib + CMake on WSL");

    // To avoid running at full capacity, 
    SetTargetFPS(60);

    int16_t x, y, dir, incr;
    x = WIDTH/2;
    y = HEIGHT/2;
    dir = 1;
    incr = 2;
    
    while (!WindowShouldClose()) {
        BeginDrawing();
            ClearBackground(BLACK);

            DrawCircle(x, y, 40, WHITE);
            x += dir * incr;

            if (x >= WIDTH - 40 || x <= 40)
                dir *= -1;
        EndDrawing();
    }

    CloseWindow();
    return 0;
}