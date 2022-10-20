#pragma once

extern inline constexpr int gridSize = 100;
extern inline constexpr int iter = 4;

constexpr int Get2DCoordinate(int x, int y)
{
    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x > gridSize - 1) x = gridSize - 1;
    if (y > gridSize - 1) y = gridSize - 1;

    return (y * gridSize) + x;
}