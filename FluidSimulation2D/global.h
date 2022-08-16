#pragma once

extern constexpr int N = 100;
extern constexpr int iter = 8;

constexpr int Get2DCoordinate(int x, int y)
{
    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x > N - 1) x = N - 1;
    if (y > N - 1) y = N - 1;

    return (y * N) + x;
}