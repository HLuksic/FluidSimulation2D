#include <cstring>
#include "fluid.h"

Fluid::Fluid(float dt, float d, float v)
{
    diffusionAmt = d;
    timeStep     = dt;
    viscosity    = v;

    this->Reset();
}

void Fluid::Step()
{
    Diffuse(1, pVelocityX, velocityX, viscosity, timeStep);
    Diffuse(2, pVelocityY, velocityY, viscosity, timeStep);

    Project(pVelocityX, pVelocityY, velocityX, velocityY);

    Advect(1, velocityX, pVelocityX, pVelocityX, pVelocityY, timeStep);
    Advect(2, velocityY, pVelocityY, pVelocityX, pVelocityY, timeStep);

    Project(velocityX, velocityY, pVelocityX, pVelocityY);

    Diffuse(0, pDensity, density, diffusionAmt, timeStep);
    Advect(0, density, pDensity, velocityX, velocityY, timeStep);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            if (density[Get2DCoordinate(i, j)] > 0.0f)
                density[Get2DCoordinate(i, j)] -= 0.017f;
            if (density[Get2DCoordinate(i, j)] > 254.9f)
                density[Get2DCoordinate(i, j)] == 254.9f;
        }
}

void Fluid::Reset()
{
    memset(this->density,    0, sizeof(this->density));
    memset(this->pDensity,   0, sizeof(this->pDensity));
    memset(this->velocityX,  0, sizeof(this->velocityX));
    memset(this->velocityY,  0, sizeof(this->velocityY));
    memset(this->pVelocityX, 0, sizeof(this->pVelocityX));
    memset(this->pVelocityY, 0, sizeof(this->pVelocityY));
}

void Fluid::AddDensity(int x, int y, float amount)
{
    int index = Get2DCoordinate(x, y);
    density[index] += amount;
}

void Fluid::AddVelocity(int x, int y, float amountX, float amountY)
{
    int index = Get2DCoordinate(x, y);
    velocityX[index] += amountX;
    velocityY[index] += amountY;
}

void Fluid::Diffuse(int b, float x[], float x0[], float factor, float timeStep)
{
    float a = timeStep * factor * (N - 2) * (N - 2);
    this->LinearSolve(b, x, x0, a, 1 + 6 * a);
}

void Fluid::LinearSolve(int b, float x[], float x0[], float a, float c)
{
    float cRecip = 1.0f / c;

    for (int k = 0; k < iter; k++)
    {
        for (int j = 1; j < N - 1; j++)
            for (int i = 1; i < N - 1; i++)
                x[Get2DCoordinate(i, j)] = (x0[Get2DCoordinate(i, j)] + a * (x[Get2DCoordinate(i + 1, j)]
                                                                           + x[Get2DCoordinate(i - 1, j)]
                                                                           + x[Get2DCoordinate(i, j + 1)]
                                                                           + x[Get2DCoordinate(i, j - 1)]
                                                                           + x[Get2DCoordinate(i, j)]
                                                                           + x[Get2DCoordinate(i, j)])) * cRecip;
        this->SetBoundary(b, x);
    }
}

void Fluid::Project(float velocX[], float velocY[], float velocX0[], float velocY0[])
{
    for (int j = 1; j < N - 1; j++)
        for (int i = 1; i < N - 1; i++)
        {
            velocY0[Get2DCoordinate(i, j)] = -0.5f * (velocX[Get2DCoordinate(i + 1, j)]
                                                    - velocX[Get2DCoordinate(i - 1, j)]
                                                    + velocY[Get2DCoordinate(i, j + 1)]
                                                    - velocY[Get2DCoordinate(i, j - 1)]) / N;
            velocX0[Get2DCoordinate(i, j)] = 0;
        }

    this->SetBoundary(0, velocY0);
    this->SetBoundary(0, velocX0);
    this->LinearSolve(0, velocX0, velocY0, 1, 6);

    for (int j = 1; j < N - 1; j++)
        for (int i = 1; i < N - 1; i++)
        {
            velocX[Get2DCoordinate(i, j)] -= 0.5f * (velocX0[Get2DCoordinate(i + 1, j)] - velocX0[Get2DCoordinate(i - 1, j)]) * N;
            velocY[Get2DCoordinate(i, j)] -= 0.5f * (velocX0[Get2DCoordinate(i, j + 1)] - velocX0[Get2DCoordinate(i, j - 1)]) * N;
        }

    this->SetBoundary(1, velocX);
    this->SetBoundary(2, velocY);
}

void Fluid::Advect(int b, float d[], float d0[], float velocX[], float velocY[], float timeStep)
{
    int i0, i1, j0, j1;
    float s0, s1, t0, t1, x, y, timeStep0;

    timeStep0 = timeStep * (N - 2);

    for (int j = 1; j < N - 1; j++)
        for (int i = 1; i < N - 1; i++)
        {
            x = float(i) - timeStep0 * velocX[Get2DCoordinate(i, j)];
            y = float(j) - timeStep0 * velocY[Get2DCoordinate(i, j)];

            if (x < 0.5f) x = 0.5f;
            if (y < 0.5f) y = 0.5f;
            if (x > float(N) + 0.5f) x = float(N) + 0.5f;
            if (y > float(N) + 0.5f) y = float(N) + 0.5f;

            i0 = int(x);
            i1 = i0 + 1;

            j0 = int(y);
            j1 = j0 + 1;

            s1 = x - i0;
            s0 = 1.0f - s1;

            t1 = y - j0;
            t0 = 1.0f - t1;

            d[Get2DCoordinate(i, j)] =
                s0 * (t0 * d0[Get2DCoordinate(int(i0), int(j0))] + t1 * d0[Get2DCoordinate(int(i0), int(j1))]) +
                s1 * (t0 * d0[Get2DCoordinate(int(i1), int(j0))] + t1 * d0[Get2DCoordinate(int(i1), int(j1))]);
        }

    this->SetBoundary(b, d);
}

void Fluid::SetBoundary(int b, float x[])
{
    for (int i = 1; i < N - 1; i++)
    {
        x[Get2DCoordinate(i, 0)]     = b == 2 ? -x[Get2DCoordinate(i, 1)]     : x[Get2DCoordinate(i, 1)];
        x[Get2DCoordinate(i, N - 1)] = b == 2 ? -x[Get2DCoordinate(i, N - 2)] : x[Get2DCoordinate(i, N - 2)];
        x[Get2DCoordinate(0, i)]     = b == 1 ? -x[Get2DCoordinate(1, i)]     : x[Get2DCoordinate(1, i)];
        x[Get2DCoordinate(N - 1, i)] = b == 1 ? -x[Get2DCoordinate(N - 2, i)] : x[Get2DCoordinate(N - 2, i)];
    }

    x[Get2DCoordinate(0, 0)]         = 0.33f * (x[Get2DCoordinate(1, 0)]         + x[Get2DCoordinate(0, 1)]         + x[Get2DCoordinate(0, 0)]);
    x[Get2DCoordinate(0, N - 1)]     = 0.33f * (x[Get2DCoordinate(0, N - 2)]     + x[Get2DCoordinate(1, N - 1)]     + x[Get2DCoordinate(0, N - 1)]);
    x[Get2DCoordinate(N - 1, 0)]     = 0.33f * (x[Get2DCoordinate(N - 2, 0)]     + x[Get2DCoordinate(N - 1, 1)]     + x[Get2DCoordinate(N - 1, 0)]);
    x[Get2DCoordinate(N - 1, N - 1)] = 0.33f * (x[Get2DCoordinate(N - 1, N - 2)] + x[Get2DCoordinate(N - 2, N - 1)] + x[Get2DCoordinate(N - 1, N - 1)]);
}