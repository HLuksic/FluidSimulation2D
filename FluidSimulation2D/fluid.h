#pragma once

#include "global.h"

class Fluid
{
public:
    Fluid(float dt, float d, float v);

public:
    float density[N * N];
private:
    float pDensity[N * N]; //previous density
    float velocityX[N * N];
    float velocityY[N * N];
    float pVelocityX[N * N]; // previous velocity
    float pVelocityY[N * N];
    float timeStep;
    float diffusionAmt;
    float viscosity; //thickness

public:
    void Step();
    void Reset();
    void AddDensity(int x, int y, float amount);
    void AddVelocity(int x, int y, float amountX, float amountY);
    void Diffuse(int b, float x[], float x0[], float factor, float timeStep);
    void LinearSolve(int b, float x[], float x0[], float a, float c);
    void Project(float velocX[], float velocY[], float velocX0[], float velocY0[]);
    void Advect(int b, float d[], float d0[], float velocX[], float velocY[], float timeStep);
    void SetBoundary(int b, float x[]);
};