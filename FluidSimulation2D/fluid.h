#pragma once

#include "global.h"

class Fluid
{
public:
    Fluid(float dt, float d, float v);

public:
    float density[gridSize * gridSize];
private:
    float pDensity[gridSize * gridSize]; //previous density
    float velocityX[gridSize * gridSize];
    float velocityY[gridSize * gridSize];
    float pVelocityX[gridSize * gridSize]; // previous velocity
    float pVelocityY[gridSize * gridSize];
    float timeStep;
    float diffusionAmt;
    float viscosity; //thickness

public:
    void Step();
    void Reset();
    void Fade();
    void AddDensity(int x, int y, float amount);
    void AddVelocity(int x, int y, float amountX, float amountY);
    void Diffuse(int b, float x[], float x0[], float factor, float timeStep);
    void LinearSolve(int b, float x[], float x0[], float a, float c);
    void Project(float velocX[], float velocY[], float velocX0[], float velocY0[]);
    void Advect(int b, float d[], float d0[], float velocX[], float velocY[], float timeStep);
    void SetBoundary(int b, float x[]);
};