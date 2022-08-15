#define OLC_PGE_APPLICATION

#include <math.h>
#include "olcPixelGameEngine.h"

constexpr int N = 100;
constexpr int iter = 4;

constexpr auto Get2DCoordinate(int x, int y) 
{
    if (x < 0) { x = 0; }
    if (x > N - 1) { x = N - 1; }

    if (y < 0) { y = 0; }
    if (y > N - 1) { y = N - 1; }

    return (y * N) + x;
}

class FluidGrid 
{
public:
    FluidGrid(float dt, float d, float v) 
    {
        timeStep = dt;
        diffusionAmt = d;
        viscosity = v;
    }

public:
    float density[N * N]    = {};
    float dDensity[N * N]   = {}; //previous density
    float velocityX[N * N]  = {};
    float velocityY[N * N]  = {};
    float velocityX0[N * N] = {};
    float velocityY0[N * N] = {};
    float timeStep; //time step
    float diffusionAmt; //diffusion amount
    float viscosity; //thickness of fluid

public:
    void step() 
    {
        Diffuse(1, velocityX0, velocityX, viscosity, timeStep);
        Diffuse(2, velocityY0, velocityY, viscosity, timeStep);

        Project(velocityX0, velocityY0, velocityX, velocityY);

        Advect(1, velocityX, velocityY0, velocityX0, velocityY0, timeStep);
        Advect(2, velocityY, velocityY0, velocityX0, velocityY0, timeStep);

        Project(velocityX, velocityY, velocityX0, velocityY0);

        Diffuse(0, dDensity, density, diffusionAmt, timeStep);
        Advect(0, density, dDensity, velocityX, velocityY, timeStep);

        for (int i = 0; i < N * N; i++)
        {
            if (density[i] > 255.0f)
                density[i] = 255.0f;
            if (dDensity[i] > 255.0f)
                dDensity[i] = 255.0f;
        }
    }

    void AddDensity(int x, int y, float amount) 
    {
        int index = Get2DCoordinate(x, y);
        density[index] += amount;
    }

    void AddVelocity(int x, int y, float amountX, float amountY) 
    {
        int index = Get2DCoordinate(x, y);
        velocityX[index] += amountX;
        velocityY[index] += amountY;
    }

    void Diffuse(int b, float x[], float x0[], float diffusionAmt, float timeStep) 
    {
        float a = timeStep * diffusionAmt * (N - 2) * (N - 2);
        this->LinearSolve(b, x, x0, a, 1 + 6 * a);
    }

    void LinearSolve(int b, float x[], float x0[], float a, float c) 
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

    void Project(float velocX[], float velocY[], float p[], float div[]) 
    {
        for (int j = 1; j < N - 1; j++)
            for (int i = 1; i < N - 1; i++) 
            {
                div[Get2DCoordinate(i, j)] = -0.5f * 
                     (velocX[Get2DCoordinate(i + 1, j)]
                    - velocX[Get2DCoordinate(i - 1, j)]
                    + velocY[Get2DCoordinate(i, j + 1)]
                    - velocY[Get2DCoordinate(i, j - 1)]) / N;

                p[Get2DCoordinate(i, j)] = 0;
            }

        this->SetBoundary(0, div);
        this->SetBoundary(0, p);
        this->LinearSolve(0, p, div, 1, 6);

        for (int j = 1; j < N - 1; j++)
            for (int i = 1; i < N - 1; i++) 
            {
                velocX[Get2DCoordinate(i, j)] -= 0.5f * (p[Get2DCoordinate(i + 1, j)] - p[Get2DCoordinate(i - 1, j)]) * N;
                velocY[Get2DCoordinate(i, j)] -= 0.5f * (p[Get2DCoordinate(i, j + 1)] - p[Get2DCoordinate(i, j - 1)]) * N;
            }

        this->SetBoundary(1, velocX);
        this->SetBoundary(2, velocY);
    }


    void Advect(int b, float d[], float d0[], float velocX[], float velocY[], float timeStep)
    {
        float i0, i1, j0, j1;

        float dtx = timeStep * (N - 2);
        float dty = timeStep * (N - 2);

        float s0, s1, t0, t1;
        float tmp1, tmp2, x, y;

        float Nfloat = N;
        float ifloat, jfloat;

        int i, j;

        for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
            for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[Get2DCoordinate(i, j)];
                tmp2 = dty * velocY[Get2DCoordinate(i, j)];
                x = ifloat - tmp1;
                y = jfloat - tmp2;

                if (x < 0.5f) x = 0.5f;
                if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
                i0 = ::floorf(x);
                i1 = i0 + 1.0f;
                if (y < 0.5f) y = 0.5f;
                if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
                j0 = ::floorf(y);
                j1 = j0 + 1.0f;

                s1 = x - i0;
                s0 = 1.0f - s1;
                t1 = y - j0;
                t0 = 1.0f - t1;

                int i0i = i0;
                int i1i = i1;
                int j0i = j0;
                int j1i = j1;

                d[Get2DCoordinate(i, j)] =
                    s0 * (t0 * d0[Get2DCoordinate(i0i, j0i)] + t1 * d0[Get2DCoordinate(i0i, j1i)]) +
                    s1 * (t0 * d0[Get2DCoordinate(i1i, j0i)] + t1 * d0[Get2DCoordinate(i1i, j1i)]);
            }
        }
        SetBoundary(b, d);
    }


    void SetBoundary(int b, float x[]) 
    {
        for (int i = 1; i < N - 1; i++) 
        {
            x[Get2DCoordinate(i, 0)]     = b == 2 ? -x[Get2DCoordinate(i, 1)] : x[Get2DCoordinate(i, 1)];
            x[Get2DCoordinate(i, N - 1)] = b == 2 ? -x[Get2DCoordinate(i, N - 2)] : x[Get2DCoordinate(i, N - 2)];
        }

        for (int j = 1; j < N - 1; j++) 
        {
            x[Get2DCoordinate(0, j)]     = b == 1 ? -x[Get2DCoordinate(1, j)] : x[Get2DCoordinate(1, j)];
            x[Get2DCoordinate(N - 1, j)] = b == 1 ? -x[Get2DCoordinate(N - 2, j)] : x[Get2DCoordinate(N - 2, j)];
        }

        x[Get2DCoordinate(0, 0)] = 0.33f * 
             (x[Get2DCoordinate(1, 0)]
            + x[Get2DCoordinate(0, 1)]
            + x[Get2DCoordinate(0, 0)]);

        x[Get2DCoordinate(0, N - 1)] = 0.33f * 
             (x[Get2DCoordinate(1, N - 1)]
            + x[Get2DCoordinate(0, N - 2)]
            + x[Get2DCoordinate(0, N - 1)]);

        x[Get2DCoordinate(N - 1, 0)] = 0.33f * 
             (x[Get2DCoordinate(N - 2, 0)]
            + x[Get2DCoordinate(N - 1, 1)]
            + x[Get2DCoordinate(N - 1, 0)]);

        x[Get2DCoordinate(N - 1, N - 1)] = 0.33f * 
             (x[Get2DCoordinate(N - 2, N - 1)]
            + x[Get2DCoordinate(N - 1, N - 2)]
            + x[Get2DCoordinate(N - 1, N - 1)]);
    }
};

class FluidSimulation2D : public olc::PixelGameEngine
{
public:
	FluidSimulation2D()
	{
		sAppName = "Fluid Simulation 2D";
	}

private:
    FluidGrid* Fluid;

    float RandFloat(float a, float b)
    {
        return ((b - a) * ((float)rand() / RAND_MAX)) + a;
    }

public:
	bool OnUserCreate() override
	{
        Fluid = new FluidGrid(0.2f, 0.0f, 0.0000001f);

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
        Clear(olc::BLACK);

        SetPixelMode(olc::Pixel::ALPHA);

        if (GetMouse(olc::Mouse::LEFT).bHeld)
        {
            Fluid->AddDensity(GetMouseX(), GetMouseY(), 200.0f);
            Fluid->AddVelocity(GetMouseX(), GetMouseY(), 10.0f, 0.0f);
        }

        Fluid->step();
        
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                Draw({ i , j }, olc::Pixel(255, 255, 255, uint8_t(Fluid->density[Get2DCoordinate(i, j)])));
		
        return true;
	}
};

int main()
{
	FluidSimulation2D _FluidSimulation2D;
	
	if (_FluidSimulation2D.Construct(N, N, 5, 5, false, true))
		_FluidSimulation2D.Start();

	return 0;
}