#define OLC_PGE_APPLICATION

#include <math.h>
#include "olcPixelGameEngine.h"

constexpr int Nx = 168;
constexpr int Ny = 105;
constexpr int iter = 4;

template<typename T1, typename T2>
constexpr auto IX(T1 x, T2  y) 
{
    if (y > Nx - 1)
        x = Nx - 1;
    else if (x < 0)
        x = 0;
    if (y > Ny - 1)
        y = Ny - 1;
    else if (y < 0)
        y = 0;
    return x + (y * Nx); 
}

class FluidGrid 
{
public:
    FluidGrid(float dt, float diffusion, float viscosity) 
    {
        dt = dt;
        diff = diffusion;
        visc = viscosity;

        /*s = new float[Nx * Ny]{};
        density = new float[Nx * Ny]{};

        Vx = new float[Nx * Ny]{};
        Vy = new float[Nx * Ny]{};

        Vx0 = new float[Nx * Ny]{};
        Vy0 = new float[Nx * Ny]{};*/
    }

public:
    float dt; //time step
    float diff; //diffusion amount
    float visc; //thickness of fluid
         
    float s[Nx * Ny] = {}; //previous density
    float density[Nx * Ny] = {};
         
    float Vx[Nx * Ny] = {};
    float Vy[Nx * Ny] = {};
         
    float Vx0[Nx * Ny] = {};
    float Vy0[Nx * Ny] = {};

public:
    void step() 
    {
        diffuse(1, Vx0, Vx, visc, dt);
        diffuse(2, Vy0, Vy, visc, dt);

        project(Vx0, Vy0, Vx, Vy);

        advect(1, Vx, Vy0, Vx0, Vy0, dt);
        advect(2, Vy, Vy0, Vx0, Vy0, dt);

        project(Vx, Vy, Vx0, Vy0);

        diffuse(0, s, density, diff, dt);
        advect(0, density, s, Vx, Vy, dt);

        for (int i = 0; i < Nx * Ny; i++)
        {
            if (density[i] > 255.0f)
                density[i] = 255.0f;
            if (s[i] > 255.0f)
                s[i] = 255.0f;
        }
    }

    void addDensity(int x, int y, float amount) 
    {
        int index = IX(x, y);
        density[index] += amount;
    }

    void addVelocity(int x, int y, float amountX, float amountY) 
    {
        int index = IX(x, y);
        Vx[index] += amountX;
        Vy[index] += amountY;
    }

    void diffuse(int b, float x[], float x0[], float diff, float dt) 
    {
        float a = dt * diff * (Nx - 2) * (Ny - 2);
        lin_solve(b, x, x0, a, 1 + 6 * a);
    }

    void lin_solve(int b, float x[], float x0[], float a, float c) {
        float cRecip = 1.0f / c;
        for (int k = 0; k < iter; k++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int i = 1; i < Nx - 1; i++) {
                    x[IX(i, j)] =
                        (x0[IX(i, j)]
                            + a *
                            (x[IX(i + 1, j)]
                                + x[IX(i - 1, j)]
                                + x[IX(i, j + 1)]
                                + x[IX(i, j - 1)]
                                )) * cRecip;
                }
            }
        }
        set_bnd(b, x);
    }

    void project(float velocX[], float velocY[], float p[], float div[]) 
    {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                div[IX(i, j)] = -0.5f * (
                      velocX[IX(i + 1, j)]
                    - velocX[IX(i - 1, j)]
                    + velocY[IX(i, j + 1)]
                    - velocY[IX(i, j - 1)]
                    ) / ((Nx + Ny) * 0.5);
                p[IX(i, j)] = 0;
            }
        }
        set_bnd(0, div);
        set_bnd(0, p);
        lin_solve(0, p, div, 1.0f, 6.0f);

        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)]
                    - p[IX(i - 1, j)]) * Nx;
                velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)]
                    - p[IX(i, j - 1)]) * Ny;
            }
        }

        set_bnd(1, velocX);
        set_bnd(2, velocY);
    }


    void advect(int b, float d[], float d0[], float velocX[], float velocY[], float dt)
    {
        float i0, i1, j0, j1;
        float s0, s1, t0, t1;
        float x, y;

        float dtx = dt * (Nx - 2);
        float dty = dt * (Ny - 2);


        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                x = float(i) - dtx * velocX[IX(i, j)];
                y = float(j) - dty * velocY[IX(i, j)];

                if (x < 0.5f) x = 0.5f;
                if (x > float(Nx) + 0.5f) x = float(Nx) + 0.5f;
                i0 = floor(x);
                i1 = i0 + 1.0f;
                if (y < 0.5f) y = 0.5f;
                if (y > float(Ny) + 0.5f) y = float(Ny) + 0.5f;
                j0 = floor(y);
                j1 = j0 + 1.0f;

                s1 = x - i0;
                s0 = 1.0f - s1;
                t1 = y - j0;
                t0 = 1.0f - t1;

                d[IX(i, j)] =
                    s0 * (t0 * d0[int(IX(i0, j0))] + t1 * d0[int(IX(i0, j1))]) +
                    s1 * (t0 * d0[int(IX(i1, j0))] + t1 * d0[int(IX(i1, j1))]);
            }
        }
        set_bnd(b, d);
    }


    void set_bnd(int b, float x[]) 
    {
        for (int i = 1; i < Ny - 1; i++) {
            x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
            x[IX(Nx - 1, i)] = b == 1 ? -x[IX(Nx - 1, i)] : x[IX(Nx - 1, i)];
        }
        for (int i = 1; i < Nx - 1; i++) {
            x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, Ny - 1)] = b == 2 ? -x[IX(i, Ny - 1)] : x[IX(i, Ny - 1)];
        }
        x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, Ny - 1)] = 0.5f * (x[IX(1, Ny - 1)] + x[IX(0, Ny - 1)]);
        x[IX(Nx - 1, 0)] = 0.5f * (x[IX(Nx - 2, 0)] + x[IX(Nx - 1, 1)]);
        x[IX(Nx - 1, Ny - 1)] = 0.5f * (x[IX(Nx - 2, Ny - 1)] + x[IX(Nx - 1, Ny - 2)]);
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
    FluidGrid* fluid;

    float RandFloat(float a, float b)
    {
        return ((b - a) * ((float)rand() / RAND_MAX)) + a;
    }

public:
	bool OnUserCreate() override
	{
        fluid = new FluidGrid(0.5f, 0.0f, 0.00001f);

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
        Clear(olc::BLACK);

        SetPixelMode(olc::Pixel::ALPHA);

        DrawStringDecal({ 10, 10 }, std::to_string(fluid->Vx[IX(GetMouseX(), GetMouseY())]) + ", " + std::to_string(fluid->Vy[IX(GetMouseX(), GetMouseY())]), olc::WHITE, { 0.5f,0.5f });
        DrawStringDecal({ 10, 20 }, std::to_string(fluid->density[IX(GetMouseX(), GetMouseY())]), olc::WHITE, {0.5f,0.5f});

        fluid->addDensity(GetMouseX(), GetMouseY(), 200.0f);
        fluid->addVelocity(GetMouseX(), GetMouseY(), 100.0f, 100.0f);

        fluid->step();
        
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                Draw({ i , j }, olc::Pixel(255, 255, 255, uint8_t(fluid->density[IX(i, j)])/*255U * (fluid->density[IX(i, j)] / 255)*/));
		
        return true;
	}
};

int main()
{
	FluidSimulation2D _FluidSimulation2D;
	
	if (_FluidSimulation2D.Construct(Nx, Ny, 5, 5, false, true))
		_FluidSimulation2D.Start();

	return 0;
}