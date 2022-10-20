#define OLC_PGE_APPLICATION

#include "olcPixelGameEngine.h"
#include "global.h"
#include "fluid.h"

class FluidSimulation2D : public olc::PixelGameEngine
{
public:
	FluidSimulation2D()
	{
		sAppName = "Fluid Simulation 2D";
	}

private:
    Fluid* _Fluid;
    olc::vi2d mousePos;
    olc::vi2d pMousePos;

    float RandFloat(float a, float b)
    {
        return ((b - a) * ((float)rand() / RAND_MAX)) + a;
    }

    void Render()
    {
        Clear(olc::BLACK);

        for (int i = 0; i < gridSize; i++)
            for (int j = 0; j < gridSize; j++)
                Draw({ i, j }, olc::Pixel(255, 255, 255, uint8_t(_Fluid->density[Get2DCoordinate(i, j)])));

        DrawRect({ 0, 0 }, { ScreenWidth() - 1, ScreenHeight() - 21 }, olc::VERY_DARK_BLUE);
    }

    void Input()
    {
        mousePos = GetMousePos();

        if (GetMouse(olc::Mouse::LEFT).bHeld)
        {
            _Fluid->AddDensity(GetMouseX(), GetMouseY(), 200.0f);
            _Fluid->AddVelocity(GetMouseX(), GetMouseY(), (mousePos.x - pMousePos.x) / 4.0f, (mousePos.y - pMousePos.y) / 4.0f);
        }

        pMousePos = mousePos;

        if (GetKey(olc::Key::R).bPressed)
        {
            _Fluid->Reset();
        }
    }

public:
	bool OnUserCreate() override
	{
        _Fluid = new Fluid(0.5f, 0.000001f, 0.0000001f);
        mousePos = GetMousePos();
        pMousePos = GetMousePos();
        SetPixelMode(olc::Pixel::ALPHA);

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
        _Fluid->Step();
        Input();
        Render();

        return true;
	}
};

int main()
{
	FluidSimulation2D _FluidSimulation2D;
	
	if (_FluidSimulation2D.Construct(gridSize, gridSize + 20, 2, 2, true))
		_FluidSimulation2D.Start();

	return 0;
}