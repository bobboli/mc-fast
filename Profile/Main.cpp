// profile.cpp : Source file for your target.
//

#include "../MarchingCubes/MarchingCubes.h"

void SetSphere(MarchingCubes& mc)
{
	const int cX = mc.resX / 2;
	const int cY = mc.resY / 2;
	const int cZ = mc.resZ / 2;
	const float scale = 1.0f/min(min(mc.resX, mc.resY), mc.resZ);
	for (int i = 0; i < mc.resX; i++)
	{
		for (int j = 0; j < mc.resY; j++)
		{
			for (int k = 0; k < mc.resZ; k++)
			{
				float val = (i - cX) * (i - cX) + (j - cY) * (j - cY) + (k - cZ) * (k - cZ);
				val *= scale;
				mc.setIsoValue(i, j, k, val);
			}
		}
	}
}

int main()
{
	MarchingCubes mc;
	mc.setup(100, 100, 100);
	SetSphere(mc);
	const float radius = 0.4;
	mc.update(radius);
	mc.exportObj("Sphere");
	return 0;
}
