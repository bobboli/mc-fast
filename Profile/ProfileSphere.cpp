// profile.cpp : Source file for your target.
//

#include "Profile.h"

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
                val = sqrt(val);
				val *= scale;
				mc.setIsoValue(i, j, k, val);
			}
		}
	}
}

int main(int argc, char** argv)
{
    if (argc != 6) { printf("usage: <res> <threshold> <blockX> <blockY> <blockZ>. \n"); return -1; }
    int res = atoi(argv[1]);
    float radius = atof(argv[2]);
    int blockX = atoi(argv[3]);
    int blockY = atoi(argv[4]);
    int blockZ = atoi(argv[5]);
    printf("resolution=%d\n", res);

    // 1. baseline result
    MarchingCubes mc;
    mc.setup(res, res, res);
    SetSphere(mc);
    //const float radius = 0.4;
    mc.update(radius);
    mc.exportObj("Sphere");
    
    // 2. Blocking result
    MarchingCubes mc_b;
    mc_b.setup(res, res, res);
    mc_b.setBlocking(blockX, blockY, blockZ);
    SetSphere(mc_b);
    //const float radius = 0.4;
    mc_b.update_block(radius);
    mc_b.exportObj("Sphere_block");


    // 3. baseline timing
    void (MarchingCubes:: * ptr_update)(float) = &MarchingCubes::update;
    double c = c_clock(mc, ptr_update, radius);
    printf("C clock() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY. \n\n", c / CLOCKS_PER_SEC, c / CLOCKS_PER_SEC * FREQUENCY);

    LARGE_INTEGER f;
    double t = gettickcount(mc, ptr_update, radius);
    printf("Windows getTickCount() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", t / 1000.0, t / 1000.0 * FREQUENCY);

    QueryPerformanceFrequency((LARGE_INTEGER*)&f);
    double p = queryperfcounter(mc, ptr_update, radius, f);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p / f.QuadPart, p / f.QuadPart * FREQUENCY);
    

    // blocking timing
    void (MarchingCubes:: * ptr_update_block)(float) = &MarchingCubes::update_block;
    double c_b = c_clock(mc_b, ptr_update_block, radius);
    printf("C clock() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY. \n\n", c_b / CLOCKS_PER_SEC, c_b / CLOCKS_PER_SEC * FREQUENCY);

    LARGE_INTEGER f_b;
    double t_b = gettickcount(mc_b, ptr_update_block, radius);
    printf("Windows getTickCount() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", t_b / 1000.0, t_b / 1000.0 * FREQUENCY);

    QueryPerformanceFrequency((LARGE_INTEGER*)&f_b);
    double p_b = queryperfcounter(mc_b, ptr_update_block, radius, f_b);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_b / f_b.QuadPart, p_b / f_b.QuadPart * FREQUENCY);


    // 4. count floating point operations
    operation_counts counts;
    mc.count_ops(radius, counts);
    counts.print();

	return 0;
}
