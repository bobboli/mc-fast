﻿// profile.cpp : Source file for your target.
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

void SetRandom(MarchingCubes &mc, float _lo=-1.0, float _hi=1.0)
{
    for (int i = 0; i < mc.resX; i++)
	{
		for (int j = 0; j < mc.resY; j++)
		{
			for (int k = 0; k < mc.resZ; k++)
			{
				float val = _lo + static_cast<float>(rand()) / (static_cast <float> (RAND_MAX/(_hi-_lo)));
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

    // baseline result
    MarchingCubes mc;
    mc.setup(res, res, res, 15 * res * res * res);
    // SetSphere(mc);
    SetRandom(mc);
    //const float radius = 0.4;
    mc.update(radius);
    mc.exportObj("Sphere");
    
    // blocking result
    MarchingCubes mc_b;
    mc_b.setup(res, res, res, 15 * res * res * res);
    mc_b.setBlocking(blockX, blockY, blockZ);
    // SetSphere(mc_b);
    SetRandom(mc_b);
    mc_b.update_block(radius);
    mc_b.exportObj("Sphere_block");
    
    // blocking (new) result
	MarchingCubes mc_b_new;
	mc_b_new.setup(res, res, res, 15 * res * res * res);
	mc_b_new.setBlocking(blockX, blockY, blockZ);
	// SetSphere(mc_b_new);
    SetRandom(mc_b_new);
	mc_b_new.update_block_new(radius);
	mc_b_new.exportObj("Sphere_block_new");


    // vectorization result
    MarchingCubes mc_v;
    mc_v.setup(res, res, res, 15 * res * res * res);
    mc_v.setBlocking(1, 1, 8);
    // SetSphere(mc_v);
    SetRandom(mc_v);
    mc_v.update_vec(radius);
    mc_v.exportObj("Sphere_vec");

    // baseline timing
    void (MarchingCubes:: * ptr_update)(float) = &MarchingCubes::update;
    LARGE_INTEGER f;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f);
    double p = queryperfcounter(mc, ptr_update, radius, f);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p / f.QuadPart, p / f.QuadPart * FREQUENCY);

    // blocking timing
    void (MarchingCubes:: * ptr_update_block)(float) = &MarchingCubes::update_block;
    LARGE_INTEGER f_b;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f_b);
    double p_b = queryperfcounter(mc_b, ptr_update_block, radius, f_b);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_b / f_b.QuadPart, p_b / f_b.QuadPart * FREQUENCY);

	// blocking (new) timing
	void (MarchingCubes:: * ptr_update_block_new)(float) = &MarchingCubes::update_block_new;
	LARGE_INTEGER f_b_new;
	QueryPerformanceFrequency((LARGE_INTEGER*)&f_b_new);
	double p_b_new = queryperfcounter(mc_b_new, ptr_update_block_new, radius, f_b_new);
	printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_b_new / f_b_new.QuadPart, p_b_new / f_b_new.QuadPart * FREQUENCY);

    // vectorization timing
    void (MarchingCubes:: * ptr_update_vec)(float) = &MarchingCubes::update_vec;
    LARGE_INTEGER f_v;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f_v);
    double p_v = queryperfcounter(mc_v, ptr_update_vec, radius, f_v);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_v / f_v.QuadPart, p_v / f_v.QuadPart * FREQUENCY);

    // count baseline floating point operations
    operation_counts counts;
    mc.count_ops(radius, counts);
    counts.print();

	return 0;
}
