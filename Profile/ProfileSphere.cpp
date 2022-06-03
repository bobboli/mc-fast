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

void SetRandom(MarchingCubes& mc, float _lo = -1.0, float _hi = 1.0)
{
	srand(1);
	for (int i = 0; i < mc.resX; i++)
	{
		for (int j = 0; j < mc.resY; j++)
		{
			for (int k = 0; k < mc.resZ; k++)
			{
				float val = _lo + static_cast<float>(rand()) / (static_cast <float> (RAND_MAX / (_hi - _lo)));
				mc.setIsoValue(i, j, k, val);
			}
		}
	}
}

int main(int argc, char** argv)
{
	//if (argc != 6) { printf("usage: <res> <threshold> <blockX> <blockY> <blockZ>. \n"); return -1; }
	//int res = atoi(argv[1]);
	//float radius = atof(argv[2]);
	//int blockX = atoi(argv[3]);
	//int blockY = atoi(argv[4]);
	//int blockZ = atoi(argv[5]);
	//printf("resolution=%d\n", res);

	float radius = 0.4;
	int resx = 100, resy = 100, resz= 100;
	int blockX = 10, blockY = 10, blockZ = 10;

    // baseline result
    MarchingCubes mc;
    mc.setup(resx, resy, resz);
    //SetSphere(mc);
	SetRandom(mc);
    //const float radius = 0.4;
    mc.update(radius);
	//mc.exportObj("Sphere");
    printf("%d\n", mc.vertexCount);
    
    // blocking result
    //MarchingCubes mc_b;
    //mc_b.setup(res, res, res);
    //mc_b.setBlocking(blockX, blockY, blockZ);
    //SetSphere(mc_b);
    //mc_b.update_block(radius);
    //mc_b.exportObj("Sphere_block");
    
    // blocking (new) result
	//MarchingCubes mc_b_new;
	//mc_b_new.setup(resx, resy, resz);
	//mc_b_new.setBlocking(blockX, blockY, blockZ);
	//SetRandom(mc_b_new);
	//mc_b_new.update_block_new(radius);
	//mc_b_new.exportObj("Sphere_block_new");

    // Level-by-level result
	MarchingCubes mc_l;
    mc_l.setup(resx, resy, resz);
	//SetSphere(mc_l);
	SetRandom(mc_l);
    mc_l.update_level(radius);
    printf("%d\n", mc_l.vertexCount);
	//mc_l.exportObj("Sphere_level");

	// Level-by-level no switch result
	MarchingCubes mc_ln;
	mc_ln.setup(resx, resy, resz);
	//SetSphere(mc_l);
	SetRandom(mc_ln);
	mc_ln.update_level_noswitch(radius);
	printf("%d\n", mc_ln.vertexCount);
	//mc_l.exportObj("Sphere_level_noswitch");

	// Level-by-level vectorization result
	MarchingCubes mc_lv;
	mc_lv.setup(resx, resy, resz);
	//SetSphere(mc_lv);
	SetRandom(mc_lv);
	mc_lv.update_level_vec(radius);
	printf("%d\n", mc_lv.vertexCount);
	//mc_lv.exportObj("Sphere_level");

    // vectorization result
    MarchingCubes mc_v;
    mc_v.setup(resx, resy, resz);
    mc_v.setBlocking(1, 1, 8);
    //SetSphere(mc_v);
	SetRandom(mc_v);
    mc_v.update_vec(radius);
	printf("%d\n", mc_v.vertexCount);
    //mc_v.exportObj("Sphere_vec");

    // baseline timing
    void (MarchingCubes:: * ptr_update)(float) = &MarchingCubes::update;
    LARGE_INTEGER f;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f);
    double p = queryperfcounter(mc, ptr_update, radius, f);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p / f.QuadPart, p / f.QuadPart * FREQUENCY);

    // blocking timing
 /*   void (MarchingCubes:: * ptr_update_block)(float) = &MarchingCubes::update_block;
	LARGE_INTEGER f_b;
	QueryPerformanceFrequency((LARGE_INTEGER*)&f_b);
	double p_b = queryperfcounter(mc_b, ptr_update_block, radius, f_b);
	printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_b / f_b.QuadPart, p_b / f_b.QuadPart * FREQUENCY);*/

	// blocking (new) timing
	//void (MarchingCubes:: * ptr_update_block_new)(float) = &MarchingCubes::update_block_new;
	//LARGE_INTEGER f_b_new;
	//QueryPerformanceFrequency((LARGE_INTEGER*)&f_b_new);
	//double p_b_new = queryperfcounter(mc_b_new, ptr_update_block_new, radius, f_b_new);
	//printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_b_new / f_b_new.QuadPart, p_b_new / f_b_new.QuadPart * FREQUENCY);

    // Level-by-level timing
	void (MarchingCubes:: * ptr_update_level)(float) = &MarchingCubes::update_level;
	LARGE_INTEGER f_l;
	QueryPerformanceFrequency((LARGE_INTEGER*)&f_l);
	double p_l = queryperfcounter(mc_l, ptr_update_level, radius, f_l);
	printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_l / f_l.QuadPart, p_l / f_l.QuadPart * FREQUENCY);

	// Level-by-level no switch timing
	void (MarchingCubes::* ptr_update_level_noswitch)(float) = &MarchingCubes::update_level_noswitch;
	LARGE_INTEGER f_ln;
	QueryPerformanceFrequency((LARGE_INTEGER*)&f_ln);
	double p_ln = queryperfcounter(mc_ln, ptr_update_level_noswitch, radius, f_ln);
	printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_ln / f_ln.QuadPart, p_ln / f_ln.QuadPart * FREQUENCY);

	// Level-by-level vectorization timing
	void (MarchingCubes:: * ptr_update_level_vec)(float) = &MarchingCubes::update_level_vec;
	LARGE_INTEGER f_lv;
	QueryPerformanceFrequency((LARGE_INTEGER*)&f_lv);
	double p_lv = queryperfcounter(mc_lv, ptr_update_level_vec, radius, f_lv);
	printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_lv / f_lv.QuadPart, p_lv / f_lv.QuadPart * FREQUENCY);


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
