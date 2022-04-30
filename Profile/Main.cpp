// profile.cpp : Source file for your target.
//

#include "../MarchingCubes/MarchingCubes.h"

#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define NUM_RUNS 10
#define CYCLES_REQUIRED 1e9
#define FREQUENCY 2.8e9
#define CALIBRATE

double c_clock(MarchingCubes& mc, float r, MarchingCubesStage stage = ALL) {
    int i, num_runs;
    double cycles;
    clock_t start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14)) {
        start = clock();
        for (i = 0; i < num_runs; ++i) {
            mc.update(r);
        }
        end = clock();

        cycles = (double)(end - start);

        // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of CLOCKS_PER_SEC
        if (cycles >= CYCLES_REQUIRED / (FREQUENCY / CLOCKS_PER_SEC)) break;

        num_runs *= 2;
    }
#endif
    start = clock();
    for (i = 0; i < num_runs; ++i) {
        mc.update(r);
    }
    end = clock();

    return (double)(end - start) / num_runs;
}

double gettickcount(MarchingCubes& mc, float r, MarchingCubesStage stage = ALL) {
    int i, num_runs;
    double cycles, start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14)) {
        start = (double)GetTickCount();
        for (i = 0; i < num_runs; ++i) {
            mc.update(r);
        }
        end = (double)GetTickCount();

        cycles = (end - start) * FREQUENCY / 1e3; // end-start provides a measurement in the order of milliseconds

        if (cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = (double)GetTickCount();
    for (i = 0; i < num_runs; ++i) {
        mc.update(r);
    }
    end = (double)GetTickCount();

    return (end - start) / num_runs;
}

double queryperfcounter(MarchingCubes& mc, float r, LARGE_INTEGER f, MarchingCubesStage stage = ALL) {
    int i, num_runs;
    double cycles;
    LARGE_INTEGER start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14)) {
        QueryPerformanceCounter(&start);
        for (i = 0; i < num_runs; ++i) {
            mc.update(r);
        }
        QueryPerformanceCounter(&end);

        cycles = (double)(end.QuadPart - start.QuadPart);

        // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of f
        if (cycles >= CYCLES_REQUIRED / (FREQUENCY / f.QuadPart)) break;

        num_runs *= 2;
    }
#endif

    QueryPerformanceCounter(&start);
    for (i = 0; i < num_runs; ++i) {
        mc.update(r);
    }
    QueryPerformanceCounter(&end);

    return (double)(end.QuadPart - start.QuadPart) / num_runs;
}


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
    // 1. baseline result
    MarchingCubes mc;
    mc.setup(100, 100, 100);
    SetSphere(mc);
    const float radius = 0.4;
    mc.update(radius);
    mc.exportObj("Sphere");

    
    // 2. timing
    MarchingCubes mc_s;
    mc_s.setup(100, 100, 100);
    SetSphere(mc_s);

    double c = c_clock(mc_s, radius);
    printf("C clock() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY. \n\n", c / CLOCKS_PER_SEC, c / CLOCKS_PER_SEC * FREQUENCY);

    LARGE_INTEGER f;
    double t = gettickcount(mc_s, radius);
    printf("Windows getTickCount() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", t / 1000.0, t / 1000.0 * FREQUENCY);

    QueryPerformanceFrequency((LARGE_INTEGER*)&f);

    double p = queryperfcounter(mc_s, radius, f);
    printf("Windows QueryPerformanceCounter()  %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p / f.QuadPart, p / f.QuadPart * FREQUENCY);
    
    mc_s.exportObj("Sphere_smoothed");

	return 0;
}
