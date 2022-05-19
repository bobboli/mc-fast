#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define NUM_RUNS 10
#define CYCLES_REQUIRED 1e8
#define CALIBRATE
#define FREQUENCY 2.8e9

double c_clock(MarchingCubes& mc, void (MarchingCubes::* ptr_update)(float), float r, MarchingCubesStage stage = ALL) {
    int i, num_runs;
    double cycles;
    clock_t start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14)) {
        start = clock();
        for (i = 0; i < num_runs; ++i) {
            (mc.*ptr_update)(r);
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
        (mc.*ptr_update)(r);
    }
    end = clock();

    return (double)(end - start) / num_runs;
}

double gettickcount(MarchingCubes& mc, void (MarchingCubes::* ptr_update)(float), float r, MarchingCubesStage stage = ALL) {
    int i, num_runs;
    double cycles, start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14)) {
        start = (double)GetTickCount();
        for (i = 0; i < num_runs; ++i) {
            (mc.*ptr_update)(r);
        }
        end = (double)GetTickCount();

        cycles = (end - start) * FREQUENCY / 1e3; // end-start provides a measurement in the order of milliseconds

        if (cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = (double)GetTickCount();
    for (i = 0; i < num_runs; ++i) {
        (mc.*ptr_update)(r);
    }
    end = (double)GetTickCount();

    return (end - start) / num_runs;
}

double queryperfcounter(MarchingCubes& mc, void (MarchingCubes::* ptr_update)(float), float r, LARGE_INTEGER f, MarchingCubesStage stage = ALL) {
    int i, num_runs;
    double cycles;
    LARGE_INTEGER start, end;

    num_runs = NUM_RUNS;
#ifdef CALIBRATE
    while (num_runs < (1 << 14)) {
        QueryPerformanceCounter(&start);
        for (i = 0; i < num_runs; ++i) {
            (mc.*ptr_update)(r);
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
        (mc.*ptr_update)(r);
    }
    QueryPerformanceCounter(&end);

    return (double)(end.QuadPart - start.QuadPart) / num_runs;
}

