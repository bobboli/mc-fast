#include "../MarchingCubes/MarchingCubes.h"
#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define CALIBRATE
#define FREQUENCY 2.8e9

double c_clock(MarchingCubes& mc, void (MarchingCubes::* ptr_update)(float), float r, MarchingCubesStage stage = ALL);

double gettickcount(MarchingCubes& mc, void (MarchingCubes::* ptr_update)(float), float r, MarchingCubesStage stage = ALL);

double queryperfcounter(MarchingCubes& mc, void (MarchingCubes::* ptr_update)(float), float r, LARGE_INTEGER f, MarchingCubesStage stage = ALL);