#include "Profile.h"

void SetupFromFile(MarchingCubes & mc, string filename)
{
    ifstream ifs;
    ifs.open(filename, ios::in);

    if (!ifs.is_open())
    {
        cout << "can't open file." << endl;
        return;
    }

    int resX, resY, resZ;
    ifs >> resX >> resY >> resZ;
    mc.setup(resX, resY, resZ);
    cout << resX << " " << resY << " " << resZ << endl;

    for (int i = 0; i < resX; i++)
    {
        for (int j = 0; j < resY; j++)
        {
            for (int k = 0; k < resZ; k++)
            {
                float val = 0.0;
                ifs >> val;
                mc.setIsoValue(i, j, k, val);
            }
        }
    }
}

int main(int argc, char** argv)
{
    if (argc != 3) { printf("usage: <filename> <threshold>. \n"); return -1; }
    string filename = argv[1];
    float threshold = atof(argv[2]);

    // 1. baseline result
    MarchingCubes mc;
    //const float radius = 0.4;
    SetupFromFile(mc, filename);
    mc.update(threshold);
    mc.exportObj("Fluid");
    
    // 2. Blocking result
    //MarchingCubes mc_b;
    //mc_b.setup(res, res, res);
    //mc_b.setBlocking(blockX, blockY, blockZ);
    ////const float radius = 0.4;
    //mc_b.update_block(radius);
    //mc_b.exportObj("Sphere_block");


    // 3. baseline timing
    void (MarchingCubes:: * ptr_update)(float) = &MarchingCubes::update;
    LARGE_INTEGER f;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f);
    double p = queryperfcounter(mc, ptr_update, threshold, f);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p / f.QuadPart, p / f.QuadPart * FREQUENCY);
    
    // blocking timing
   /* void (MarchingCubes:: * ptr_update_block)(float) = &MarchingCubes::update_block;
    LARGE_INTEGER f_b;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f_b);
    double p_b = queryperfcounter(mc_b, ptr_update_block, radius, f_b);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_b / f_b.QuadPart, p_b / f_b.QuadPart * FREQUENCY);*/

	return 0;
}
