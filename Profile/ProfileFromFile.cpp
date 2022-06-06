#include "Profile.h"
#include <intrin.h>

// regular text file
// first 3 int as resolution x y z, then float data
void SetupFromFile(MarchingCubes & mc, char* filename)
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

// vol file from https://web.cs.ucdavis.edu/~okreylos/PhDStudies/Spring2000/ECS277/
// big_endian binary
void SetupFromVolFile(MarchingCubes& mc, char* filename)
{
    ifstream ifs;
    ifs.open(filename, ios::in | ios::binary);
    uint32_t resX, resY, resZ;
    uint32_t big_endian;
    ifs.read((char *)&big_endian, sizeof(uint32_t));
    // gcc use uint32_t __builtin_bswap32 (uint32_t x)
    resX = _byteswap_ulong(big_endian);
    ifs.read((char*)&big_endian, sizeof(uint32_t));
    resY = _byteswap_ulong(big_endian);
    ifs.read((char*)&big_endian, sizeof(uint32_t));
    resZ = _byteswap_ulong(big_endian);
    mc.setup(resX, resY, resZ);

    // useless data
    ifs.read((char*)&big_endian, sizeof(uint32_t));
    uint32_t dummy = _byteswap_ulong(big_endian);

    float a;
    ifs.read((char*)&a, sizeof(float));
    ifs.read((char*)&a, sizeof(float));
    ifs.read((char*)&a, sizeof(float));

    
    cout << resX << " " << resY << " " << resZ << endl;

    for (int i = 0; i < resX; i++)
    {
        for (int j = 0; j < resY; j++)
        {
            for (int k = 0; k < resZ; k++)
            {
                int c;
                ifs.read((char*)&c, 1);
                mc.setIsoValue(i, j, k, c);
            }
        }
    }
    ifs.close();
}

int main(int argc, char** argv)
{
    if (argc != 3) { printf("usage: <filepath> <threshold>. \n"); return -1; }
    char* filepath = argv[1];
    float threshold = atof(argv[2]);

    // baseline result
    MarchingCubes mc;
    MarchingCubes mc_l;
    MarchingCubes mc_lv;
    int len = strlen(filepath);
    char ext[4];
    ext[0] = filepath[len - 3];
    ext[1] = filepath[len - 2];
    ext[2] = filepath[len - 1];
    ext[3] = '\0';
   
    bool isVolFile = (strcmp(ext, "vol") == 0);
    if(isVolFile)
    {
        SetupFromVolFile(mc, filepath);
        SetupFromVolFile(mc_l, filepath);
        SetupFromVolFile(mc_lv, filepath);
    }
    else
    {
        SetupFromFile(mc, filepath);
        SetupFromFile(mc_l, filepath);
        SetupFromFile(mc_lv, filepath);
    }
    mc.update(threshold);
    mc.exportObj("OutputMesh");

    mc_l.update_level(threshold);
    mc_l.exportObj("OutputMesh_level"); 

    mc_lv.update_level_vec(threshold);
    mc_lv.exportObj("OutputMesh_level_vec");

    // baseline timing
    void (MarchingCubes:: * ptr_update)(float) = &MarchingCubes::update;
    LARGE_INTEGER f;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f);
    double p = queryperfcounter(mc, ptr_update, threshold, f);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p / f.QuadPart, p / f.QuadPart * FREQUENCY);
    
    // Level-by-level timing
    void (MarchingCubes:: * ptr_update_level)(float) = &MarchingCubes::update_level;
    LARGE_INTEGER f_l;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f_l);
    double p_l = queryperfcounter(mc_l, ptr_update_level, threshold, f_l);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_l / f_l.QuadPart, p_l / f_l.QuadPart * FREQUENCY);

    // Level-by-level vectorization timing
    void (MarchingCubes:: * ptr_update_level_vec)(float) = &MarchingCubes::update_level_vec;
    LARGE_INTEGER f_lv;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f_lv);
    double p_lv = queryperfcounter(mc_lv, ptr_update_level_vec, threshold, f_lv);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_lv / f_lv.QuadPart, p_lv / f_lv.QuadPart * FREQUENCY);

	return 0;
}
