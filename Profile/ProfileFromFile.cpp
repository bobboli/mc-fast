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
    mc.setup(resX, resY, resZ, 2000000);

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
    if (argc != 6) { printf("usage: <filepath> <threshold> <blockX> <blockY> <blockZ>. \n"); return -1; }
    char* filepath = argv[1];
    float threshold = atof(argv[2]);
    int blockX = atoi(argv[3]);
    int blockY = atoi(argv[4]);
    int blockZ = atoi(argv[5]);

    // baseline result
    MarchingCubes mc;
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
    }
    else
    {
        SetupFromFile(mc, filepath);
    }
    mc.update(threshold);
    mc.exportObj("OutputMesh");
    
    // Blocking result
    MarchingCubes mc_b;
    mc_b.setBlocking(blockX, blockY, blockZ);
    if (isVolFile)
    {
        SetupFromVolFile(mc_b, filepath);
    }
    else
    {
        SetupFromFile(mc_b, filepath);
    }
    mc_b.update_block(threshold);
    mc_b.exportObj("OutputMesh_block");

    // vectorization result
    MarchingCubes mc_v;
    mc_v.setBlocking(1, 1, 8);
    if (isVolFile)
    {
        SetupFromVolFile(mc_v, filepath);
    }
    else
    {
        SetupFromFile(mc_v, filepath);
    }
    mc_v.update_vec(threshold);
    mc_v.exportObj("OutputMesh_vec");

    //MarchingCubes mc_v_16;
    //mc_v_16.setBlocking(1, 1, 16);
    //if (isVolFile)
    //{
    //    SetupFromVolFile(mc_v_16, filepath);
    //}
    //else
    //{
    //    SetupFromFile(mc_v_16, filepath);
    //}
    //mc_v_16.update_vec_16bit(threshold);
    //mc_v_16.exportObj("OutputMesh_vec16");

    // baseline timing
    void (MarchingCubes:: * ptr_update)(float) = &MarchingCubes::update;
    LARGE_INTEGER f;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f);
    double p = queryperfcounter(mc, ptr_update, threshold, f);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p / f.QuadPart, p / f.QuadPart * FREQUENCY);
    
    // blocking timing
    void (MarchingCubes:: * ptr_update_block)(float) = &MarchingCubes::update_block;
    LARGE_INTEGER f_b;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f_b);
    double p_b = queryperfcounter(mc_b, ptr_update_block, threshold, f_b);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_b / f_b.QuadPart, p_b / f_b.QuadPart * FREQUENCY);

    void (MarchingCubes:: * ptr_update_vec)(float) = &MarchingCubes::update_vec;
    LARGE_INTEGER f_v;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f_v);
    double p_v = queryperfcounter(mc_v, ptr_update_vec, threshold, f_v);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_v / f_v.QuadPart, p_v / f_v.QuadPart * FREQUENCY);

   /* void (MarchingCubes:: * ptr_update_vec_16bit)(float) = &MarchingCubes::update_vec_16bit;
    LARGE_INTEGER f_v_16;
    QueryPerformanceFrequency((LARGE_INTEGER*)&f_v_16);
    double p_v_16 = queryperfcounter(mc_v_16, ptr_update_vec_16bit, threshold, f_v_16);
    printf("Windows QueryPerformanceCounter() timing: %lf seconds. ==> %lf cycles based on FRENQUENCY.\n\n", p_v_16 / f_v_16.QuadPart, p_v_16 / f_v_16.QuadPart * FREQUENCY);*/

	return 0;
}
