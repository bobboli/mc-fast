
// http://paulbourke.net/geometry/polygonise/
//

// Implemented based on ofxMarchingCubes:
// https://github.com/larsberg/ofxMarchingCubes
// Modifications:
// Replace ofxVec3f with Vector3f to get rid of the dependency of the ofx.
// Ignore transformation matrices.
// Remove GUI related stuff.

/*
TODO::
 -get worldposition in grid
 -add iso value at world position
 */

#pragma once

#include "mcTables.h"
#include "Geometry.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include <immintrin.h>

using namespace std;

enum MarchingCubesStage
{
	ALL,
	STAGE1,
	STAGE2,
	STAGE3,
};

class operation_counts
{
public:
	int fl_cmp;
	int fl_add;
	int fl_mul;
	int fl_div;
	int fl_abs;
	int empty_cells;

	operation_counts()
	{
		fl_cmp = 0;
		fl_add = 0;
		fl_mul = 0;
		fl_div = 0;
		fl_abs = 0;
		empty_cells = 0;
	}

	void print()
	{
		printf("floating point cmp: %d, \nadd: %d,\nmul: %d,\ndiv: %d,\nfabs: %d,\ntotal number of flops: %d\n",
			   fl_cmp,
			   fl_add,
			   fl_mul,
			   fl_div,
			   fl_abs,
			   fl_cmp + fl_add + fl_mul + fl_div + fl_abs);
		printf("total number of empty cells: %d", empty_cells);
	}
};

class MarchingCubes
{
public:
	MarchingCubes();
	~MarchingCubes();

	void setMaxVertexCount(int _maxVertexCount = 100000);

	void setup(int resX = 30, int resY = 20, int resZ = 30, int _maxVertexCount = 200000);
	void setBlocking(int blockX, int blockY, int blockZ);
	void update() { update(threshold); }
	void update_block() { update_block(threshold); }
	void update_vec() { update_vec(threshold); }
	void update(float _threshold);
	void update_block(float _threshold);
	void update_vec(float _threshold);
	void update_vec_16bit(float _threshold);
	void update_block_new(float _threshold);
	void update_block_new() { update_block_new(threshold); }

	void MarchingCubes::update_level(float _threshold)
	{
		threshold = _threshold;

		vertexCount = 0;

		for (int x = 0; x < sx; ++x)
		{
			polygonise_level(x);
		}

	}
	void update_level() { update_level(threshold); }

	void count_ops(float _threshold, operation_counts &counts);

	// void draw( GLenum renderType = GL_TRIANGLES );
	// void drawWireframe();
	// void drawArrays( vector<Vector3f>* _vertices, vector<Vector3f>* _normals=NULL );
	// void drawGrid( bool drawGridPoints=true);

	void flipNormals() { flipNormalsValue *= -1; }
	void setResolution(int _x = 10, int _y = 10, int _z = 10);
	void polygonise(int i, int j, int k);
	void polygonise_block(int i, int j, int k, int bX, int bY, int bZ);
	void polygonise_block_new(int ibx, int iby, int ibz);
	void polygonise_level(int level);


	void polygonise_vec(int i, int j, int k, int bX, int bY, int bZ);
	void polygonise_vec_16bit(int i, int j, int k, int bX, int bY, int bZ);
	void polygonise_count_ops(int i, int j, int k, operation_counts &counts);
	void computeNormal(int i, int j, int k);
	inline void vertexInterp(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, Vector3f &v, Vector3f &n);
	inline void vertexInterp_X(float threshold, int x1, int x2, int y, int z, Vector3f &v, Vector3f &n);
	inline void vertexInterp_Y(float threshold, int x, int y1, int y2, int z, Vector3f &v, Vector3f &n);
	inline void vertexInterp_Z(float threshold, int x, int y, int z1, int z2, Vector3f &v, Vector3f &n);

	void vertexInterp_vec(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, Vector3f &v, Vector3f &n);
	void vertexInterp_count_ops(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, operation_counts &counts);

	void setIsoValue(int x, int y, int z, float value);

	void encodeIsoValsMorton();

	inline int getEdgeIdxX(int x, int y, int z)
	{
		return x * by1 * bz1 + y * bz1 + z;
	}

	inline int getEdgeIdxY(int x, int y, int z)
	{
		return x * by * bz1 + y * bz1 + z;
	}

	inline int getEdgeIdxZ(int x, int y, int z)
	{
		return x * by1 * bz + y * bz + z;
	}

	bool getSmoothing() { return bSmoothed; }
	void setSmoothing(bool _bSmooth) { bSmoothed = _bSmooth; }

	//void wipeIsoValues(float value = 0.f);

	void clear(); // deletes all the data. use whip

	void setGridPoints(float _x, float _y, float _z);

	// void updateTransformMatrix(){
	//	transform.makeScaleMatrix( scale );
	//	transform.rotate( orientation );
	//	transform.setTranslation( position );
	// }

	inline float MarchingCubes::getIsoValue(int x, int y, int z)
	{
		//return isoValsMorton[libmorton::morton3D_32_encode(x, y, z)];
		return isoVals[x * resY * resZ + y * resZ + z];
	}

	//inline float MarchingCubes::getIsoValueMorton(int x, int y, int z)
	//{
	//	return isoValsMorton[libmorton::morton3D_32_encode(x, y, z)];
	//}

	inline Vector3f getGridPoint(int x, int y, int z)
	{
		return gridPoints[x * resY * resZ + y * resZ + z];
	}

	void exportObj(string fileName);

	// private:
	// ofMatrix4x4 transform;
	// Number of grid points in each axis:
	int resX, resY, resZ;
	// Number of cubes in each axis, which is 1 less than number of grid points:
	int resXm1, resYm1, resZm1;

	int sx, sy, sz;
	int sx1, sy1, sz1;
	float dx, dy, dz;
	int bx, by, bz;
	int bx1, by1, bz1;
	int nbx, nby, nbz, nb;

	float flipNormalsValue;
	Vector3f cellDim;
	//vector<float> isoVals;
	vector<Vector3f> gridPoints;

	float* isoVals = nullptr;
	float* thresCmpArray = nullptr;
	int* thresCmpIntArray = nullptr;
	short* thresCmpShortArray = nullptr;
	float* edgeInterpVal = nullptr;
	Vector3f* verticesOffset = nullptr;

	vector<Vector3f> vertices;
	vector<Vector3f> normals;
	// ofVbo boundaryVbo;
	// ofVbo gridPointsVbo;
	int vertexCount, maxVertexCount;


	Vector3f vertList[12], normList[12];

	float threshold;
	bool bSmoothed, beenWarned;

	// transforms
	// Vector3f position, scale, up;
	// ofQuaternion orientation;
	//
	// ofMatrix3x3 normalMatrix;

	// ofVbo vbo;
	bool bUpdateMesh;

	// Blocking intermediate results
	int bX = 1, bY = 1, bZ = 1;  // Block size in number of cells.
	bool *thresCmp = nullptr;
	short *cubeIndices = nullptr;
	int *cubeIndicesInt = nullptr;
	Vector3f *bVertList = nullptr;

	// Exclude left boundary and include right boundary.
	Vector3f* interpVerticesX = nullptr;
	Vector3f* interpVerticesY = nullptr;
	Vector3f* interpVerticesZ = nullptr;


	// Level-by-level 
	bool* thresCmpLevel;
	int* cubeIndexLevel;
	Vector3f* vertInterpX;
	Vector3f* vertInterpY;
	Vector3f* vertInterpZ;
	float* isoValsMorton;

};
