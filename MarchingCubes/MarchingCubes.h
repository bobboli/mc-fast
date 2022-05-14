
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
	~MarchingCubes()
	{
		if (thresCmp != nullptr) delete[] thresCmp;
		if (cubeIndices != nullptr) delete[] cubeIndices;
		if (bVertList != nullptr) delete[] bVertList;
	}

	void setMaxVertexCount(int _maxVertexCount = 100000);

	void setup(int resX = 30, int resY = 20, int resZ = 30, int _maxVertexCount = 200000);
	void setBlocking(int blockX, int blockY, int blockZ);
	void update() { update(threshold); }
	void update_block() { update_block(threshold); }
	void update(float _threshold);
	void update_block(float _threshold);
	void count_ops(float _threshold, operation_counts &counts);

	// void draw( GLenum renderType = GL_TRIANGLES );
	// void drawWireframe();
	// void drawArrays( vector<Vector3f>* _vertices, vector<Vector3f>* _normals=NULL );
	// void drawGrid( bool drawGridPoints=true);

	void flipNormals() { flipNormalsValue *= -1; }
	void setResolution(int _x = 10, int _y = 10, int _z = 10);
	void polygonise(int i, int j, int k);
	void polygonise_block(int i, int j, int k, int bX, int bY, int bZ);
	void polygonise_count_ops(int i, int j, int k, operation_counts &counts);
	void computeNormal(int i, int j, int k);
	void vertexInterp(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, Vector3f &v, Vector3f &n);
	void vertexInterp_count_ops(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, operation_counts &counts);

	void setIsoValue(int x, int y, int z, float value);
	void addToIsoValue(int x, int y, int z, float value)
	{
		getIsoValue(x, y, z) += value;
		bUpdateMesh = true;
	}

	bool getSmoothing() { return bSmoothed; }
	void setSmoothing(bool _bSmooth) { bSmoothed = _bSmooth; }

	void wipeIsoValues(float value = 0.f);

	void clear(); // deletes all the data. use whip

	void setGridPoints(float _x, float _y, float _z);

	// void updateTransformMatrix(){
	//	transform.makeScaleMatrix( scale );
	//	transform.rotate( orientation );
	//	transform.setTranslation( position );
	// }

	inline float &getIsoValue(int x, int y, int z)
	{
		return isoVals[x * resY * resZ + y * resZ + z];
	}
	inline Vector3f &getGridPoint(int x, int y, int z)
	{
		return gridPoints[x * resY * resZ + y * resZ + z];
	}
	inline Vector3f &getNormalVal(int x, int y, int z)
	{
		return normalVals[x * resY * resZ + y * resZ + z];
	}
	inline unsigned int &getGridPointComputed(int x, int y, int z)
	{
		return gridPointComputed[x * resY * resZ + y * resZ + z];
	}

	void exportObj(string fileName);

	// private:
	// ofMatrix4x4 transform;
	int resX, resY, resZ;
	int resXm1, resYm1, resZm1;
	float flipNormalsValue;
	Vector3f cellDim;
	vector<float> isoVals;
	vector<Vector3f> gridPoints;
	vector<Vector3f> normalVals;
	vector<unsigned int> gridPointComputed;

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
	int bX = 1, bY = 1, bZ = 1;
	bool *thresCmp;
	short *cubeIndices;
	Vector3f *bVertList;
};
