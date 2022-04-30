
//http://paulbourke.net/geometry/polygonise/
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

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>

using namespace std;

enum MarchingCubesStage
{
	ALL,
	STAGE1,
	STAGE2,
	STAGE3,
};

class MarchingCubes{
public:
	MarchingCubes();
	~MarchingCubes();
	
	void setMaxVertexCount( int _maxVertexCount = 100000 );
	
	void setup( int resX=30, int resY=20, int resZ=30, int _maxVertexCount=150000);
	void update(){		update( threshold );}
	void update(float _threshold);
	
	//void draw( GLenum renderType = GL_TRIANGLES );
	//void drawWireframe();
	//void drawArrays( vector<Vector3f>* _vertices, vector<Vector3f>* _normals=NULL );
	//void drawGrid( bool drawGridPoints=true);
	
	void flipNormals(){	flipNormalsValue *= -1;}
	void setResolution( int _x=10, int _y=10, int _z=10 );
	void polygonise( int i, int j, int k );
	void computeNormal( int i, int j, int k );
	void vertexInterp(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, Vector3f& v, Vector3f& n);
	
	void setIsoValue( int x, int y, int z, float value);
	void addToIsoValue( int x, int y, int z, float value){
		getIsoValue(x,y,z) += value;
        bUpdateMesh = true;
	}
	
	bool getSmoothing(){	return bSmoothed;}
	void setSmoothing( bool _bSmooth ){		bSmoothed = _bSmooth;}
	
	
	void wipeIsoValues( float value=0.f);
	
	void clear();//deletes all the data. use whip
	
	
	void setGridPoints( float _x, float _y, float _z);
	
	//void updateTransformMatrix(){
	//	transform.makeScaleMatrix( scale );
	//	transform.rotate( orientation );
	//	transform.setTranslation( position );
	//}
	
	inline float& getIsoValue( int x, int y, int z){
		return isoVals[ x*resY*resZ+ y*resZ + z ];
	}
	inline Vector3f& getGridPoint( int x, int y, int z){
		return gridPoints[ x*resY*resZ+ y*resZ + z ];
	}
	inline Vector3f& getNormalVal( int x, int y, int z){
		return normalVals[ x*resY*resZ+ y*resZ + z ];
	}
	inline unsigned int& getGridPointComputed( int x, int y, int z){
		return gridPointComputed[ x*resY*resZ+ y*resZ + z ];
	}
	
	void exportObj( string fileName );
	
	//private:
	//ofMatrix4x4 transform;
	int	resX, resY, resZ;
	int resXm1, resYm1, resZm1;
	float flipNormalsValue;
	Vector3f cellDim;
	vector<float> isoVals;
	vector<Vector3f> gridPoints;
	vector<Vector3f> normalVals;
	vector<unsigned int> gridPointComputed;

	
	vector< Vector3f > vertices;
	vector< Vector3f > normals;
	//ofVbo boundaryVbo;
	//ofVbo gridPointsVbo;
	int vertexCount, maxVertexCount;
	
	Vector3f vertList[12], normList[12];
	
	float threshold;
	bool bSmoothed, beenWarned;
	
	//transforms
	//Vector3f position, scale, up;
	//ofQuaternion orientation;
	//
	//ofMatrix3x3 normalMatrix;
	
	//ofVbo vbo;
	bool bUpdateMesh;
};
