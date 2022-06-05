#include "MarchingCubes.h"
#include "libmorton/morton.h"

MarchingCubes::MarchingCubes(){
	threshold = .5;
	bSmoothed = true;
	flipNormalsValue = -1;	
}

MarchingCubes::~MarchingCubes() {
	if (isoVals != nullptr) delete[] isoVals;
	if (thresCmpArray != nullptr) delete[] thresCmpArray;
	if (thresCmpIntArray != nullptr) delete[] thresCmpIntArray;
	if (edgeInterpVal != nullptr) delete[] edgeInterpVal;
}

void MarchingCubes::setMaxVertexCount( int _maxVertexCount ){
	maxVertexCount = _maxVertexCount;
}


void MarchingCubes::setup(int resX, int resY, int resZ, int _maxVertexCount)
{
	clear();


	setResolution(resX, resY, resZ);
	setMaxVertexCount(_maxVertexCount);

	vertexCount = 0;


	// Trick: make them all same size so it's easy to calculate offsets
	vertIndexX = new int[sy1 * sz1 * 5];
	vertIndexY = vertIndexX + (sy1 * sz1);
	vertIndexZ = vertIndexY + (sy1 * sz1 * 2);
	int sz_xy = sy1 * sz1 * 3;
	offsetLookUp[1] = sy1 * sz1 * 2; 		offsetLookUp[13] = sy1 * sz1;
	offsetLookUp[2] = sz1;					offsetLookUp[14] = sz1;
	offsetLookUp[3] = sy1 * sz1; 			offsetLookUp[15] = sy1 * sz1 * 2;
	offsetLookUp[4] = 1; 					offsetLookUp[16] = 1;
	offsetLookUp[5] = sy1 * sz1 * 2 + 1;	offsetLookUp[17] = sy1 * sz1 + 1;
	offsetLookUp[6] = sz1 + 1;				offsetLookUp[18] = sz1 + 1;
	offsetLookUp[7] = sy1 * sz1 + 1;		offsetLookUp[19] = sy1 * sz1 * 2 + 1;
	offsetLookUp[8] = sz_xy;				offsetLookUp[20] = sz_xy + sy1 * sz1;
	offsetLookUp[9] = sz_xy + sy1 * sz1;	offsetLookUp[21] = sz_xy;
	offsetLookUp[10] = sz_xy + (sy1 + 1) * sz1; offsetLookUp[22] = sz_xy + sz1;
	offsetLookUp[11] = sz_xy + sz1;			offsetLookUp[23] = sz_xy + (sy1 + 1) * sz1;

	reset();
}

void MarchingCubes::reset()
{
	vertices.clear();
	indices.clear();
	int numEdges = resXm1 * (resYm1 + 1) * (resZm1 + 1) + resYm1 * (resXm1 + 1) * (resZm1 + 1) + resZm1 * (resXm1 + 1) * (resYm1 + 1);
	const long long int pctInterpEdge = 10;  // Might be overflow here.
	int numEdgesInterp = numEdges * pctInterpEdge / 100;

	vertices.reserve(numEdgesInterp);
	indices.reserve(4 * numEdgesInterp);
	//normals.reserve(numEdgesInterp);
}

void MarchingCubes::update(float _threshold){
	// update every time
			
		threshold = _threshold;
		
		//std::fill( normalVals.begin(), normalVals.end(), Vector3f());
		//std::fill( gridPointComputed.begin(), gridPointComputed.end(), 0 );
		vertexCount = 0;
		for (int x = 0; x < resX - 1; x++) {
			for (int y = 0; y < resY - 1; y++) {
				for (int z = 0; z < resZ - 1; z++) {
					polygonise(x, y, z);
				}
			}
		}
}


void MarchingCubes::update_level(float _threshold)
{
	threshold = _threshold;

	vertexCount = 0;

	for (int x = 0; x < sx; ++x)
	{
		polygonise_level(x);
	}
}



void MarchingCubes::count_ops(float _threshold, operation_counts& counts) {
	threshold = _threshold;

	vertexCount = 0;
	for (int x = 0; x < resX; x++) {
		for (int y = 0; y < resY; y++) {
			for (int z = 0; z < resZ; z++) {
				polygonise_count_ops(x, y, z, counts);
			}
		}
	}
}

//void ofxMarchingCubes::draw( GLenum renderType )
//{
//	ofPushMatrix();
//	ofMultMatrix( transform.getPtr() );
//	
//	vbo.draw( renderType, 0, vertexCount );
//
//	ofPopMatrix();
//}

void MarchingCubes::polygonise(int i, int j, int k){
	
	//if( vertexCount+3 < maxVertexCount ){
		bUpdateMesh = true;
		/*
		 Determine the index into the edge table which
		 tells us which vertices are inside of the surface
		 */
		int cubeindex = 0;
		int i1 = i + 1, j1 = j + 1, k1 = k + 1;
		cubeindex |= getIsoValue(i,j,k) > threshold ?   1 : 0;
		cubeindex |= getIsoValue(i1,j,k) > threshold ?   2 : 0;
		cubeindex |= getIsoValue(i1,j1,k) > threshold ?   4 : 0;
		cubeindex |= getIsoValue(i,j1,k) > threshold ?   8 : 0;
		cubeindex |= getIsoValue(i,j,k1) > threshold ?  16 : 0;
		cubeindex |= getIsoValue(i1,j,k1) > threshold ?  32 : 0;
		cubeindex |= getIsoValue(i1,j1,k1) > threshold ?  64 : 0;
		cubeindex |= getIsoValue(i,j1,k1) > threshold ? 128 : 0;
		
		/* Cube is entirely in/out of the surface */
		if (edgeTable[cubeindex] == 0)		return;
		
		/* Find the vertices where the surface intersects the cube */
		
		if (edgeTable[cubeindex] & 1)		vertexInterp(threshold, i,j,k, i1,j,k, vertList[0], normList[0]);
		if (edgeTable[cubeindex] & 2)		vertexInterp(threshold, i1,j,k, i1,j1,k, vertList[1], normList[1]);
		if (edgeTable[cubeindex] & 4)		vertexInterp(threshold, i1,j1,k, i,j1,k, vertList[2], normList[2]);
		if (edgeTable[cubeindex] & 8)		vertexInterp(threshold, i,j1,k, i,j,k, vertList[3], normList[3]);
		if (edgeTable[cubeindex] & 16)		vertexInterp(threshold, i,j,k1, i1,j,k1, vertList[4], normList[4]);
		if (edgeTable[cubeindex] & 32)		vertexInterp(threshold, i1,j,k1, i1,j1,k1, vertList[5], normList[5]);
		if (edgeTable[cubeindex] & 64)		vertexInterp(threshold, i1,j1,k1, i,j1,k1, vertList[6], normList[6]);
		if (edgeTable[cubeindex] & 128)		vertexInterp(threshold, i,j1,k1, i,j,k1, vertList[7], normList[7]);
		if (edgeTable[cubeindex] & 256)		vertexInterp(threshold, i,j,k, i,j,k1, vertList[8], normList[8]);
		if (edgeTable[cubeindex] & 512)		vertexInterp(threshold, i1,j,k, i1,j,k1, vertList[9], normList[9]);
		if (edgeTable[cubeindex] & 1024)	vertexInterp(threshold, i1,j1,k, i1,j1,k1, vertList[10], normList[10]);
		if (edgeTable[cubeindex] & 2048)	vertexInterp(threshold, i,j1,k, i,j1,k1, vertList[11], normList[11]);

		for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
			//if(bSmoothed){
			//	//smoothed normals
			//	normals[vertexCount] = normList[triTable[cubeindex][i]];
			//	normals[vertexCount+1] = normList[triTable[cubeindex][i+1]];
			//	normals[vertexCount+2] = normList[triTable[cubeindex][i+2]];
			//	
			//}
			//else{
			//	//faceted
			//	Vector3f a = vertList[triTable[cubeindex][i+1]] - vertList[triTable[cubeindex][i]];
			//	Vector3f b = vertList[triTable[cubeindex][i+2]] - vertList[triTable[cubeindex][i+1]];
			//	
			//	normals[vertexCount] = normals[vertexCount+1] = normals[vertexCount+2] = a.crossed(b).normalize() * flipNormalsValue;
			//}

			// Vector3f a = vertList[triTable[cubeindex][i + 1]] - vertList[triTable[cubeindex][i]];
			// Vector3f b = vertList[triTable[cubeindex][i+2]] - vertList[triTable[cubeindex][i+1]];

			//vertices[vertexCount++] = vertList[triTable[cubeindex][i]];
			//vertices[vertexCount++] = vertList[triTable[cubeindex][i+1]];
			//vertices[vertexCount++] = vertList[triTable[cubeindex][i+2]];

			vertices.push_back(vertList[triTable[cubeindex][i]]);
			vertices.push_back(vertList[triTable[cubeindex][i+1]]);
			vertices.push_back(vertList[triTable[cubeindex][i+2]]);
			vertexCount += 3;
		}
	//}
	//else if(!beenWarned){
	//	std::cerr << "ofxMarhingCubes: maximum vertex("+to_string(maxVertexCount)+") count exceded. try increasing the maxVertexCount with setMaxVertexCount()";
	//	beenWarned = true;
	//}
}



void MarchingCubes::polygonise_level(int level)
{

	bUpdateMesh = true;
	/*
	 Determine the index into the edge table which
	 tells us which vertices are inside of the surface
	 */

	bool* thresCmpOld, *thresCmpNew;

	int* vertIndexYOld, * vertIndexYNew;
	int* vertIndexZOld, * vertIndexZNew;

	// Trick: make them all same size so it's easy to calculate offset
	if (level % 2 == 0)
	{
		thresCmpOld = thresCmpLevel;
		thresCmpNew = thresCmpLevel + sy1 * sz1;

		vertIndexYOld = vertIndexY;
		vertIndexYNew = vertIndexY + sy1 * sz1;
		vertIndexZOld = vertIndexZ;
		vertIndexZNew = vertIndexZ + sy1 * sz1;
	}
	else
	{
		thresCmpNew = thresCmpLevel;
		thresCmpOld = thresCmpLevel + sy1 * sz1;

		vertIndexYNew = vertIndexY;
		vertIndexYOld = vertIndexY + sy1 * sz1;
		vertIndexZNew = vertIndexZ;
		vertIndexZOld = vertIndexZ + sy1 * sz1;
	}

	int x = level;
	int x1 = x + 1;

	// Threshold computing
	if (level == 0)
	{
		int iGrid = 0;
		for (int y = 0; y < sy1; ++y)
		{
			for (int z = 0; z < sz1; ++z)
			{
				thresCmpOld[iGrid] = getIsoValue(0, y, z) > threshold;//isoVals[iGrid] > threshold;
				++iGrid;
			}
		}
	}

	{
		int iGrid = 0;
		int iIsoVal = x1*sy1*sz1;
		for (int y = 0; y < sy1; ++y)
		{
			for (int z = 0; z < sz1; ++z)
			{
				thresCmpNew[iGrid] = getIsoValue(x1, y, z) > threshold;// isoVals[iIsoVal] > threshold;
				//if (level == 0)  // It seems that you would rather iterate over y and z once more like above, than using a nested if here
				//{
				//	thresCmpOld[iGrid] = getIsoValue(0, y, z) > threshold;//isoVals[iGrid] > threshold;
				//}
				++iGrid;
				++iIsoVal;
			}
		}
	}

	// Cube index computation
	{
		int iCube = 0;
		for (int y = 0; y < sy; ++y)
		{
			int y1 = y + 1;
			for (int z = 0; z < sz; ++z)
			{
				int z1 = z + 1;

				int cubeIndexOld = 0;
				int cubeIndexNew = 0;
				int base = y * sz1 + z;
				// todo: accessing of isoVals could be further optimized without using getIsoValue
				int thresCmpOld0 = thresCmpOld[base];
				int thresCmpNew0 = thresCmpNew[base];
				int thresCmpOld1 = thresCmpNew[base + sz1];
				int thresCmpNew1 = thresCmpOld[base + sz1];
				int thresCmpOld2 = thresCmpOld[base + 1];
				int thresCmpNew2 = thresCmpNew[base + 1];
				int thresCmpOld3 = thresCmpNew[base + sz1 + 1];
				int thresCmpNew3 = thresCmpOld[base + sz1 + 1];
				cubeIndexOld |= thresCmpOld0 | (thresCmpOld1 << 2) | (thresCmpOld2 << 4) | (thresCmpOld3 << 6);
				cubeIndexNew |= (thresCmpNew0 << 1) | (thresCmpNew1 << 3) | (thresCmpNew2 << 5) | (thresCmpNew3 << 7);
				cubeIndexLevel[iCube++] = cubeIndexOld | cubeIndexNew;
			}
		}
	}

	// Vertex interpolation
	Vector3f dummyN;
	Vector3f vert;  // todo: use multiple variables to increase ILP?


	// Y and Z edges
	if (level == 0)
	{
		for (int z = 0; z < sz; ++z)
		{
			if (edgeTable[cubeIndexLevel[z]] & 256)
			{
				vertIndexZOld[z] = vertices.size();
				vertexInterp_Z(threshold, 0, 0, z, z + 1, vert, dummyN);
				vertices.push_back(vert);
			}
		}

		for (int y = 0; y < sy; ++y)
		{
			if (edgeTable[cubeIndexLevel[y * sz]] & 8)
			{
				vertIndexYOld[y * sz1] = vertices.size();
				vertexInterp_Y(threshold, 0, y, y + 1, 0, vert, dummyN);
				vertices.push_back(vert);
			}
		}

		for (int y = 0; y < sy; ++y)
		{
			for (int z = 0; z < sz; ++z)
			{
				int edgeIndex = edgeTable[cubeIndexLevel[y * sz + z]];
				if (edgeIndex & 128)
				{
					vertIndexYOld[y * sz1 + (z + 1)] = vertices.size();
					vertexInterp_Y(threshold, 0, y, y + 1, z + 1, vert, dummyN);
					vertices.push_back(vert);
				}
				if (edgeIndex & 2048)
				{
					vertIndexZOld[(y + 1) * sz1 + z] = vertices.size();
					vertexInterp_Z(threshold, 0, y + 1, z, z + 1, vert, dummyN);
					vertices.push_back(vert);
				}
			}
		}
	}

	// X, Y and Z edges
	{
		{
			if (edgeTable[cubeIndexLevel[0]] & 1)
			{
				vertIndexX[0] = vertices.size();
				vertexInterp_X(threshold, x, x + 1, 0, 0, vert, dummyN);
				vertices.push_back(vert);
			}
		}

		for (int z = 0; z < sz; ++z)
		{
			if (edgeTable[cubeIndexLevel[z]] & 16)
			{
				vertIndexX[z + 1] = vertices.size();
				vertexInterp_X(threshold, x, x + 1, 0, z + 1, vert, dummyN);
				vertices.push_back(vert);	
			}
			if (edgeTable[cubeIndexLevel[z]] & 512)
			{
				vertIndexZNew[z] = vertices.size();
				vertexInterp_Z(threshold, x + 1, 0, z, z + 1, vert, dummyN);
				vertices.push_back(vert);
			}
		}


		for (int y = 0; y < sy; ++y)
		{
			if (edgeTable[cubeIndexLevel[y * sz]] & 4)
			{
				vertIndexX[(y + 1) * sz1] = vertices.size();
				vertexInterp_X(threshold, x, x + 1, y + 1, 0, vert, dummyN);
				vertices.push_back(vert);
			}
			if (edgeTable[cubeIndexLevel[y * sz]] & 2)
			{
				vertIndexYNew[y * sz1] = vertices.size();
				vertexInterp_Y(threshold, x + 1, y, y + 1, 0, vert, dummyN);
				vertices.push_back(vert);
			}
		}

		for (int y = 0; y < sy; ++y)
		{
			for (int z = 0; z < sz; ++z)
			{
				int iCube = y * sz + z;
				int cubeIndex = cubeIndexLevel[iCube];
				int edgeIndex = edgeTable[cubeIndex];
				if (edgeIndex & 64)
				{
					vertIndexX[(y + 1) * sz1 + (z + 1)] = vertices.size();
					vertexInterp_X(threshold, x, x + 1, y + 1, z + 1, vert, dummyN);
					vertices.push_back(vert);
				}
				if (edgeIndex & 32)
				{
					vertIndexYNew[y * sz1 + (z + 1)] = vertices.size();
					vertexInterp_Y(threshold, x + 1, y, y + 1, z + 1, vert, dummyN);
					vertices.push_back(vert);
				}
				if (edgeIndex & 1024)
				{
					vertIndexZNew[(y + 1) * sz1 + z] = vertices.size();
					vertexInterp_Z(threshold, x + 1, y + 1, z, z + 1, vert, dummyN);
					vertices.push_back(vert);
				}

				// Assembly triangles
				int base = iCube + y; // y * sz1 + z;
				int* triTableEntry = triTable[cubeIndex];
				for (int i = 0; triTableEntry[i] != -1; ++i)
				{
					int val = triTableEntry[i];
					int idx = base + offsetLookUp[val + (level & 1) * 12];
					int vertIndex = vertIndexX[idx];
					indices.push_back(vertIndex);
					++vertexCount;
				}
			}
		}
	}

	// Assembly triangles
	//for (int y = 0; y < sy; ++y)
	//{
	//	for (int z = 0; z < sz; ++z)
	//	{
	//		int iCube = y * sz + z;
	//		int base = iCube + y; // y * sz1 + z;
	//		int *triTableEntry = triTable[cubeIndexLevel[iCube]];
	//		for (int ti = 0; triTableEntry[ti] != -1; ti += 3)
	//		{
	//			for (int tj = 0; tj < 3; tj++)
	//			{
	//				int val = triTable[cubeIndexLevel[iCube]][ti + tj];
	//				int idx = base + offsetLookUp[val + (level & 1) * 12];
	//				int vertIndex = vertIndexX[idx];
	//				indices.push_back(vertIndex);
	//				++vertexCount;
	//			}
	//		}
	//	}
	//}
}


void MarchingCubes::polygonise_level_vec(int level)
{

	bUpdateMesh = true;
	/*
	 Determine the index into the edge table which
	 tells us which vertices are inside of the surface
	 */

	int* thresCmpOld, * thresCmpNew;

	int* vertIndexYOld, * vertIndexYNew;
	int* vertIndexZOld, * vertIndexZNew;

	// Trick: make them all same size so it's easy to calculate offset
	if (level % 2 == 0)
	{
		thresCmpOld = thresCmpLevelInt;
		thresCmpNew = thresCmpLevelInt + sy1 * sz1;

		vertIndexYOld = vertIndexY;
		vertIndexYNew = vertIndexY + sy1 * sz1;
		vertIndexZOld = vertIndexZ;
		vertIndexZNew = vertIndexZ + sy1 * sz1;
	}
	else
	{
		thresCmpNew = thresCmpLevelInt;
		thresCmpOld = thresCmpLevelInt + sy1 * sz1;

		vertIndexYNew = vertIndexY;
		vertIndexYOld = vertIndexY + sy1 * sz1;
		vertIndexZNew = vertIndexZ;
		vertIndexZOld = vertIndexZ + sy1 * sz1;
	}

	int x = level;
	int x1 = x + 1;

	__m256 vt = _mm256_set1_ps(threshold);
	__m256 c1 = _mm256_set1_ps(1);
	int DX = sz1 * sy1, DY = sz1;

	// Threshold computing
	if (level == 0)
	{
		for (int y = 0; y < sy1; ++y)
		{
			int z;
			for (z = 0; z < sz1 - 7; z += 8)
			{
				int idx = y * DY + z;
				__m256 vals = _mm256_loadu_ps(isoVals + idx);
				__m256 cmp = _mm256_cmp_ps(vals, vt, _CMP_GT_OQ);
				__m256 res = _mm256_and_ps(cmp, c1);
				_mm256_storeu_epi32(thresCmpOld + idx, _mm256_cvtps_epi32(res));
			}
			for (; z < sz1; ++z)
			{
				int idx = y * DY + z;
				thresCmpOld[idx] = getIsoValue(0, y, z) > threshold;//isoVals[iGrid] > threshold;
			}
		}
	}

	{
		int base = x1 * DX;
		for (int y = 0; y < sy1; ++y)
		{
			int z;
			for (z = 0; z < sz1 - 7; z += 8)
			{
				int idx = base + y * DY + z;
				int level_idx = y * DY + z;

				__m256 vals = _mm256_loadu_ps(isoVals + idx);
				__m256 cmp = _mm256_cmp_ps(vals, vt, _CMP_GT_OQ);
				__m256 res = _mm256_and_ps(cmp, c1);
				_mm256_storeu_epi32(thresCmpNew + level_idx, _mm256_cvtps_epi32(res));
			}
			for (; z < sz1; ++z)
			{
				int idx = y * DY + z;
				thresCmpNew[idx] = getIsoValue(x1, y, z) > threshold;//isoVals[iGrid] > threshold;
			}
		}
	}

	// Cube index computation
	{
		for (int y = 0; y < sy; ++y)
		{
			int y1 = y + 1;
			int z = 0;
			for (z = 0; z < sz - 7; z += 8)
			{
				int base = y * DY + z;
				int idx = y * sz + z;
				__m256i b_000 = _mm256_loadu_epi32(thresCmpOld + base);
				__m256i b_001 = _mm256_loadu_epi32(thresCmpOld + base + 1);
				__m256i b_010 = _mm256_loadu_epi32(thresCmpOld + base + DY);
				__m256i b_011 = _mm256_loadu_epi32(thresCmpOld + base + DY + 1);
				__m256i b_100 = _mm256_loadu_epi32(thresCmpNew + base);
				__m256i b_101 = _mm256_loadu_epi32(thresCmpNew + base + 1);
				__m256i b_110 = _mm256_loadu_epi32(thresCmpNew + base + DY);
				__m256i b_111 = _mm256_loadu_epi32(thresCmpNew + base + DY + 1);
				__m256i bs_100 = _mm256_slli_epi32(b_100, 1);
				__m256i bs_110 = _mm256_slli_epi32(b_110, 2);
				__m256i bs_010 = _mm256_slli_epi32(b_010, 3);
				__m256i bs_001 = _mm256_slli_epi32(b_001, 4);
				__m256i bs_101 = _mm256_slli_epi32(b_101, 5);
				__m256i bs_111 = _mm256_slli_epi32(b_111, 6);
				__m256i bs_011 = _mm256_slli_epi32(b_011, 7);
				__m256i cube_index = _mm256_set1_epi32(0);
				cube_index = _mm256_add_epi32(cube_index, b_000);
				cube_index = _mm256_add_epi32(cube_index, bs_001);
				cube_index = _mm256_add_epi32(cube_index, bs_010);
				cube_index = _mm256_add_epi32(cube_index, bs_011);
				cube_index = _mm256_add_epi32(cube_index, bs_100);
				cube_index = _mm256_add_epi32(cube_index, bs_101);
				cube_index = _mm256_add_epi32(cube_index, bs_110);
				cube_index = _mm256_add_epi32(cube_index, bs_111);
				_mm256_storeu_epi32(cubeIndexLevel + idx, cube_index);
			}

			for(; z < sz; ++z)
			{
				int z1 = z + 1;
				int idx = sz * y + z;
				int cubeIndexOld = 0;
				int cubeIndexNew = 0;
				int base = y * sz1 + z;
				// todo: accessing of isoVals could be further optimized without using getIsoValue
				int thresCmpOld0 = thresCmpOld[base];
				int thresCmpNew0 = thresCmpNew[base];
				int thresCmpOld1 = thresCmpNew[base + sz1];
				int thresCmpNew1 = thresCmpOld[base + sz1];
				int thresCmpOld2 = thresCmpOld[base + 1];
				int thresCmpNew2 = thresCmpNew[base + 1];
				int thresCmpOld3 = thresCmpNew[base + sz1 + 1];
				int thresCmpNew3 = thresCmpOld[base + sz1 + 1];
				cubeIndexOld |= thresCmpOld0 | (thresCmpOld1 << 2) | (thresCmpOld2 << 4) | (thresCmpOld3 << 6);
				cubeIndexNew |= (thresCmpNew0 << 1) | (thresCmpNew1 << 3) | (thresCmpNew2 << 5) | (thresCmpNew3 << 7);
				cubeIndexLevel[idx] = cubeIndexOld | cubeIndexNew;
			}
		}
	}

	// Vertex interpolation
	Vector3f dummyN;
	Vector3f vert;  // todo: use multiple variables to increase ILP?


	// Y and Z edges
	if (level == 0)
	{
		for (int z = 0; z < sz; ++z)
		{
			if (edgeTable[cubeIndexLevel[z]] & 256)
			{
				vertIndexZOld[z] = vertices.size();
				vertexInterp_Z(threshold, 0, 0, z, z + 1, vert, dummyN);
				vertices.push_back(vert);
			}
		}

		for (int y = 0; y < sy; ++y)
		{
			if (edgeTable[cubeIndexLevel[y * sz]] & 8)
			{
				vertIndexYOld[y * sz1] = vertices.size();
				vertexInterp_Y(threshold, 0, y, y + 1, 0, vert, dummyN);
				vertices.push_back(vert);
			}
		}

		for (int y = 0; y < sy; ++y)
		{
			for (int z = 0; z < sz; ++z)
			{
				int edgeIndex = edgeTable[cubeIndexLevel[y * sz + z]];
				if (edgeIndex & 128)
				{
					vertIndexYOld[y * sz1 + (z + 1)] = vertices.size();
					vertexInterp_Y(threshold, 0, y, y + 1, z + 1, vert, dummyN);
					vertices.push_back(vert);
				}
				if (edgeIndex & 2048)
				{
					vertIndexZOld[(y + 1) * sz1 + z] = vertices.size();
					vertexInterp_Z(threshold, 0, y + 1, z, z + 1, vert, dummyN);
					vertices.push_back(vert);
				}
			}
		}
	}

	// X, Y and Z edges
	{
		{
			if (edgeTable[cubeIndexLevel[0]] & 1)
			{
				vertIndexX[0] = vertices.size();
				vertexInterp_X(threshold, x, x + 1, 0, 0, vert, dummyN);
				vertices.push_back(vert);
			}
		}

		for (int z = 0; z < sz; ++z)
		{
			if (edgeTable[cubeIndexLevel[z]] & 16)
			{
				vertIndexX[z + 1] = vertices.size();
				vertexInterp_X(threshold, x, x + 1, 0, z + 1, vert, dummyN);
				vertices.push_back(vert);
			}
			if (edgeTable[cubeIndexLevel[z]] & 512)
			{
				vertIndexZNew[z] = vertices.size();
				vertexInterp_Z(threshold, x + 1, 0, z, z + 1, vert, dummyN);
				vertices.push_back(vert);
			}
		}


		for (int y = 0; y < sy; ++y)
		{
			if (edgeTable[cubeIndexLevel[y * sz]] & 4)
			{
				vertIndexX[(y + 1) * sz1] = vertices.size();
				vertexInterp_X(threshold, x, x + 1, y + 1, 0, vert, dummyN);
				vertices.push_back(vert);
			}
			if (edgeTable[cubeIndexLevel[y * sz]] & 2)
			{
				vertIndexYNew[y * sz1] = vertices.size();
				vertexInterp_Y(threshold, x + 1, y, y + 1, 0, vert, dummyN);
				vertices.push_back(vert);
			}
		}

		int numEdgeX = 0, numEdgeY = 0, numEdgeZ = 0;
		int curNumVertices = vertices.size();
		int oldNumVertices = curNumVertices;

		// First pass: Generate vertex indices. Generate a list of edges to be interpolated (CSR-like format).
		// Does not really interpolate vertices.
		
		for (int y = 0; y < sy; ++y)
		{
			yStart_EdgeX[y] = numEdgeX;
			yStart_EdgeY[y] = numEdgeY;
			yStart_EdgeZ[y] = numEdgeZ;

			{
				int iCube = y * sz;
				for (int z = 0; z < sz; ++z)
				{
					edgeIndexArray[z] = edgeTable[cubeIndexLevel[iCube++]];
				}
			}

			// todo: this step (especially & operation) should also be vectorized.
			for (int z = 0; z < sz; ++z)
			{
				if (edgeIndexArray[z] & 64)
				{
					vertIndexX[(y + 1) * sz1 + (z + 1)] = curNumVertices++;
					zIndex_EdgeX[numEdgeX++] = z;
					//vertexInterp_X(threshold, x, x + 1, y + 1, z + 1, vert, dummyN);
					//vertices.push_back(vert);
				}
			}

			for (int z = 0; z < sz; ++z)
			{
				if (edgeIndexArray[z] & 32)
				{
					vertIndexYNew[y * sz1 + (z + 1)] = curNumVertices++;
					zIndex_EdgeY[numEdgeY++] = z;
					//vertexInterp_Y(threshold, x + 1, y, y + 1, z + 1, vert, dummyN);
					//vertices.push_back(vert);
				}
			}

			for (int z = 0; z < sz; ++z)
			{
				if (edgeIndexArray[z] & 1024)
				{
					vertIndexZNew[(y + 1) * sz1 + z] = curNumVertices++;
					zIndex_EdgeZ[numEdgeZ++] = z;
					//vertexInterp_Z(threshold, x + 1, y + 1, z, z + 1, vert, dummyN);
					//vertices.push_back(vert);
				}

				// Vertex indices could already be assembled at this time:
				// (Though a little bit strange to put it here).
				{
					int iCube = y * sz + z;
					int base = iCube + y;
					int cubeIndex = cubeIndexLevel[iCube];
					int* triTableEntry = triTable[cubeIndex];
					for (int i = 0; triTableEntry[i] != -1; ++i)
					{
						int val = triTableEntry[i];
						int idx = base + offsetLookUp[val + (level & 1) * 12];
						int vertIndex = vertIndexX[idx];
						indices.push_back(vertIndex);
						++vertexCount;
					}
				}


			}
		}
		yStart_EdgeX[sy] = numEdgeX;
		yStart_EdgeY[sy] = numEdgeY;
		yStart_EdgeZ[sy] = numEdgeZ;



		// Second pass: Interpolate vertices and assemble triangles.
		vertices.resize(curNumVertices);

	#if 0
		// Non-SIMD vertex interpolation.
		for (int y = 0; y < sy; ++y)
		{

			// Interpolate vertices on Edge X. 
			for (int i = yStart_EdgeX[y]; i < yStart_EdgeX[y + 1]; ++i)
			{
				int z = zIndex_EdgeX[i];  // All z's that should be interpolated
				vertexInterp_X(threshold, x, x + 1, y + 1, z + 1, vert, dummyN);
				vertices[oldNumVertices++] = vert;
			}

			// Interpolate vertices on Edge Y. 
			for (int i = yStart_EdgeY[y]; i < yStart_EdgeY[y + 1]; ++i)
			{
				int z = zIndex_EdgeY[i];  // All z's that should be interpolated
				vertexInterp_Y(threshold, x + 1, y, y + 1, z + 1, vert, dummyN);
				vertices[oldNumVertices++] = vert;
			}

			// Interpolate vertices on Edge X. 
			for (int i = yStart_EdgeZ[y]; i < yStart_EdgeZ[y + 1]; ++i)
			{
				int z = zIndex_EdgeZ[i];  // All z's that should be interpolated
				vertexInterp_Z(threshold, x + 1, y + 1, z, z + 1, vert, dummyN);
				vertices[oldNumVertices++] = vert;
			}
		}


	#else
		__m256 cdx = _mm256_set1_ps(dx);
		__m256 cdy = _mm256_set1_ps(dy);
		__m256 cdz = _mm256_set1_ps(dz);

		__m256 cx = _mm256_set1_ps(x);
		__m256 p1x = _mm256_mul_ps(cdx, cx);
		__m256 c0 = _mm256_set1_ps(0);
		__m256 c1 = _mm256_set1_ps(1);
		__m256i c1i = _mm256_set1_epi32(1);
		float vert_x1 = x * dx;
		float vert_x2 = vert_x1 + dx;
		for (int y = 0; y < sy; ++y)
		{
			float vert_y1 = dy * y;
			float vert_y2 = vert_y1 + dy;
			int basex = x * DX + (y+1) * DY;
			int basey = (x+1) * DX + y * DY;
			int basez = (x+1) * DX + (y+1) * DY;
			__m256 p1y = _mm256_mul_ps(cdy, _mm256_set1_ps(y));

			// Interpolate vertices on Edge X. 
			int i = yStart_EdgeX[y];
			for (i = yStart_EdgeX[y]; i < yStart_EdgeX[y + 1] - 7; i += 8)
			{
				__m256i idx = _mm256_loadu_epi32(zIndex_EdgeX + i); // z
				idx = _mm256_add_epi32(idx, c1i);  // z+1
				__m256 val_start = _mm256_i32gather_ps(isoVals + basex, idx, 4);
				__m256 val_end = _mm256_i32gather_ps(isoVals + basex + DX, idx, 4);
				__m256 denominator = _mm256_sub_ps(val_end, val_start);
				__m256 numerator = _mm256_sub_ps(vt, val_start);
				__m256 tmp = _mm256_sub_ps(vt, val_end);
				__m256 lerp = _mm256_div_ps(numerator, denominator);
				__m256 dp = _mm256_mul_ps(cdx, lerp);
				__m256 px = _mm256_add_ps(p1x, dp);

				vert.y = vert_y2;
				for (int j = 0; j < 8; ++j)
				{
					vert.z = (idx.m256i_i32[j]) * dz;
					// there is no m256 abs ps
					if (abs(denominator.m256_f32[j]) < 0.0001)
					{
						vert.x = vert_x1;
						vertices[oldNumVertices++] = vert;
						continue;
					}
					if (abs(numerator.m256_f32[j]) < 0.0001)
					{
						vert.x = vert_x1;
						vertices[oldNumVertices++] = vert;
						continue;
					}
					if (abs(tmp.m256_f32[j]) < 0.0001)
					{
						vert.x = vert_x2;
						vertices[oldNumVertices++] = vert;
						continue;
					}
										
					vert.x = px.m256_f32[j];
					vertices[oldNumVertices++] = vert;
				}
			}
			for (; i < yStart_EdgeX[y + 1]; i++)
			{
				int z = zIndex_EdgeX[i];  // All z's that should be interpolated
				vertexInterp_X(threshold, x, x + 1, y + 1, z + 1, vert, dummyN);
				vertices[oldNumVertices++] = vert;
			}

			// Interpolate vertices on Edge Y. 
			for (i = yStart_EdgeY[y]; i < yStart_EdgeY[y + 1] - 7; i += 8)
			{
				__m256i idx = _mm256_loadu_epi32(zIndex_EdgeY + i); // z
				idx = _mm256_add_epi32(idx, c1i);  // z+1
				__m256 val_start = _mm256_i32gather_ps(isoVals + basey, idx, 4); // x+1 y z+1
				__m256 val_end = _mm256_i32gather_ps(isoVals + basey + DY, idx, 4); // x+1 y+1 z+1
				__m256 denominator = _mm256_sub_ps(val_end, val_start);
				__m256 numerator = _mm256_sub_ps(vt, val_start);
				__m256 tmp = _mm256_sub_ps(vt, val_end);
				__m256 lerp = _mm256_div_ps(numerator, denominator);
				__m256 dp = _mm256_mul_ps(cdy, lerp);
				__m256 py = _mm256_add_ps(p1y, dp);

				vert.x = vert_x2;
				for (int j = 0; j < 8; ++j)
				{
					vert.z = (idx.m256i_i32[j]) * dz;
					// there is no m256 abs ps
					if (abs(denominator.m256_f32[j]) < 0.0001)
					{
						vert.y = vert_y1;
						vertices[oldNumVertices++] = vert;
						continue;
					}
					if (abs(numerator.m256_f32[j]) < 0.0001)
					{
						vert.y = vert_y1;
						vertices[oldNumVertices++] = vert;
						continue;
					}
					if (abs(tmp.m256_f32[j]) < 0.0001)
					{
						vert.y = vert_y2;
						vertices[oldNumVertices++] = vert;
						continue;
					}

					vert.y = py.m256_f32[j];
					vertices[oldNumVertices++] = vert;
				}
			}

			for (; i < yStart_EdgeY[y + 1]; ++i)
			{
				int z = zIndex_EdgeY[i];
				vertexInterp_Y(threshold, x + 1, y, y + 1, z + 1, vert, dummyN);
				vertices[oldNumVertices++] = vert;
			}

			// Interpolate vertices on Edge Z. 
			for (i = yStart_EdgeZ[y]; i < yStart_EdgeZ[y + 1] - 7; i += 8)
			{
				__m256i idx = _mm256_loadu_epi32(zIndex_EdgeZ + i); // z
				__m256i idx1 = _mm256_add_epi32(idx, c1i);  // z+1
				__m256 p1z = _mm256_mul_ps(cdz, _mm256_cvtepi32_ps(idx));
				__m256 p2z = _mm256_mul_ps(cdz, _mm256_cvtepi32_ps(idx1));
				__m256 val_start = _mm256_i32gather_ps(isoVals + basez, idx, 4); // x+1 y+1 z
				__m256 val_end = _mm256_i32gather_ps(isoVals + basez, idx1, 4); // x+1 y+1 z+1
				__m256 denominator = _mm256_sub_ps(val_end, val_start);
				__m256 numerator = _mm256_sub_ps(vt, val_start);
				__m256 tmp = _mm256_sub_ps(vt, val_end);
				__m256 lerp = _mm256_div_ps(numerator, denominator);
				__m256 dp = _mm256_mul_ps(cdz, lerp);
				__m256 pz = _mm256_add_ps(p1z, dp);

				vert.x = vert_x2;
				vert.y = vert_y2;
				for (int j = 0; j < 8; ++j)
				{
					// there is no m256 abs ps
					if (abs(denominator.m256_f32[j]) < 0.0001)
					{
						vert.z = p1z.m256_f32[j];
						vertices[oldNumVertices++] = vert;
						continue;
					}
					if (abs(numerator.m256_f32[j]) < 0.0001)
					{
						vert.z = p1z.m256_f32[j];
						vertices[oldNumVertices++] = vert;
						continue;
					}
					if (abs(tmp.m256_f32[j]) < 0.0001)
					{
						vert.z = p2z.m256_f32[j];
						vertices[oldNumVertices++] = vert;
						continue;
					}

					vert.z = pz.m256_f32[j];
					vertices[oldNumVertices++] = vert;
				}
			}

			for (; i < yStart_EdgeZ[y + 1]; ++i)
			{
				int z = zIndex_EdgeZ[i];
				vertexInterp_Z(threshold, x + 1, y + 1, z, z + 1, vert, dummyN);
				vertices[oldNumVertices++] = vert;
			}
		}

	#endif
	}
}

void MarchingCubes::polygonise_count_ops(int i, int j, int k, operation_counts& counts)
{

	int cubeindex = 0;
	int i1 = min(i + 1, resXm1), j1 = min(j + 1, resYm1), k1 = min(k + 1, resZm1);
	cubeindex |= getIsoValue(i, j, k) > threshold ? 1 : 0;
	cubeindex |= getIsoValue(i1, j, k) > threshold ? 2 : 0;
	cubeindex |= getIsoValue(i1, j1, k) > threshold ? 4 : 0;
	cubeindex |= getIsoValue(i, j1, k) > threshold ? 8 : 0;
	cubeindex |= getIsoValue(i, j, k1) > threshold ? 16 : 0;
	cubeindex |= getIsoValue(i1, j, k1) > threshold ? 32 : 0;
	cubeindex |= getIsoValue(i1, j1, k1) > threshold ? 64 : 0;
	cubeindex |= getIsoValue(i, j1, k1) > threshold ? 128 : 0;
	// 8 float cmps
	counts.fl_cmp += 8;

	if (edgeTable[cubeindex] == 0)
	{
		counts.empty_cells += 1;
		return;
	}

	if (edgeTable[cubeindex] & 1)		vertexInterp_count_ops(threshold, i, j, k, i1, j, k, counts);
	if (edgeTable[cubeindex] & 2)		vertexInterp_count_ops(threshold, i1, j, k, i1, j1, k, counts);
	if (edgeTable[cubeindex] & 4)		vertexInterp_count_ops(threshold, i1, j1, k, i, j1, k, counts);
	if (edgeTable[cubeindex] & 8)		vertexInterp_count_ops(threshold, i, j1, k, i, j, k, counts);
	if (edgeTable[cubeindex] & 16)		vertexInterp_count_ops(threshold, i, j, k1, i1, j, k1, counts);
	if (edgeTable[cubeindex] & 32)		vertexInterp_count_ops(threshold, i1, j, k1, i1, j1, k1, counts);
	if (edgeTable[cubeindex] & 64)		vertexInterp_count_ops(threshold, i1, j1, k1, i, j1, k1, counts);
	if (edgeTable[cubeindex] & 128)		vertexInterp_count_ops(threshold, i, j1, k1, i, j, k1, counts);
	if (edgeTable[cubeindex] & 256)		vertexInterp_count_ops(threshold, i, j, k, i, j, k1, counts);
	if (edgeTable[cubeindex] & 512)		vertexInterp_count_ops(threshold, i1, j, k, i1, j, k1, counts);
	if (edgeTable[cubeindex] & 1024)	vertexInterp_count_ops(threshold, i1, j1, k, i1, j1, k1, counts);
	if (edgeTable[cubeindex] & 2048)	vertexInterp_count_ops(threshold, i, j1, k, i, j1, k1, counts);

	for (i = 0; triTable[cubeindex][i] != -1; i += 3)
	{
		// Vector3f a = vertList[triTable[cubeindex][i + 1]] - vertList[triTable[cubeindex][i]];
		// Vector3f b = vertList[triTable[cubeindex][i + 2]] - vertList[triTable[cubeindex][i + 1]];

		//counts.fl_add += 6;

		/*vertices[vertexCount] = vertList[triTable[cubeindex][i]];
		vertices[vertexCount + 1] = vertList[triTable[cubeindex][i + 1]];
		vertices[vertexCount + 2] = vertList[triTable[cubeindex][i + 2]];*/

		vertexCount += 3;
	}
}

inline void MarchingCubes::vertexInterp(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, Vector3f& v, Vector3f& n){
	
	Vector3f p1 = Vector3f(i1 * dx, j1 * dy, k1 * dz);
	Vector3f p2 = Vector3f(i2 * dx, j2 * dy, k2 * dz);
	
	float iso1 = getIsoValue(i1,j1,k1);
	float iso2 = getIsoValue(i2,j2,k2);
	
	//if(bSmoothed){
	//	//we need to interpolate/calculate the normals
	//	computeNormal( i1, j1, k1 );
	//	computeNormal( i2, j2, k2 );
	//	
	//	Vector3f& n1 = getNormalVal(i1,j1,k1);
	//	Vector3f& n2 = getNormalVal(i2,j2,k2);
	//	
	//	if (abs(threshold-iso1) < 0.00001){
	//		v = p1;
	//		n = n1;
	//		return;
	//	}
	//	if (abs(threshold-iso2) < 0.00001){
	//		v = p2;
	//		n = n2;
	//		return;
	//	}
	//	if (abs(iso1-iso2) < 0.00001){
	//		v = p1;
	//		n = n1;
	//		return;
	//	}
	//	
	//	//lerp
	//	float t = (threshold - iso1) / (iso2 - iso1);
	//	v = p1 + (p2-p1) * t;
	//	n = n1 + (n2-n1) * t;
	//}
	//
	//else{
	//	//we'll calc the normal later
	//	if (abs(threshold-iso1) < 0.00001){
	//		v = p1;
	//		return;
	//	}
	//	if (abs(threshold-iso2) < 0.00001){
	//		v = p2;
	//		return;
	//	}
	//	if (abs(iso1-iso2) < 0.00001){
	//		v = p1;
	//		return;
	//	}
	//	
	//	//lerp
	//	v = p1 + (p2-p1) * (threshold - iso1) / (iso2 - iso1);
	//}
	if (abs(threshold - iso1) < 0.00001) {
		v = p1;
		return;
	}
	if (abs(threshold-iso2) < 0.00001){
		v = p2;
		return;
	}
	if (abs(iso1-iso2) < 0.00001){
		v = p1;
		return;
	}
			
	//lerp
	v = p1 + (p2-p1) * ((threshold - iso1) / (iso2 - iso1));
}

void MarchingCubes::vertexInterp_X(float threshold, int x1, int x2, int y, int z, Vector3f& v, Vector3f& n)
{
	v.y = y * dy;
	v.z = z * dz;
	float p1 = x1 * dx, p2 = x2 * dx;
	float iso1 = getIsoValue(x1, y, z), iso2 = getIsoValue(x2, y, z);
	if (abs(threshold - iso1) < 0.00001)
	{
		v.x = p1;
		return;
	}
	if (abs(threshold - iso2) < 0.00001)
	{
		v.x = p2;
		return;
	}
	if (abs(iso1 - iso2) < 0.00001)
	{
		v.x = p1;
		return;
	}

	//lerp
	v.x = p1 + (p2 - p1) * ((threshold - iso1) / (iso2 - iso1));
}

void MarchingCubes::vertexInterp_Y(float threshold, int x, int y1, int y2, int z, Vector3f& v, Vector3f& n)
{
	v.x = x * dx;
	v.z = z * dz;
	float p1 = y1 * dy, p2 = y2 * dy;
	float iso1 = getIsoValue(x, y1, z), iso2 = getIsoValue(x, y2, z);
	if (abs(threshold - iso1) < 0.00001)
	{
		v.y = p1;
		return;
	}
	if (abs(threshold - iso2) < 0.00001)
	{
		v.y = p2;
		return;
	}
	if (abs(iso1 - iso2) < 0.00001)
	{
		v.y = p1;
		return;
	}

	//lerp
	v.y = p1 + (p2 - p1) * ((threshold - iso1) / (iso2 - iso1));
}

void MarchingCubes::vertexInterp_Z(float threshold, int x, int y, int z1, int z2, Vector3f& v, Vector3f& n)
{
	v.x = x * dx;
	v.y = y * dy;
	float p1 = z1 * dz, p2 = z2 * dz;
	float iso1 = getIsoValue(x, y, z1), iso2 = getIsoValue(x, y, z2);
	if (abs(threshold - iso1) < 0.00001)
	{
		v.z = p1;
		return;
	}
	if (abs(threshold - iso2) < 0.00001)
	{
		v.z = p2;
		return;
	}
	if (abs(iso1 - iso2) < 0.00001)
	{
		v.z = p1;
		return;
	}

	//lerp
	v.z = p1 + (p2 - p1) * ((threshold - iso1) / (iso2 - iso1));
}

void MarchingCubes::vertexInterp_vec(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, Vector3f& v, Vector3f& n) {

	Vector3f p1 = Vector3f(i1 * dx, j1 * dy, k1 * dz);
	Vector3f p2 = Vector3f(i2 * dx, j2 * dy, k2 * dz);
	float iso1 = getIsoValue(i1, j1, k1);
	float iso2 = getIsoValue(i2, j2, k2);

	int idx1 = i1 * resZ * resY + j1 * resZ + k1;
	int idx2 = i2 * resZ * resY + j2 * resZ + k2;
	
	if (abs(threshold - iso1) < 0.00001) {
		v = p1;
		return;
	}
	if (abs(threshold - iso2) < 0.00001) {
		v = p2;
		return;
	}
	if (abs(iso1 - iso2) < 0.00001) {
		v = p1;
		return;
	}
	int num = resX * resY * resZ;

	if (i2 == i1 + 1)
	{
		//lerp
		v = p1 + (p2 - p1) * edgeInterpVal[idx1];
		return;
	}
	else if (i2 == i1 - 1)
	{
		//lerp
		v = p2 + (p1 - p2) * edgeInterpVal[idx2];
		return;
	}
	else if (j2 == j1 + 1)
	{
		v = p1 + (p2 - p1) * edgeInterpVal[idx1 + num];
		return;
	}
	else if (j2 == j1 - 1) 
	{
		v = p2 + (p1 - p2) * edgeInterpVal[idx2 + num];
		return;
	}
	else if(k2 == k1 + 1)
	{
		v = p1 + (p2 - p1) * edgeInterpVal[idx1 + num * 2];
		return;
	}
	else if (k2 == k1 - 1)
	{
		v = p2 + (p1 - p2) * edgeInterpVal[idx2 + num * 2];
		return;
	}
}

void MarchingCubes::vertexInterp_count_ops(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, operation_counts& counts) {

	Vector3f p1 = Vector3f(i1 * dx, j1 * dy, k1 * dz);
	Vector3f p2 = Vector3f(i2 * dx, j2 * dy, k2 * dz);

	float iso1 = getIsoValue(i1, j1, k1);
	float iso2 = getIsoValue(i2, j2, k2);

	// 1 float add, 1 fabs, 1 fcmp
	counts.fl_cmp += 1;
	counts.fl_add += 1;
	counts.fl_abs += 1;
	if (abs(threshold - iso1) < 0.00001) {
		return;
	}
	counts.fl_cmp += 1;
	counts.fl_add += 1;
	counts.fl_abs += 1;
	if (abs(threshold - iso2) < 0.00001) {
		return;
	}
	counts.fl_cmp += 1;
	counts.fl_add += 1;
	counts.fl_abs += 1;
	if (abs(iso1 - iso2) < 0.00001) {
		return;
	}

	//lerp 
	// 3*2+2=8 adds, 1 div, 3 muls
	//v = p1 + (p2 - p1) * ((threshold - iso1) / (iso2 - iso1));
	counts.fl_add += 8;
	counts.fl_div += 1;
	counts.fl_mul += 3;
}

void MarchingCubes::computeNormal( int i, int j, int k ) {
	

		//Vector3f& n = getNormalVal(i, j, k);// normalVals[i][j][k];
		//n.set(getIsoValue(min(resXm1, i+1), j, k) - getIsoValue(max(0,i-1),j,k),
		//	  getIsoValue(i,min(resYm1, j+1),k) - getIsoValue(i,max(0,j-1),k),
		//	  getIsoValue(i,j,min(resZm1, k+1)) - getIsoValue(i,j,max(0,k-1)));
		//
		//n.normalize();
		//n = n*flipNormalsValue;
		//getGridPointComputed(i,j,k) = 1;
};


void MarchingCubes::setIsoValue( int x, int y, int z, float value){
	isoVals[min(sx, x) * resY * resZ + min(sy, y) * resZ + min(sz, z)] = value;
	bUpdateMesh = true;
}

void MarchingCubes::encodeIsoValsMorton()
{
	for (int x = 0; x < sx1; ++x)
	{
		for (int y = 0; y < sy1; ++y)
		{
			for (int z = 0; z < sz1; ++z)
			{
				int idx = libmorton::morton3D_32_encode(x, y, z);
				isoValsMorton[idx] = isoVals[x * sy1 * sz1 + y * sz1 + z];//getIsoValue(x, y, z);
			}
		}
	}
}

void MarchingCubes::setResolution( int _x, int _y, int _z ){
	
	resX = _x;
	resY = _y;
	resZ = _z;
	resXm1 = resX-1;
	resYm1 = resY-1;
	resZm1 = resZ-1;

	sx = resXm1;
	sy = resYm1;
	sz = resZm1;
	sx1 = resX;
	sy1 = resY;
	sz1 = resZ;
	dx = 1.0 / sx;
	dy = 1.0 / sy;
	dz = 1.0 / sz;
	
	if (isoVals != nullptr) delete[] isoVals;
	isoVals = new float[resX * resY * resZ];
	if (thresCmpArray != nullptr) delete[] thresCmpArray;
	thresCmpArray = new float[resX * resY * resZ];
	if (thresCmpIntArray != nullptr) delete[] thresCmpIntArray;
	thresCmpIntArray = new int[resX * resY * resZ];
	if (thresCmpShortArray != nullptr) delete[] thresCmpShortArray;
	thresCmpShortArray = new short[resX * resY * resZ];
	if (edgeInterpVal != nullptr) delete[] edgeInterpVal;
	edgeInterpVal = new float[resX * resY * resZ * 3];

	// Level-by-level
	thresCmpLevel = new bool[2 * sy1 * sz1];
	thresCmpLevelInt = new int[2 * sy1 * sz1];
	cubeIndexLevel = new int[sy * sz];

	zIndex_EdgeX = new int[3 * sy1 * sz1];
	zIndex_EdgeY = zIndex_EdgeX + sy1 * sz1;
	zIndex_EdgeZ = zIndex_EdgeY + sy1 * sz1;

	yStart_EdgeX = new int[3 * sy1];
	yStart_EdgeY = yStart_EdgeX + sy1;
	yStart_EdgeZ = yStart_EdgeY + sy1;

	edgeIndexArray = new int[sz];

	//isoValsMorton = new float[sx1 * sy1 * sz1];

	//encodeIsoValsMorton();
}

//void MarchingCubes::wipeIsoValues( float value){
//	
//	std::fill(isoVals.begin(), isoVals.end(), value);
//
//}


void MarchingCubes::clear(){
}

void MarchingCubes::exportObj( string fileName ){
	//super simple obj export. doesn;t have any texture coords, materials, or anything super special.
	//just vertices and normals and not optimized either...
	
	//write file
	string title = fileName;
	fileName = /*"../output/" +*/ title + ".obj";
	char *fn = (char*)fileName.c_str();
	
	ofstream outfile (fn);
	
	
	Vector3f v, n;
	
	for( int i=0; i<vertices.size(); ++i){
		v = vertices[i];
		outfile<< "v ";
		outfile<< v.x << " ";
		outfile<< v.y << " ";
		outfile<< v.z << endl;
	}
	outfile << endl;
	
	//for( int i=0; i<vertexCount; i++){
	//	//n = normMat * normals[i];
	//	n = normals[i];
	//	outfile<< "vn ";
	//	outfile<< n.x << " ";
	//	outfile<< n.y << " ";
	//	outfile<< n.z << endl;
	//}
	//outfile << endl;
	
	//this only works for triangulated meshes
	if (indices.empty())
	{
		for (int i = 0; i < vertexCount; i += 3)
		{
			outfile << "f ";
			for (int j = 0; j < 3; j++)
			{
				outfile << i + j + 1;  // Vertex index starts from 1 in obj format.
				//outfile << "/" << i + j + 1;  // This is the normal index.
				outfile << " ";
			}
			outfile << "\n";
		}
	}
	else
	{
		int nFaces = indices.size() / 3;
		for (int i = 0; i < nFaces; ++i)
		{
			outfile << "f ";
			for (int j = 0; j < 3; j++)
			{
				outfile << indices[3*i + j] + 1;
				//outfile << "/" << indices[i + j] + 1;
				outfile << " ";
			}
			outfile << "\n";
		}
	}

	outfile << endl;
	
	outfile.close();
	
	cout <<"model exported as: "<< title <<".obj"<< endl;

}
