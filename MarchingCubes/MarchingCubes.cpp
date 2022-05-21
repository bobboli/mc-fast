#include "MarchingCubes.h"

MarchingCubes::MarchingCubes(){
	threshold = .5;
	bSmoothed = true;
	flipNormalsValue = -1;	
}

MarchingCubes::~MarchingCubes() {
	if (isoValArray != nullptr) delete[] isoValArray;
	if (thresCmpArray != nullptr) delete[] thresCmpArray;
	if (thresCmpIntArray != nullptr) delete[] thresCmpIntArray;
	if (edgeInterpVal != nullptr) delete[] edgeInterpVal;
	if (thresCmp != nullptr) delete[] thresCmp;
	if (cubeIndices != nullptr) delete[] cubeIndices;
	if (bVertList != nullptr) delete[] bVertList;
}

void MarchingCubes::setMaxVertexCount( int _maxVertexCount ){
	maxVertexCount = _maxVertexCount;
	beenWarned = false;
}


void MarchingCubes::setup( int resX, int resY, int resZ, int _maxVertexCount){
	
	clear();

	//up.set(0,1,0);
	
	//float boxVerts[] = {-.5, -.5, -.5, .5, -.5, -.5, -.5, .5, -.5, .5, .5, -.5, -.5, -.5, .5, .5, -.5, .5, -.5, .5, .5, .5, .5, .5, -.5, -.5, .5, -.5, -.5, -.5, -.5, .5, .5, -.5, .5, -.5, .5, -.5, .5, .5, -.5, -.5, .5, .5, .5, .5, .5, -.5, -.5, .5, -.5, -.5, -.5, -.5, -.5, .5, .5, -.5, -.5, .5, .5, .5, -.5, .5, -.5, -.5, .5, .5, .5, .5, -.5, .5,};
	//boundaryVbo.setVertexData( boxVerts, 3, 24, GL_STATIC_DRAW );
	
	setResolution( resX, resY, resZ );
	setMaxVertexCount( _maxVertexCount );

	vertexCount = 0;
	
	vertices.resize( maxVertexCount );
	normals.resize( maxVertexCount );
	
	//vbo.setVertexData( &vertices[0], vertices.size(),GL_DYNAMIC_READ );
	//vbo.setNormalData( &normals[0], normals.size(), GL_DYNAMIC_READ );
}

void MarchingCubes::setBlocking(int blockX, int blockY, int blockZ) {
	bX = blockX; 
	bY = blockY;
	bZ = blockZ;

	if (thresCmp != nullptr) delete[] thresCmp;
	if (cubeIndices != nullptr) delete[] cubeIndices;
	if (bVertList != nullptr) delete[] bVertList;

	thresCmp = new bool[(bX + 1) * (bY + 1) * (bZ + 1)];
	cubeIndices = new short[bX * bY * bZ];
	bVertList = new Vector3f[(bX + 1) * (bY + 1) * (bZ + 1) * 3];
}

void MarchingCubes::update(float _threshold){
	// update every time
	
	//if( bUpdateMesh || threshold != _threshold ){
			
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
		
		//updateTransformMatrix();
		
		//vbo.updateVertexData( &vertices[0], vertexCount );
		//vbo.updateNormalData( &normals[0], vertexCount );
		
	/*	bUpdateMesh = false;
	}*/
}

void MarchingCubes::update_block(float _threshold) {
	threshold = _threshold;

	//std::fill(gridPointComputed.begin(), gridPointComputed.end(), 0);
	vertexCount = 0;

	int x, y, z;
	for (x = 0; x < resX - bX; x += bX) {
		for (y = 0; y < resY - bY; y += bY) {
			for (z = 0; z < resZ - bZ; z += bZ) {
				polygonise_block(x, y, z, bX, bY, bZ);
			}
			for (; z < resZ - 1; z++) {
				polygonise_block(x, y, z, bX, bY, 1);
			}
		}
		for (; y < resY - 1; y++) {
			for (z = 0; z < resZ - bZ; z += bZ) {
				polygonise_block(x, y, z, bX, 1, bZ);
			}
			for (; z < resZ - 1; z++) {
				polygonise_block(x, y, z, bX, 1, 1);
			}
		}
	}

	for (; x < resX - 1; x += bX) {
		for (y = 0; y < resY - bY; y += bY) {
			for (z = 0; z < resZ - bZ; z += bZ) {
				polygonise_block(x, y, z, 1, bY, bZ);
			}
			for (; z < resZ - 1; z++) {
				polygonise_block(x, y, z, 1, bY, 1);
			}
		}
		for (; y < resY - 1; y++) {
			for (z = 0; z < resZ - bZ; z += bZ) {
				polygonise_block(x, y, z, 1, 1, bZ);
			}
			for (; z < resZ - 1; z++) {
				polygonise_block(x, y, z, 1, 1, 1);
			}
		}
	}
}

void MarchingCubes::update_vec(float _threshold) {
	threshold = _threshold;

	vertexCount = 0;

	// if we use 1-byte integer iso value and threshold (as in .vol files), we can do better than this

	// compare stage: build global cmp array
	__m256 vt = _mm256_set1_ps(threshold);
	__m256 c1 = _mm256_set1_ps(1);
	int num = resX * resY * resZ, n;
	int dx = resY * resZ, dy = resZ;
	for (n = 0; n < num - 7; n += 8) {
		__m256 vals = _mm256_loadu_ps(isoValArray+n);
		__m256 cmp = _mm256_cmp_ps(vals, vt, _CMP_GT_OQ);
		__m256 res = _mm256_and_ps(cmp, c1);
		// don't have _mm256_store_epi32, only avaiable in avx512
		_mm256_storeu_epi32(thresCmpIntArray+n, _mm256_cvtps_epi32(res));
	}
	for (; n < num; ++n)
	{
		thresCmpIntArray[n] = isoValArray[n] > threshold ? 1 : 0;
	}

	for (int i = 0; i < num; ++i)
	{
		int val = isoValArray[i] > threshold ? 1 : 0;
		if (val != thresCmpIntArray[i])
			cout << "!!!!!" << endl;
	}
	
	// global intersection
	int x, y, z;
	// check if edge is active for three directions
	for (x = 0; x < resX-1; x++) {
		for (y = 0; y < resY-1; y++) {
			for (z = 0; z < resZ - 8; z += 8) {
				int idx = x * dx + y * dy + z;
				// i j k
				__m256i b = _mm256_loadu_epi32(thresCmpIntArray+idx);

				// i+1 j k
				__m256i b1 = _mm256_loadu_epi32(thresCmpIntArray+idx+dx);
				__m256i cmp1 = _mm256_cmpeq_epi32(b1, b);
				unsigned int mask1 = _mm256_movemask_epi8(cmp1);
				if (mask1 != 0xffffffff)
				{
					// should gather the active edges into vector
					__m256 val_start = _mm256_loadu_ps(isoValArray+idx);
					__m256 val_end = _mm256_loadu_ps(isoValArray+idx+dx);
					__m256 denominator = _mm256_sub_ps(val_end, val_start);
					__m256 numerator = _mm256_sub_ps(vt, val_start);
					__m256 res = _mm256_div_ps(numerator, denominator);
					_mm256_storeu_ps(edgeInterpVal + idx, res);
				}

				// i j+1 k
				__m256i b2 = _mm256_loadu_epi32(thresCmpIntArray + idx + dy);
				__m256i cmp2 = _mm256_cmpeq_epi32(b2, b);
				unsigned int mask2 = _mm256_movemask_epi8(cmp2);
				if (mask2 != 0xffffffff)
				{
					// should gather the active edges into vector
					__m256 val_start = _mm256_loadu_ps(isoValArray + idx);
					__m256 val_end = _mm256_loadu_ps(isoValArray + idx + dy);
					__m256 denominator = _mm256_sub_ps(val_end, val_start);
					__m256 numerator = _mm256_sub_ps(vt, val_start);
					__m256 res = _mm256_div_ps(numerator, denominator);
					_mm256_storeu_ps(edgeInterpVal + num + idx, res);
				}

				// i j k+1
				__m256i b3 = _mm256_loadu_epi32(thresCmpIntArray + idx + 1);
				__m256i cmp3 = _mm256_cmpeq_epi32(b3, b);
				unsigned int mask3 = _mm256_movemask_epi8(cmp3);
				if (mask3 != 0xffffffff)
				{
					// should gather the active edges into vector
					__m256 val_start = _mm256_loadu_ps(isoValArray + idx);
					__m256 val_end = _mm256_loadu_ps(isoValArray + idx + 1);
					__m256 denominator = _mm256_sub_ps(val_end, val_start);
					__m256 numerator = _mm256_sub_ps(vt, val_start);
					__m256 res = _mm256_div_ps(numerator, denominator);
					_mm256_storeu_ps(edgeInterpVal + num + num + idx, res);
				}
			}
			for (; z < resZ - 1; ++z)
			{
				int idx = x * dx + y * dy + z;
				float start = isoValArray[idx];
				float end_x = isoValArray[idx+dx];
				float end_y = isoValArray[idx+dy];
				float end_z = isoValArray[idx+1];
				float numerator = threshold - start;
				edgeInterpVal[idx] = numerator / (end_x - start);
				edgeInterpVal[idx+num] = numerator / (end_y - start);
				edgeInterpVal[idx+num*2] = numerator / (end_z - start);
			}
		}
	}
	
	for (x = 0; x < resX - 1; x += 1) {
		for (y = 0; y < resY - 1; y += 1) {
			for (z = 0; z < resZ - bZ; z += bZ) {
				polygonise_vec(x, y, z, 1, 1, bZ);
			}
			for (; z < resZ - 1; z++) {
				polygonise(x, y, z);
			}
		}
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
	
	if( vertexCount+3 < maxVertexCount ){
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
			vertices[vertexCount++] = vertList[triTable[cubeindex][i]];
			vertices[vertexCount++] = vertList[triTable[cubeindex][i+1]];
			vertices[vertexCount++] = vertList[triTable[cubeindex][i+2]];
		}
	}
	else if(!beenWarned){
		std::cerr << "ofxMarhingCubes: maximum vertex("+to_string(maxVertexCount)+") count exceded. try increasing the maxVertexCount with setMaxVertexCount()";
		beenWarned = true;
	}
}

void MarchingCubes::polygonise_block(int i, int j, int k, int bX, int bY, int bZ) {
	if (vertexCount >= maxVertexCount) return;

	bUpdateMesh = true;

	Vector3f dummyN;
	int idx, x, y, z;

	idx = 0;
	for (x = i; x <= i + bX; x++) {
		for (y = j; y <= j + bY; y++) {
			for (z = k; z <= k + bZ; z++) {
				thresCmp[idx] = getIsoValue(x, y, z) > threshold;
				idx++;
			}
		}
	}

	idx = 0;
	for (x = 0; x < bX; x++) {
		for (y = 0; y < bY; y++) {
			for (z = 0; z < bZ; z++) {
				int grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
				cubeIndices[idx] = 0;
				cubeIndices[idx] |= thresCmp[grid_idx] ? 1 : 0;
				cubeIndices[idx] |= thresCmp[grid_idx + (bY+1) * (bZ+1)] ? 2 : 0;
				cubeIndices[idx] |= thresCmp[grid_idx + (bY+1) * (bZ+1) + (bZ+1)] ? 4 : 0;
				cubeIndices[idx] |= thresCmp[grid_idx + (bZ+1)] ? 8 : 0;
				cubeIndices[idx] |= thresCmp[grid_idx + 1] ? 16 : 0;
				cubeIndices[idx] |= thresCmp[grid_idx + (bY+1) * (bZ+1) + 1] ? 32 : 0;
				cubeIndices[idx] |= thresCmp[grid_idx + (bY+1) * (bZ+1) + (bZ+1) + 1] ? 64 : 0;
				cubeIndices[idx] |= thresCmp[grid_idx + (bZ+1) + 1] ? 128 : 0;
				idx++;
			}
		}
	}
	
	int grid_idx;
	for (x = 0; x < bX; x++) {
		for (y = 0; y < bY; y++) {
			for (z = 0; z < bZ; z++) {
				idx = x * bY * bZ + y * bZ + z;
				grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
				if (edgeTable[cubeIndices[idx]] & 1) vertexInterp(threshold, x+i, y+j, z+k, x+i+1, y+j, z+k, bVertList[grid_idx * 3], dummyN);
				if (edgeTable[cubeIndices[idx]] & 8) vertexInterp(threshold, x+i, y+j, z+k, x+i, y+j+1, z+k, bVertList[grid_idx * 3 + 1], dummyN);
				if (edgeTable[cubeIndices[idx]] & 256) vertexInterp(threshold, x+i, y+j, z+k, x+i, y+j, z+k+1, bVertList[grid_idx * 3 + 2], dummyN);
			}
			// Boundary z
			idx = x * bY * bZ + y * bZ + (z-1);
			grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
			if (edgeTable[cubeIndices[idx]] & 16) vertexInterp(threshold, x+i, y+j, z+k, x+i+1, y+j, z+k, bVertList[grid_idx * 3], dummyN);
			if (edgeTable[cubeIndices[idx]] & 128) vertexInterp(threshold, x+i, y+j, z+k, x+i, y+j+1, z+k, bVertList[grid_idx * 3 + 1], dummyN);
		}
		// Boundary y
		for (z = 0; z < bZ; z++) {
			idx = x * bY * bZ + (y-1) * bZ + z;
			grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
			if (edgeTable[cubeIndices[idx]] & 4) vertexInterp(threshold, x+i, y+j, z+k, x+i+1, y+j, z+k, bVertList[grid_idx * 3], dummyN);
			if (edgeTable[cubeIndices[idx]] & 2048) vertexInterp(threshold, x+i, y+j, z+k, x+i, y+j, z+k+1, bVertList[grid_idx * 3 + 2], dummyN);
		}
		// Boundary y,z
		idx = x * bY * bZ + (y-1) * bZ + (z-1);
		grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
		if (edgeTable[cubeIndices[idx]] & 64) vertexInterp(threshold, x+i, y+j, z+k, x+i+1, y+j, z+k, bVertList[grid_idx * 3], dummyN);
	}
	// Boundary x
	for (y = 0; y < bY; y++) {
		for (z = 0; z < bZ; z++) {
			idx = (x-1) * bY * bZ + y * bZ + z;
			grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
			if (edgeTable[cubeIndices[idx]] & 2) vertexInterp(threshold, x+i, y+j, z+k, x+i, y+j+1, z+k, bVertList[grid_idx * 3 + 1], dummyN);
			if (edgeTable[cubeIndices[idx]] & 512) vertexInterp(threshold, x+i, y+j, z+k, x+i, y+j, z+k+1, bVertList[grid_idx * 3 + 2], dummyN);
		}
		// Boundary x,z
		idx = (x-1) * bY * bZ + y * bZ + (z-1);
		grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
		if (edgeTable[cubeIndices[idx]] & 32) vertexInterp(threshold, x+i, y+j, z+k, x+i, y+j+1, z+k, bVertList[grid_idx * 3 + 1], dummyN);
	}
	// Boundary x,y
	for (z = 0; z < bZ; z++) {
		idx = (x-1) * bY * bZ + (y-1) * bZ + z;
		grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
		if (edgeTable[cubeIndices[idx]] & 1024) vertexInterp(threshold, x+i, y+j, z+k, x+i, y+j, z+k+1, bVertList[grid_idx * 3 + 2], dummyN);
	}
	
	idx = 0;
	for (x = 0; x < bX; x++) {
		for (y = 0; y < bY; y++) {
			for (z = 0; z < bZ; z++) {
				grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
				for (int ti = 0; triTable[cubeIndices[idx]][ti] != -1; ti += 3) {
					for (int tj = 0; tj < 3; tj++) {
						switch (triTable[cubeIndices[idx]][ti+tj]) {
						case 0: // i,j,k - i1,j,k
							vertices[vertexCount++] = bVertList[grid_idx * 3]; break;
						case 1: // i1,j,k - i1,j1,k
							vertices[vertexCount++] = bVertList[(grid_idx + (bY+1) * (bZ+1)) * 3 + 1]; break;
						case 2: // i,j1,k - i1,j1,k
							vertices[vertexCount++] = bVertList[(grid_idx + (bZ+1)) * 3]; break;
						case 3: // i,j,k - i,j1,k
							vertices[vertexCount++] = bVertList[grid_idx * 3 + 1]; break;
						case 4: // i,j,k1 - i1,j,k1
							vertices[vertexCount++] = bVertList[(grid_idx + 1) * 3]; break;
						case 5: // i1,j,k1 - i1,j1,k1
							vertices[vertexCount++] = bVertList[(grid_idx + (bY+1) * (bZ+1) + 1) * 3 + 1]; break;
						case 6: // i,j1,k1 - i1,j1,k1
							vertices[vertexCount++] = bVertList[(grid_idx + (bZ+1) + 1) * 3]; break;
						case 7: // i,j,k1 - i,j1,k1
							vertices[vertexCount++] = bVertList[(grid_idx + 1) * 3 + 1]; break;
						case 8: // i,j,k - i,j,k1
							vertices[vertexCount++] = bVertList[grid_idx * 3 + 2]; break;
						case 9: // i1,j,k - i1,j,k1
							vertices[vertexCount++] = bVertList[(grid_idx + (bY+1) * (bZ+1)) * 3 + 2]; break;
						case 10: // i1,j1,k - i1,j1,k1
							vertices[vertexCount++] = bVertList[(grid_idx + (bY+1) * (bZ+1) + (bZ+1)) * 3 + 2]; break;
						case 11: // i,j1,k - i,j1,k1
							vertices[vertexCount++] = bVertList[(grid_idx + (bZ+1)) * 3 + 2]; break;
						}
					}
					if (vertexCount >= maxVertexCount && !beenWarned) {
						std::cerr << "ofxMarhingCubes: maximum vertex("+to_string(maxVertexCount)+") count exceded. try increasing the maxVertexCount with setMaxVertexCount()";
						beenWarned = true;
						return;
					}
				}
				idx++;
			}
		}
	}
}

void MarchingCubes::polygonise_vec(int i, int j, int k, int bX, int bY, int bZ) {
	if (vertexCount >= maxVertexCount) return;

	Vector3f dummyN;
	int idx, x, y, z;

	int dx = resY * resZ, dy = resZ;
	i = min(i, resXm1);
	j = min(j, resYm1);
	k = min(k, resZm1);
	int base = i * dx + j * dy + k;
	idx = 0;

	for (z = 0; z < bZ; z++) {
		int grid_idx = z + base;
		cubeIndices[idx] = 0;
		cubeIndices[idx] |= thresCmpIntArray[grid_idx] ? 1 : 0;
		cubeIndices[idx] |= thresCmpIntArray[grid_idx + dx] ? 2 : 0;
		cubeIndices[idx] |= thresCmpIntArray[grid_idx + dx + dy] ? 4 : 0;
		cubeIndices[idx] |= thresCmpIntArray[grid_idx + dy] ? 8 : 0;
		cubeIndices[idx] |= thresCmpIntArray[grid_idx + 1] ? 16 : 0;
		cubeIndices[idx] |= thresCmpIntArray[grid_idx + dx + 1] ? 32 : 0;
		cubeIndices[idx] |= thresCmpIntArray[grid_idx + dx + dy + 1] ? 64 : 0;
		cubeIndices[idx] |= thresCmpIntArray[grid_idx + dy + 1] ? 128 : 0;
		idx++;
	}
	int i1 = i + 1, j1= j + 1; 
	for (int n = 0; n < bZ; ++n)
	{
		int cubeindex = cubeIndices[n];
		if (edgeTable[cubeindex] == 0) continue;

		/* Find the vertices where the surface intersects the cube */
		int kn = k + n;
		int kn1 = k + n + 1;
		if (edgeTable[cubeindex] & 1)		vertexInterp_vec(threshold, i, j, kn, i1, j, kn, vertList[0], normList[0]);
		if (edgeTable[cubeindex] & 2)		vertexInterp_vec(threshold, i1, j, kn, i1, j1, kn, vertList[1], normList[1]);
		if (edgeTable[cubeindex] & 4)		vertexInterp_vec(threshold, i1, j1, kn, i, j1, kn, vertList[2], normList[2]);
		if (edgeTable[cubeindex] & 8)		vertexInterp_vec(threshold, i, j1, kn, i, j, kn, vertList[3], normList[3]);
		if (edgeTable[cubeindex] & 16)		vertexInterp_vec(threshold, i, j, kn1, i1, j, kn1, vertList[4], normList[4]);
		if (edgeTable[cubeindex] & 32)		vertexInterp_vec(threshold, i1, j, kn1, i1, j1, kn1, vertList[5], normList[5]);
		if (edgeTable[cubeindex] & 64)		vertexInterp_vec(threshold, i1, j1, kn1, i, j1, kn1, vertList[6], normList[6]);
		if (edgeTable[cubeindex] & 128)		vertexInterp_vec(threshold, i, j1, kn1, i, j, kn1, vertList[7], normList[7]);
		if (edgeTable[cubeindex] & 256)		vertexInterp_vec(threshold, i, j, kn, i, j, kn1, vertList[8], normList[8]);
		if (edgeTable[cubeindex] & 512)		vertexInterp_vec(threshold, i1, j, kn, i1, j, kn1, vertList[9], normList[9]);
		if (edgeTable[cubeindex] & 1024)	vertexInterp_vec(threshold, i1, j1, kn, i1, j1, kn1, vertList[10], normList[10]);
		if (edgeTable[cubeindex] & 2048)	vertexInterp_vec(threshold, i, j1, kn, i, j1, kn1, vertList[11], normList[11]);

		for (int ii = 0; triTable[cubeindex][ii] != -1; ii += 3) {
			vertices[vertexCount++] = vertList[triTable[cubeindex][ii]];
			vertices[vertexCount++] = vertList[triTable[cubeindex][ii + 1]];
			vertices[vertexCount++] = vertList[triTable[cubeindex][ii + 2]];
		}
	}
}


void MarchingCubes::polygonise_count_ops(int i, int j, int k, operation_counts& counts) {

	if (vertexCount + 3 < maxVertexCount) {
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

		for (i = 0; triTable[cubeindex][i] != -1; i += 3) {
			// Vector3f a = vertList[triTable[cubeindex][i + 1]] - vertList[triTable[cubeindex][i]];
			// Vector3f b = vertList[triTable[cubeindex][i + 2]] - vertList[triTable[cubeindex][i + 1]];

			//counts.fl_add += 6;

			/*vertices[vertexCount] = vertList[triTable[cubeindex][i]];
			vertices[vertexCount + 1] = vertList[triTable[cubeindex][i + 1]];
			vertices[vertexCount + 2] = vertList[triTable[cubeindex][i + 2]];*/

			vertexCount += 3;
		}
	}
}

void MarchingCubes::vertexInterp(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, Vector3f& v, Vector3f& n){
	
	Vector3f& p1 = getGridPoint(i1,j1,k1);
	Vector3f& p2 = getGridPoint(i2,j2,k2);
	
	float& iso1 = getIsoValue(i1,j1,k1);
	float& iso2 = getIsoValue(i2,j2,k2);
	
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

void MarchingCubes::vertexInterp_vec(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, Vector3f& v, Vector3f& n) {

	Vector3f& p1 = getGridPoint(i1, j1, k1);
	Vector3f& p2 = getGridPoint(i2, j2, k2);

	int idx1 = i1 * resZ * resY + j1 * resZ + k1;
	int idx2 = i2 * resZ * resY + j2 * resZ + k2;
	float& iso1 = isoValArray[idx1];
	float& iso2 = isoValArray[idx2];

	
	/*if (abs(threshold - iso1) < 0.00001) {
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
	}*/
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
		v = p1 + (p2 - p1) * (1.0 - edgeInterpVal[idx2]);
		return;
	}
	else if (j2 == j1 + 1)
	{
		v = p1 + (p2 - p1) * edgeInterpVal[idx1 + num];
		return;
	}
	else if (j2 == j1 - 1) 
	{
		v = p1 + (p2 - p1) * (1.0 - edgeInterpVal[idx2 + num]);
		return;
	}
	else if(k2 == k1 + 1)
	{
		v = p1 + (p2 - p1) * edgeInterpVal[idx1 + num * 2];
		return;
	}
	else if (k2 == k1 - 1)
	{
		v = p1 + (p2 - p1) * (1.0 - edgeInterpVal[idx2 + num * 2]);
		return;
	}
}

void MarchingCubes::vertexInterp_count_ops(float threshold, int i1, int j1, int k1, int i2, int j2, int k2, operation_counts& counts) {

	Vector3f& p1 = getGridPoint(i1, j1, k1);
	Vector3f& p2 = getGridPoint(i2, j2, k2);

	float& iso1 = getIsoValue(i1, j1, k1);
	float& iso2 = getIsoValue(i2, j2, k2);

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
	

	if(getGridPointComputed(i,j,k) == 0){
		Vector3f& n = getNormalVal(i, j, k);// normalVals[i][j][k];
		n.set(getIsoValue(min(resXm1, i+1), j, k) - getIsoValue(max(0,i-1),j,k),
			  getIsoValue(i,min(resYm1, j+1),k) - getIsoValue(i,max(0,j-1),k),
			  getIsoValue(i,j,min(resZm1, k+1)) - getIsoValue(i,j,max(0,k-1)));
		
		n.normalize();
		n = n*flipNormalsValue;
		getGridPointComputed(i,j,k) = 1;
	}
};

//void ofxMarchingCubes::drawArrays( vector<Vector3f>* _vertices, vector<Vector3f>* _normals){
//	
//	glEnableClientState(GL_VERTEX_ARRAY);
//	glVertexPointer(3, GL_FLOAT, sizeof((*_vertices)[0]), &(*_vertices)[0].x);
//	
//	if(_normals != NULL){
//		glEnableClientState(GL_NORMAL_ARRAY);
//		glNormalPointer( GL_FLOAT, sizeof((*_normals)[0]), &(*_normals)[0].x);
//	}
//	
//	glDrawArrays(GL_TRIANGLES, 0, vertexCount);
//	
//	glDisableClientState(GL_VERTEX_ARRAY);
//	if(_normals != NULL)	glDisableClientState(GL_NORMAL_ARRAY);
//	
//}

void MarchingCubes::setGridPoints( float _x, float _y, float _z){
	cellDim = Vector3f( 1.f/resXm1, 1.f/resYm1, 1.f/resZm1 );
	
	for(int i=0; i<resX; i++){
		for(int j=0; j<resY; j++){
			for(int k=0; k<resZ; k++){
				getGridPoint( i, j, k ).set(float(i)*cellDim.x-.5f,
											float(j)*cellDim.y-.5f,
											float(k)*cellDim.z-.5f);
			}
		}
	}

	//gridPointsVbo.setVertexData( gridPoints[0].getPtr(), 3, gridPoints.size(), GL_STATIC_DRAW );
}



//void ofxMarchingCubes::drawWireframe(){
//	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//	draw();
//	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//}
//
//void ofxMarchingCubes::drawGrid( bool drawGridPoints){
//	
//	ofPushMatrix();
//	ofMultMatrix( transform.getPtr() );
//	
//	if(drawGridPoints){
//		ofSetColor(128);
//		glPointSize(1);
//		gridPointsVbo.draw(GL_POINTS, 0, gridPointsVbo.getNumVertices());
//	}
//	
//	ofColor(255);
//	ofSetLineWidth(1.5);
//	boundaryVbo.draw(GL_LINES, 0, boundaryVbo.getNumVertices());
//	
//	ofPopMatrix();
//}

void MarchingCubes::setIsoValue( int x, int y, int z, float value){
	getIsoValue(min(resXm1,x), min(resYm1,y), min(resZm1,z)) = value;
	isoValArray[min(resXm1, x) * resY * resZ + min(resYm1, y) * resZ + min(resZm1, z)] = value;
	getGridPointComputed(x,y,z) = 0;
	bUpdateMesh = true;
}

void MarchingCubes::setResolution( int _x, int _y, int _z ){
	
	resX = _x;
	resY = _y;
	resZ = _z;
	resXm1 = resX-1;
	resYm1 = resY-1;
	resZm1 = resZ-1;
	
	isoVals.resize( resX*resY*resZ );
	gridPoints.resize( resX*resY*resZ );
	normalVals.resize( resX*resY*resZ );
	gridPointComputed.resize( resX*resY*resZ );
	
	if (isoValArray != nullptr) delete[] isoValArray;
	isoValArray = new float[resX * resY * resZ];
	if (thresCmpArray != nullptr) delete[] thresCmpArray;
	thresCmpArray = new float[resX * resY * resZ];
	if (thresCmpIntArray != nullptr) delete[] thresCmpIntArray;
	thresCmpIntArray = new int[resX * resY * resZ];
	if (edgeInterpVal != nullptr) delete[] edgeInterpVal;
	edgeInterpVal = new float[resX * resY * resZ * 3];

	setGridPoints( resX*10, resY*10, resZ*10 );
}

void MarchingCubes::wipeIsoValues( float value){
	
	std::fill(gridPointComputed.begin(), gridPointComputed.end(), 0);
	std::fill(isoVals.begin(), isoVals.end(), value);

}


void MarchingCubes::clear(){
	isoVals.clear();
	gridPoints.clear();
	normalVals.clear();
}

void MarchingCubes::exportObj( string fileName ){
	//super simple obj export. doesn;t have any texture coords, materials, or anything super special.
	//just vertices and normals and not optimized either...
	
	//write file
	string title = fileName;
	fileName = /*"../output/" +*/ title + ".obj";
	char *fn = (char*)fileName.c_str();
	
	ofstream outfile (fn);
	
	outfile << "# This file uses centimeters as units for non-parametric coordinates."<< endl;
	
	
	//updateTransformMatrix();
	//ofMatrix4x4 normMat;
	//normMat.makeOrthoNormalOf( transform );
	Vector3f v, n;
	
	for( int i=0; i<vertexCount; i++){
		//v = transform * vertices[i];
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
	for( int i=0; i<vertexCount; i+=3){
		outfile<< "f ";
		for(int j=0; j<3; j++){
			outfile<< i+j+1 ;
			outfile<<"/"<<i+j+1 ;
			outfile<< " ";
		}
		outfile<< endl;
	}
	outfile << endl;
	
	outfile.close();
	
	cout <<"model exported as: "<< title <<".obj"<< endl;
}
