#include "MarchingCubes.h"

MarchingCubes::MarchingCubes(){
	threshold = .5;
	bSmoothed = true;
	flipNormalsValue = -1;	
}

MarchingCubes::~MarchingCubes() {
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
		std::fill( gridPointComputed.begin(), gridPointComputed.end(), 0 );
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

	std::fill(gridPointComputed.begin(), gridPointComputed.end(), 0);
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
			for (z = k; z <= z + bZ; z++) {
				thresCmp[idx] = getIsoValue(x, y, z) > threshold;
				idx++;
			}
		}
	}

	idx = 0;
	for (x = i; x < i + bX; x++) {
		for (y = j; y < j + bY; y++) {
			for (z = k; z < k + bZ; z++) {
				int grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
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
	for (x = i; x < i + bX; x++) {
		for (y = j; y < j + bY; y++) {
			for (z = k; z < k + bZ; z++) {
				idx = x * bY * bZ + y * bZ + z;
				grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
				if (edgeTable[cubeIndices[idx]] & 1) vertexInterp(threshold, x, y, z, x+1, y, z, bVertList[grid_idx * 3], dummyN);
				if (edgeTable[cubeIndices[idx]] & 8) vertexInterp(threshold, x, y, z, x, y+1, z, bVertList[grid_idx * 3 + 1], dummyN);
				if (edgeTable[cubeIndices[idx]] & 256) vertexInterp(threshold, x, y, z, x, y, z+1, bVertList[grid_idx * 3 + 2], dummyN);
			}
			// Boundary z
			idx = x * bY * bZ + y * bZ + (z-1);
			grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
			if (edgeTable[cubeIndices[idx]] & 16) vertexInterp(threshold, x, y, z, x+1, y, z, bVertList[grid_idx * 3], dummyN);
			if (edgeTable[cubeIndices[idx]] & 128) vertexInterp(threshold, x, y, z, x, y+1, z, bVertList[grid_idx * 3 + 1], dummyN);
		}
		// Boundary y
		for (z = k; z < k + bZ; z++) {
			idx = x * bY * bZ + (y-1) * bZ + z;
			grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
			if (edgeTable[cubeIndices[idx]] & 4) vertexInterp(threshold, x, y, z, x+1, y, z, bVertList[grid_idx * 3], dummyN);
			if (edgeTable[cubeIndices[idx]] & 2048) vertexInterp(threshold, x, y, z, x, y, z+1, bVertList[grid_idx * 3 + 2], dummyN);
		}
		// Boundary y,z
		idx = x * bY * bZ + (y-1) * bZ + (z-1);
		grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
		if (edgeTable[cubeIndices[idx]] & 64) vertexInterp(threshold, x, y, z, x+1, y, z, bVertList[grid_idx * 3], dummyN);
	}
	// Boundary x
	for (y = j; y < j + bY; y++) {
		for (z = k; z < k + bZ; z++) {
			idx = (x-1) * bY * bZ + y * bZ + z;
			grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
			if (edgeTable[cubeIndices[idx]] & 2) vertexInterp(threshold, x, y, z, x, y+1, z, bVertList[grid_idx * 3 + 1], dummyN);
			if (edgeTable[cubeIndices[idx]] & 512) vertexInterp(threshold, x, y, z, x, y, z+1, bVertList[grid_idx * 3 + 2], dummyN);
		}
		// Boundary x,z
		idx = (x-1) * bY * bZ + y * bZ + (z-1);
		grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
		if (edgeTable[cubeIndices[idx]] & 32) vertexInterp(threshold, x, y, z, x, y+1, z, bVertList[grid_idx * 3 + 1], dummyN);
	}
	// Boundary x,y
	for (z = k; z < k + bZ; z++) {
		idx = (x-1) * bY * bZ + (y-1) * bZ + z;
		grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
		if (edgeTable[cubeIndices[idx]] & 1024) vertexInterp(threshold, x, y, z, x, y, z+1, bVertList[grid_idx * 3 + 2], dummyN);
	}
	
	idx = 0;
	for (x = i; x < i + bX; x++) {
		for (y = j; y < j + bY; y++) {
			for (z = k; z < k + bZ; z++) {
				grid_idx = x * (bY + 1) * (bZ + 1) + y * (bZ + 1) + z;
				for (int ti = 0; triTable[cubeIndices[idx]][ti] != -1; ti += 3) {
					for (int tj = 0; tj < 3; tj++) {
						switch (triTable[cubeIndices[idx]][ti+tj]) {
						case 0: // i,j,k - i1,j,k
							vertices[vertexCount++] = vertList[grid_idx * 3]; break;
						case 1: // i1,j,k - i1,j1,k
							vertices[vertexCount++] = vertList[(grid_idx + (bY+1) * (bZ+1)) * 3 + 1]; break;
						case 2: // i,j1,k - i1,j1,k
							vertices[vertexCount++] = vertList[(grid_idx + (bZ+1)) * 3]; break;
						case 3: // i,j,k - i,j1,k
							vertices[vertexCount++] = vertList[grid_idx * 3 + 1]; break;
						case 4: // i,j,k1 - i1,j,k1
							vertices[vertexCount++] = vertList[(grid_idx + 1) * 3]; break;
						case 5: // i1,j,k1 - i1,j1,k1
							vertices[vertexCount++] = vertList[(grid_idx + (bY+1) * (bZ+1) + 1) * 3 + 1]; break;
						case 6: // i,j1,k1 - i1,j1,k1
							vertices[vertexCount++] = vertList[(grid_idx + (bZ+1) + 1) * 3]; break;
						case 7: // i,j,k1 - i,j1,k1
							vertices[vertexCount++] = vertList[(grid_idx + 1) * 3 + 1]; break;
						case 8: // i,j,k - i,j,k1
							vertices[vertexCount++] = vertList[grid_idx * 3 + 2]; break;
						case 9: // i1,j,k - i1,j,k1
							vertices[vertexCount++] = vertList[(grid_idx + (bY+1) * (bZ+1)) * 3 + 2]; break;
						case 10: // i1,j1,k - i1,j1,k1
							vertices[vertexCount++] = vertList[(grid_idx + (bY+1) * (bZ+1) + (bZ+1)) * 3 + 2]; break;
						case 11: // i,j1,k - i,j1,k1
							vertices[vertexCount++] = vertList[(grid_idx + (bZ+1)) * 3 + 2]; break;
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
