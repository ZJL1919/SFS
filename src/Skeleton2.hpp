/*
 * Skeleton2.hpp
 *
 *  Created on: Oct 9, 2020
 *      Author: atabb
 */

#ifndef SKELETON2_HPP_
#define SKELETON2_HPP_


#include "ImageData.hpp"
#include "Includes.hpp"
#include "ReconstructionStructure.hpp"
#include "Helper0.hpp"

//#include <opencv/cv.h>
#include <opencv/cxcore.h>
#include <opencv/highgui.h>

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include <list>
#include <ctime>
#include <time.h>

#include <fstream>
#include <sstream>
#include <sys/stat.h>
//#include <sys/time.h>
#include <inttypes.h>

class SkelGraph{
public:
	int_type_t voxel_id; // sparse numbering
	int_type_t path_id; // for skeleton and path
	int_type_t grid_id; // in the voxel grid
	int_type_t cc_id;
	vector<SkelGraph*> neighbors;
	vector<int8_t> neighbor_distance; // this can be much smaller ....

	SkelGraph();

	SkelGraph(int_type_t vi, int_type_t gi);

	~SkelGraph();

	void Print();
};


void CreateSurfaceGraphFromGrid(ReconstructionStructure& RS, vector<SkelGraph>& SG);

int Return26ConnectedNeighbors(ReconstructionStructure& RS, int_type_t start_voxel, int_type_t* member_array, int_type_t* distance);

#endif /* SKELETON2_HPP_ */
