/*
 * Skeleton2.cpp
 *
 *  Created on: Oct 9, 2020
 *      Author: Jinglong
 */

#include "Skeleton2.hpp"


SkelGraph::SkelGraph(){
	voxel_id = 0;
	path_id = 0;
	grid_id = 0;
	cc_id = 0;
}


SkelGraph::SkelGraph(int_type_t vi, int_type_t gi){
	voxel_id = vi;
	grid_id = gi;
	path_id = 0;
	cc_id = 0;
}

SkelGraph::~SkelGraph(){}

void SkelGraph::Print(){
	cout << "Voxel id " << voxel_id << endl;
	cout << "Grid id " << grid_id << endl;
	cout << "Path id " << path_id << endl;
	cout << "Number neighbors " << neighbors.size();
}

void CreateSurfaceGraphFromGrid(ReconstructionStructure& RS, vector<SkelGraph>& SG){

	int_type_t arr[26];
	int_type_t distance[26];
	int_type_t number_ns;
	int_type_t n_index;
	int_type_t d;
	int_type_t j;
	int_type_t number_ns_2;


	int_type_t* dense_to_sparse = new int_type_t[RS.number_voxels_grid]; // will be dealloc'd after this function.
	int_type_t sparse_counts = 0;

	for (int_type_t self_index = 0; self_index < RS.number_voxels_grid; self_index++){


		if (RS.configuration_grid[self_index] == 0){
			number_ns = Return26ConnectedNeighbors(RS, self_index, &arr[0], &distance[0]);

			number_ns_2 = 0;

			for (j = 0; j < number_ns; j++){
				n_index = arr[j];
				if (RS.configuration_grid[n_index] == false){
					number_ns_2++;
				}
			}

			if (number_ns_2 < 26)
			{
				//cout << "Number " << number_ns << endl;
				dense_to_sparse[self_index] = sparse_counts;
				SG.push_back(SkelGraph(sparse_counts, self_index));
				sparse_counts++;
			}
			else {
				dense_to_sparse[self_index] = RS.number_voxels_grid;
			}
		}	else {
			dense_to_sparse[self_index] = RS.number_voxels_grid;
		}
	}



	for (int_type_t voxel_id = 0; voxel_id < sparse_counts; voxel_id++){

		number_ns = Return26ConnectedNeighbors(RS, SG[voxel_id].grid_id, &arr[0], &distance[0]);

		for (j = 0; j < number_ns; j++){
			n_index = arr[j];
			d = distance[j];

			if (dense_to_sparse[n_index] != RS.number_voxels_grid){
				// link in
				SG[voxel_id].neighbors.push_back(&SG.at(dense_to_sparse[n_index]));
				SG[voxel_id].neighbor_distance.push_back(d);
			}	else {
				// not a surface node .... ignore
			}
		}

	}

	cout << "Number of connected nodes .... " << SG.size() << endl;
	delete [] dense_to_sparse;

}

int Return26ConnectedNeighbors(ReconstructionStructure& RS, int_type_t start_voxel, int_type_t* member_array, int_type_t* distance){

	int number_neighbors = 0;
	int_type_t nindex;

	int_type_t x_index, y_index, z_index;

	RS.ReturnXYZIndicesFromIndex(start_voxel, x_index, y_index, z_index);


	int_type_t minx, maxx, miny, maxy, minz, maxz;

	if (x_index == 0){
		minx = 0;
	}	else {
		minx = x_index - 1;
	}

	if (y_index == 0){
		miny = 0;
	}	else {
		miny = y_index - 1;
	}

	if (z_index == 0){
		minz = 0;
	}	else {
		minz = z_index - 1;
	}

	if (z_index == RS.number_voxels_per_dim[2] - 1){
		maxz = z_index;
	}	else {
		maxz = z_index + 1;
	}

	if (y_index == RS.number_voxels_per_dim[1] - 1){
		maxy = y_index;
	}	else {
		maxy = y_index + 1;
	}

	if (x_index == RS.number_voxels_per_dim[0] - 1){
		maxx = x_index;
	}	else {
		maxx = x_index + 1;
	}

	for (int_type_t x0 = minx; x0 <= maxx; x0++){
		for (int_type_t y0 = miny; y0 <= maxy; y0++){
			for (int_type_t z0 = minz; z0 <= maxz; z0++){
				nindex = RS.ReturnIndexFromXYZIndices(x0, y0, z0);
				if (nindex != start_voxel){
					member_array[number_neighbors] = nindex;
					// b/c these are only at most separated by 1 ...this will give us squared distance
					distance[number_neighbors] = (x0 != x_index) + (y0 != y_index) + (z0 != z_index);
					number_neighbors++;
				}
			}
		}
	}


	return number_neighbors;


}




