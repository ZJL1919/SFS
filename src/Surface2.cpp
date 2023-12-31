/*
 * Surface2.cpp
 *
 *  Created on: June 7, 2020
 *      Author: Jinglong
 */

#include "Helper0.hpp"
#include "Skeleton2.hpp"




void PrintVector(vector<double>& p){

	for (int i = 0; i < 3; i++){
		cout << p[i] << " ";
	}
	cout << endl;
}


int_type_t LookUpAndPlaceAndReturn(int_type_t self_index, vector<SkelGraph>& SG, ReconstructionStructure* RS, int_type_t* voxel_points_map, vector<double>& points,
		int_type_t x_index, int_type_t y_index, int_type_t z_index, int_type_t sub_x, int_type_t sub_y, int_type_t sub_z, int_type_t& point_counter,
		int_type_t big_number){

	// voxel points map  ....
	// relative to self ...
	// x = [0, 1], y =[0, 1], z = [0, 1]
	// then index is x*4 + y*2 + z
	int_type_t nx, ny, nz;
	int_type_t bx, by, bz; // bin x y z
	int_type_t point_index;
	int_type_t point_index_neighbor;
	int_type_t neighbor_index;
	bool contained;

	//	cout << "subx , y, z " << sub_x << ", " << sub_y << ", " << sub_z << endl;
	//	cout << "x y z " << x_index << ", " << y_index << ", " << z_index << endl;

	bx = (sub_x == x_index + 1); // 0 if sub_x = x_index*2 b/c the sub_x%x_index =  0;
	by = (sub_y == y_index + 1);
	bz = (sub_z == z_index + 1);

	//cout << "bx, by, bz " << bx << ", " << by << ", " << bz << endl;


	point_index = self_index*8 + 4*bx + 2*by + bz;
	// first, is this item already placed in self?
	//	cout << "Should be size of " << SG.size()*8 << endl;;
	//	cout << "Big number " << big_number << endl;
	//	cout << "Point index " << point_index << endl;
	//	cout << "Index " << voxel_points_map[point_index] << endl;
	//char ch; cin >> ch;
	if (voxel_points_map[point_index] < big_number){


	}	else {
		voxel_points_map[point_index] = point_counter;

		points.push_back(RS->BB[0][0] + (double(2*sub_x))*RS->inter_voxel_distance/2.0);
		points.push_back(RS->BB[0][1] + (double(2*sub_y))*RS->inter_voxel_distance/2.0);
		points.push_back(RS->BB[0][2] + (double(2*sub_z))*RS->inter_voxel_distance/2.0);

		// need to go through the neighbors and see if they have this as a point.
		for (int j =0, jn = SG[self_index].neighbors.size(); j < jn; j++){

			neighbor_index = SG[self_index].neighbors[j]->voxel_id;
			RS->ReturnXYZIndicesFromIndex(SG[self_index].neighbors[j]->grid_id, nx, ny, nz);

			contained = ((sub_x == nx) || (sub_x == nx + 1)) && ((sub_y == ny) || (sub_y == ny + 1)) && ((sub_z == nz) || (sub_z == nz + 1));


			if (contained){
				//cout << "Contained ! " << contained << endl;
				bx = (sub_x == nx + 1); // 0 if sub_x = x_index*2 b/c the sub_x%x_index =  0;
				by = (sub_y == ny + 1);
				bz = (sub_z == nz + 1);

				point_index_neighbor = neighbor_index*8 + 4*bx + 2*by + bz;
				voxel_points_map[point_index_neighbor] = point_counter;
			}
		}




		point_counter++;
	}

	//cout << "Completed!" << endl;
	return voxel_points_map[point_index] ;
}

void ReconstructionStructure::GenerateAndWriteSurfaceInPlyFormat(string outdir, int iteration, string prefix, int* cs){

	//cout << "Line 236 " << endl;cin >> ch;

	char ch;
	string filename = outdir + "smoothed_files/" + prefix + ToString<int>(iteration) + ".ply";
	// we'll refer to the points on the mesh by a subsampled version of the
	vector<double> Xhm(3,0);
	vector<int_type_t> sub_number_per_dim(3, 0);

	vector<int_type_t> color(3, 100);
	color[1] = 200;  // was green -- too dark?
	color[2] = 200;
	color[0] = 200;

	if (prefix == "model"){
		color[0] = 255;
		color[1] = 0;
		color[2] = 0;
	}	else {
		color[0] = 0;
		color[1] = 0;
		color[2] = 255;
	}

	//	if (prefix == "initial"){
	//		color[0] = 0;
	//		color[1] = 255;
	//		color[2] = 0;
	//	};

	if (cs != 0){
		color[0] = cs[0];
		color[1] = cs[1];
		color[2] = cs[2];
	}


	for (int i = 0; i < 3; i++){
		sub_number_per_dim[i] = number_voxels_per_dim[i]*2 + 1;
	}

	vector<SkelGraph> SurfaceGraph;
	CreateSurfaceGraphFromGrid(*this, SurfaceGraph);
	//CreateGraphFromGrid(*this, SurfaceGraph);

	cout << "After SG ... " << SurfaceGraph.size() << endl;
	//cin >> ch;

	// now process the points to create unique indices
	int_type_t big_number = sub_number_per_dim[0]* sub_number_per_dim[1]*sub_number_per_dim[2];
	//vector<int_type_t> voxel_points_map;

	int_type_t* voxel_points_map = new int_type_t[SurfaceGraph.size()*8];
	for (int_type_t i = 0, in = SurfaceGraph.size()*8; i < in; i++){
		voxel_points_map[i] = big_number;
	}

	int_type_t point_counter = 0;




	//cout << "Line 256 " << endl;cin >> ch;

	int_type_t above_x, below_x, above_y, below_y, above_z, below_z;
	int_type_t x_index, y_index, z_index;
	bool outer_shell;
	int_type_t voxel_index = 0;
	int_type_t sub_index;
	int_type_t number_points = 0;

	//	vector<int_type_t> start_index_for_face;
	// faces are all 4 sided
	vector<int_type_t> faces;
	vector<double> points;
	vector<int_type_t> current_face(4, 0);
	outer_shell = false;
	int_type_t current_point_index;
	int_type_t sub_x, sub_y, sub_z;

	//cout << "Line 273 " << endl;cin >> ch;
	for (int_type_t i = 0, in = SurfaceGraph.size(); i < in; i++){
		// Setting up
		voxel_index = SurfaceGraph[i].grid_id;
		outer_shell = false;

		ReturnXYZIndicesFromIndex(voxel_index, x_index, y_index, z_index);

		// this done in order .... so if the current index is smaller than all neighbors,
		//		for (int j = 0; j < 8; j++){
		//			voxel_points_map.push_back(big_number);
		//		}


		// TESTING X DIRECTION ..... ////

		if (x_index == 0){
			outer_shell = true;
		}	else {
			below_x = ReturnIndexFromXYZIndices(x_index - 1, y_index, z_index);

			if (configuration_grid[below_x] == true){
				outer_shell = true;
			}
		}

		if (outer_shell == true){
			// below x is a face ....


			// Do I share points on this face with any of the neighbors?
			// come up with a list using the subsample indices.  Look up -- increment and place if so.

			sub_x = x_index;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			///////////////////////////////////

			sub_x = x_index;
			sub_y = y_index;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			///////////////////////////

			sub_x = x_index;
			sub_y = y_index;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			///////////////////////////

			sub_x = x_index;
			sub_y = y_index + 1;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			///////////////////////////



			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}

			//i = in;

		}

		// + x
		outer_shell = false;
		if (x_index == number_voxels_per_dim[0] - 1){
			outer_shell =true;
		}	else {
			above_x = ReturnIndexFromXYZIndices(x_index + 1, y_index, z_index);

			if (configuration_grid[above_x] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}

		}


		/// -y
		//					if (y_index == 0){
		//						outer_shell =true;
		//					}	else {
		//						if (y_index == number_voxels_per_dim[1] - 1){
		//							outer_shell =true;
		//						}	else {
		//							below_y = ReturnIndexFromXYZIndices(x_index, y_index - 1, z_index);
		//
		//							if (configuration_grid[below_y] == true){
		//								outer_shell = true;
		//							}	else {
		//								above_y = ReturnIndexFromXYZIndices(x_index, y_index + 1, z_index);
		//
		//								if (configuration_grid[above_y] == true){
		//									outer_shell = true;
		//								}
		//
		//							}
		//						}
		//					}
		outer_shell = false;
		if (y_index == 0){
			outer_shell =true;
		}	else {
			below_y = ReturnIndexFromXYZIndices(x_index, y_index - 1, z_index);

			if (configuration_grid[below_y] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index;
			sub_y = y_index;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}
		}

		// + y
		outer_shell = false;
		if (y_index == number_voxels_per_dim[1] - 1){
			outer_shell =true;
		}	else {
			above_y = ReturnIndexFromXYZIndices(x_index, y_index + 1, z_index);

			if (configuration_grid[above_y] == true){
				outer_shell = true;
			}
		}

		if (outer_shell == true){
			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}

		}
		// -z
		outer_shell = false;

		if (z_index == 0){
			outer_shell =true;
		}	else {
			below_z = ReturnIndexFromXYZIndices(x_index, y_index, z_index - 1);

			if (configuration_grid[below_z] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){
			sub_x = x_index + 0;
			sub_y = y_index + 0;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}



		}

		// + z
		outer_shell = false;

		if (z_index == number_voxels_per_dim[2] - 1){
			outer_shell =true;
		}	else {
			above_z = ReturnIndexFromXYZIndices(x_index, y_index, z_index + 1);

			if (configuration_grid[above_z] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){
			sub_x = x_index + 0;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}
		}

	}

	cout << "After processing points " << endl;
	cout << "Number points " << point_counter << endl;
	//cin >> ch;

	delete [] voxel_points_map;
	//
	//	//cout << "Line 896 " << endl;cin >> ch;
	//	delete [] dense_to_sparse_map;
	//
	cout << "Before actual write ... " << endl;
	//	char fg; cin >> fg;
	if (points.size() > 0){
		WritePlyFile(filename, sub_number_per_dim, points, faces, color);
	}
	cout << "After actual write" << endl;

	//exit(1);
	//	//	cin >> fg;


}

void ReconstructionStructure::GenerateAndWriteSurfaceInPlyFormat(string outdir, int iteration, string prefix, vector<int>& cs){

	//cout << "Line 236 " << endl;cin >> ch;

	char ch;
	string filename = outdir + "smoothed_files/" + prefix + ToString<int>(iteration) + ".ply";
	// we'll refer to the points on the mesh by a subsampled version of the
	vector<double> Xhm(3,0);
	vector<int_type_t> sub_number_per_dim(3, 0);

	vector<int_type_t> color(3, 100);
	color[1] = 200;  // was green -- too dark?
	color[2] = 200;
	color[0] = 200;

	if (prefix == "model"){
		color[0] = 255;
		color[1] = 0;
		color[2] = 0;
	}	else {
		color[0] = 0;
		color[1] = 0;
		color[2] = 255;
	}

	//	if (prefix == "initial"){
	//		color[0] = 0;
	//		color[1] = 255;
	//		color[2] = 0;
	//	};

	if (cs.size() != 0){
		color[0] = cs[0];
		color[1] = cs[1];
		color[2] = cs[2];
	}


	for (int i = 0; i < 3; i++){
		sub_number_per_dim[i] = number_voxels_per_dim[i]*2 + 1;
	}

	vector<SkelGraph> SurfaceGraph;
	CreateSurfaceGraphFromGrid(*this, SurfaceGraph);
	//CreateGraphFromGrid(*this, SurfaceGraph);

	cout << "After SG ... " << SurfaceGraph.size() << endl;
	//cin >> ch;

	// now process the points to create unique indices
	int_type_t big_number = sub_number_per_dim[0]* sub_number_per_dim[1]*sub_number_per_dim[2];
	//vector<int_type_t> voxel_points_map;

	int_type_t* voxel_points_map = new int_type_t[SurfaceGraph.size()*8];
	for (int_type_t i = 0, in = SurfaceGraph.size()*8; i < in; i++){
		voxel_points_map[i] = big_number;
	}

	int_type_t point_counter = 0;




	//cout << "Line 256 " << endl;cin >> ch;

	int_type_t above_x, below_x, above_y, below_y, above_z, below_z;
	int_type_t x_index, y_index, z_index;
	bool outer_shell;
	int_type_t voxel_index = 0;
	int_type_t sub_index;
	int_type_t number_points = 0;

	//	vector<int_type_t> start_index_for_face;
	// faces are all 4 sided
	vector<int_type_t> faces;
	vector<double> points;
	vector<int_type_t> current_face(4, 0);
	outer_shell = false;
	int_type_t current_point_index;
	int_type_t sub_x, sub_y, sub_z;

	//cout << "Line 273 " << endl;cin >> ch;
	for (int_type_t i = 0, in = SurfaceGraph.size(); i < in; i++){
		// Setting up
		voxel_index = SurfaceGraph[i].grid_id;
		outer_shell = false;

		ReturnXYZIndicesFromIndex(voxel_index, x_index, y_index, z_index);

		// this done in order .... so if the current index is smaller than all neighbors,
		//		for (int j = 0; j < 8; j++){
		//			voxel_points_map.push_back(big_number);
		//		}


		// TESTING X DIRECTION ..... ////

		if (x_index == 0){
			outer_shell = true;
		}	else {
			below_x = ReturnIndexFromXYZIndices(x_index - 1, y_index, z_index);

			if (configuration_grid[below_x] == true){
				outer_shell = true;
			}
		}

		if (outer_shell == true){
			// below x is a face ....


			// Do I share points on this face with any of the neighbors?
			// come up with a list using the subsample indices.  Look up -- increment and place if so.

			sub_x = x_index;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			///////////////////////////////////

			sub_x = x_index;
			sub_y = y_index;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			///////////////////////////

			sub_x = x_index;
			sub_y = y_index;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			///////////////////////////

			sub_x = x_index;
			sub_y = y_index + 1;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			///////////////////////////



			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}

			//i = in;

		}

		// + x
		outer_shell = false;
		if (x_index == number_voxels_per_dim[0] - 1){
			outer_shell =true;
		}	else {
			above_x = ReturnIndexFromXYZIndices(x_index + 1, y_index, z_index);

			if (configuration_grid[above_x] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}

		}


		/// -y
		//					if (y_index == 0){
		//						outer_shell =true;
		//					}	else {
		//						if (y_index == number_voxels_per_dim[1] - 1){
		//							outer_shell =true;
		//						}	else {
		//							below_y = ReturnIndexFromXYZIndices(x_index, y_index - 1, z_index);
		//
		//							if (configuration_grid[below_y] == true){
		//								outer_shell = true;
		//							}	else {
		//								above_y = ReturnIndexFromXYZIndices(x_index, y_index + 1, z_index);
		//
		//								if (configuration_grid[above_y] == true){
		//									outer_shell = true;
		//								}
		//
		//							}
		//						}
		//					}
		outer_shell = false;
		if (y_index == 0){
			outer_shell =true;
		}	else {
			below_y = ReturnIndexFromXYZIndices(x_index, y_index - 1, z_index);

			if (configuration_grid[below_y] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index;
			sub_y = y_index;
			sub_z = z_index;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}
		}

		// + y
		outer_shell = false;
		if (y_index == number_voxels_per_dim[1] - 1){
			outer_shell =true;
		}	else {
			above_y = ReturnIndexFromXYZIndices(x_index, y_index + 1, z_index);

			if (configuration_grid[above_y] == true){
				outer_shell = true;
			}
		}

		if (outer_shell == true){
			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}

		}
		// -z
		outer_shell = false;

		if (z_index == 0){
			outer_shell =true;
		}	else {
			below_z = ReturnIndexFromXYZIndices(x_index, y_index, z_index - 1);

			if (configuration_grid[below_z] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){
			sub_x = x_index + 0;
			sub_y = y_index + 0;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 0;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}



		}

		// + z
		outer_shell = false;

		if (z_index == number_voxels_per_dim[2] - 1){
			outer_shell =true;
		}	else {
			above_z = ReturnIndexFromXYZIndices(x_index, y_index, z_index + 1);

			if (configuration_grid[above_z] == true){
				outer_shell = true;
			}
		}

		if (outer_shell){
			sub_x = x_index + 0;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[0] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 0;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[1] = current_point_index;

			////////////////

			sub_x = x_index + 1;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[2] = current_point_index;

			////////////////

			sub_x = x_index + 0;
			sub_y = y_index + 1;
			sub_z = z_index + 1;
			current_point_index = LookUpAndPlaceAndReturn(i, SurfaceGraph, this, voxel_points_map, points, x_index, y_index, z_index,
					sub_x, sub_y, sub_z, point_counter, big_number);

			current_face[3] = current_point_index;

			////////////////

			for (int f = 0; f < 4; f++){
				faces.push_back(current_face[f]);
			}
		}

	}

	cout << "After processing points " << endl;
	cout << "Number points " << point_counter << endl;
	//cin >> ch;

	delete [] voxel_points_map;
	//
	//	//cout << "Line 896 " << endl;cin >> ch;
	//	delete [] dense_to_sparse_map;
	//
	cout << "Before actual write ... " << endl;
	//	char fg; cin >> fg;
	if (points.size() > 0){
		WritePlyFile(filename, sub_number_per_dim, points, faces, color);
	}
	cout << "After actual write" << endl;

	//exit(1);
	//	//	cin >> fg;


}


void ReconstructionStructure::WritePlyFile(string outfile,
		vector<int_type_t>& subdims,
		vector<double>& points,
		vector<int_type_t>& faces,
		vector<int_type_t>& color){

	// each vertex needs a color ....

	cout << "Writing to " << outfile << endl;
	std::ofstream out;
	out.open(outfile.c_str());

	out << "ply" << endl;
	out << "format ascii 1.0" << endl;
	out << "element vertex " << points.size()/3 << endl;
	out << "property float x" << endl;
	out << "property float y" << endl;
	out << "property float z" << endl;
	out << "property uchar red" << endl;
	out << "property uchar green" << endl;
	out << "property uchar blue" << endl;
	out << "property uchar alpha" << endl;
	out << "element face " << faces.size()/2 << endl;
	out << "property list uchar int vertex_indices"<< endl;
	out << "end_header" << endl;


	for (int_type_t i = 0; i < (points.size()/3); i++){
		out << points[3*i] << " " <<  points[3*i + 1]  << " " <<  points[3*i + 2]  << " ";
		out << color[0] << " " <<  color[1] << " " << color[2] << " 255" << endl;
	}

	for (int_type_t i = 0; i < (faces.size()/4); i++){
		out << "3 " << faces[4*i]  << " " <<  faces[4*i + 1]   << " " <<  faces[4*i + 2]  << endl; //<< " " << faces[4*i + 3] << endl;;
		out << "3 " << faces[4*i]  << " " <<  faces[4*i + 2]   << " " <<  faces[4*i + 3]  << endl;
	}

	out << endl;

	out.close();
}


