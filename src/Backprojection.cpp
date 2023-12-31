/*
 * Backprojection.cpp
 *
 *  Created on: May 1, 2020
 *      Author: Jinglong
 */
#include "ReconstructionStructure.hpp"
//#include "GraphClass.hpp"
#include "Helper0.hpp"
#include <omp.h>





void ReconstructionStructure::ComputeProjectionInformationSpheresArrRepresentationParallel(){

	char ch;

	// allocate all of the vectors ahead of time.
	vector< vector<double> > Xhm(number_divisions, vector<double>(4, 1));
	vector< vector<double> > x(number_divisions, vector<double>(3, 1));
	vector< vector<double> > Xhmp(number_divisions, vector<double>(4, 1));
	vector< vector<double> > xp(number_divisions, vector<double>(3, 1));

	// make these private
	int_type_t thread_ID;
	int r_index, c_index, r_max, c_max;
	double row_squared_distance_from_center;
	int_type_t last_term_position = 0;
	int_type_t number_terms = term_increment;
	int_type_t voxel_index = 0;
	double r_squared, r_squared_test, r;
	int_type_t interacting_p, interacting_p_within_image;
	bool in_image, original_in_image;
	int_type_t* arr_ptr;
	int_type_t counter;

	int_type_t arr0, arr1;
	bool start_from_nothing = start_indices_terms_per_div.size() == 0;
	// global to all threads
	int_type_t max_per_thread  = number_voxels_grid/number_divisions;

	if (start_from_nothing){
		start_indices_terms_per_div.resize(number_divisions, 0);
		terms_arr_per_div.resize(number_divisions, vector<int_type_t*>());

#pragma omp parallel for
		for (int i = 0; i < number_divisions; i++){
#pragma omp critical
			{
			start_indices_terms_per_div[i] = new int_type_t[max_per_thread + 2]; // to take care fo the ends ...
			terms_arr_per_div[i].push_back(new int_type_t[term_increment]);
			}
		}

	}


	cout << "Line 852 " << endl;

	//vector<int_type_t*> arr_ptrs(number_divisions, 0);


	vector< vector<double> > V(number_cameras, vector<double>(3, 0));
	for (int_type_t inner_cam = 0;
			inner_cam < number_cameras; inner_cam++){
		if (plane_equations[inner_cam][0] != 0 || plane_equations[inner_cam][1] != 0){
			V[inner_cam][1] = plane_equations[inner_cam][0];
			V[inner_cam][0] = -plane_equations[inner_cam][1];
			V[inner_cam][2] = 0;
		}	else {
			V[inner_cam][0] = 0; //plane_equations[inner_cam][0];
			V[inner_cam][1] = -plane_equations[inner_cam][2];
			V[inner_cam][2] = plane_equations[inner_cam][1];
		}//ʲô  ˼    

		for (int r = 0; r < 3; r++){
			for (int c = 0; c < 4; c++){
				cout << projection_matrices[inner_cam][r][c] << " ";
			}
			cout << endl;
		}

		cout << endl;
		NormalizeVector(V[inner_cam]);
		//MultiplyVectorByScalar(V, inter_voxel_distance/2.0, V, 3);
		MultiplyVectorByScalar(V[inner_cam], inter_voxel_distance/2.0, V[inner_cam], 3);
		
		cout << "the " << inner_cam << "V is" << V[inner_cam][0]<<"  " << V[inner_cam][1]<<"  " << V[inner_cam][2] << endl;//0th:ȫ  0  1-18th  0 -4 0  19-35th  0 4 0
	}

	//char fg; cin >> fg;

	//	//int some_found
	//
	//	// change the indices so that start index for terms incorporates the +1 so we don't have to test all of the time.
	//
	// this can be in parallel if we allocate terms ahead of time ... and vey nice b/c then each has its own slot.  Do this later for speed.
	//cout << "After preliminary alloc ... " << endl; cin >> ch;
	omp_set_num_threads(number_divisions);


	#pragma omp parallel private(arr0, arr1, arr_ptr, thread_ID, counter, r_index, c_index, r_max, c_max, row_squared_distance_from_center, last_term_position, number_terms, voxel_index, r_squared, r_squared_test, r, interacting_p, interacting_p_within_image, in_image, original_in_image)
	//for(thread_ID = 0; thread_ID < number_divisions; thread_ID++)
	{
		// preliminaries
		thread_ID = omp_get_thread_num();
		last_term_position = 0;
		number_terms = term_increment;//pow(2,24)
		voxel_index = 0;
		//		arr_ptr = terms_arr_per_div[thread_ID][0];
		counter = 0;
		arr0 = 0;
		arr1 = 0;

		//	for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
		{

			//			last_term_position = 0;
			//			number_terms = term_increment;
			//			voxel_index = 0;
			//			arr_ptr = terms_arr_per_div[thread_ID][0];
			//			counter = 0;
			//			arr0 = 0;
			//			arr1 = 0;

			// these probably made this fail ... r_index switching valuesall of the time ...
			for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

#pragma omp critical
				{
					if (x_index % 20 == 0){
						cout << "Computing projection information " << x_index << " out of " << number_voxels_per_dim[0] << " and thread id " << thread_ID << endl;
					}
				}
				for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
					for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){

						if (voxel_index % number_divisions == thread_ID){
							// can do a ptr here ...
							start_indices_terms_per_div[thread_ID][counter] = last_term_position;
							counter++;//ָ     ؿ   

							// test here is if we place no restrictions on which ones get b-projected, do we get rid of the problem?
							//    set_by_first_d.resize(number_voxels_grid, false)    ʼ  ȫFalse
							//    multi_dim_voxels_to_process[i][j] = false;
							if (set_by_first_d[voxel_index] == false || multi_dim_voxels_to_process[thread_ID][counter - 1] == true)
							{//    ж       ʲô  ˼    

								Xhm[thread_ID][0] = initial_offset[0] + x_index*inter_voxel_distance;
								Xhm[thread_ID][1] = initial_offset[1] + y_index*inter_voxel_distance;
								Xhm[thread_ID][2] = initial_offset[2] + z_index*inter_voxel_distance;

								//cout << "Xhm " << Xhm[thread_ID][0] << ", " <<Xhm[thread_ID][1] << ", " << Xhm[thread_ID][2] << endl;


								//char fg; cin >> fg;
								for (int_type_t inner_cam = 0;
										inner_cam < number_cameras; inner_cam++){

									original_in_image = ProjectPointAndReturnIndex(projection_matrices[inner_cam],
											Xhm[thread_ID], x[thread_ID], rows, cols, interacting_p_within_image);

									//cout << "x " << x[thread_ID][0] << " " << x[thread_ID][0] << endl;

//									if (original_in_image){
//										cout << "An original!" << endl;
//									}
									// now, project a point -- this orthogonal try is not quite right ...
									// should be multiplying by an orthogonal vector will give us zero.


									for (int h = 0; h < 3; h++){
										Xhmp[thread_ID][h] = Xhm[thread_ID][h] + V[inner_cam][h];
									}

									in_image = ProjectPointAndReturnIndex(projection_matrices[inner_cam],
											Xhmp[thread_ID], xp[thread_ID], rows, cols, interacting_p_within_image);

									r_squared = SquaredDistance(x[thread_ID], xp[thread_ID], 2);
									// temp for space considerations ....
									r = sqrt(r_squared); ///2.0;
									//some_found = 0;

									if (r <= 0.5){
										x[thread_ID][0] = round(x[thread_ID][0]);
										x[thread_ID][1] = round(x[thread_ID][1]);
										r = 0.5;
									}
									// now walk through ....
									for (r_index = max(0.0, ceil(x[thread_ID][1] - r)), r_max = min(floor(x[thread_ID][1] + r), rows-1.0); r_index <= r_max; r_index++ ){
										// set up the bounds ...
										row_squared_distance_from_center = pow(r_index - x[thread_ID][1], 2);

										if (row_squared_distance_from_center <= r_squared){
											//col_squared_distance_from_center = r_squared - row_squared_distance_from_center;
											//temp_r = sqrt(col_squared_distance_from_center) + 0.25;

											//temp_r = r;
											for (c_index = max(0.0, ceil(x[thread_ID][0] - r)), c_max = min(floor(x[thread_ID][0] + r), cols-1.0); c_index <= c_max; c_index++ ){
												// make a vector ...
												xp[thread_ID][0] = c_index;  xp[thread_ID][1] = r_index;
												//r_squared_test = SquaredDistance(x, xp, 2);
												r_squared_test = pow(x[thread_ID][0] - xp[thread_ID][0], 2) + pow(x[thread_ID][1] - xp[thread_ID][1], 2);

												if (r_squared_test <= r_squared){
													//some_found++;

													interacting_p_within_image = r_index*cols + c_index;
													interacting_p = rows*cols*inner_cam + interacting_p_within_image;
													//is_in_undistorted_image  ά  Ϊnumber_pixels  bool ͵     
													if (is_in_undistorted_image[interacting_p]){//is_in_undistorted_image  Ϊtrue


														{
															terms_arr_per_div[thread_ID][arr0][arr1] = interacting_p;//terms_arr_per_div  ŵ      ر      

															arr1++;
															last_term_position++;//  һ      last_term_position  ++  arr1Ҳ++.

															if (arr1 == term_increment){
																arr1 = 0;
																arr0++;

#pragma omp critical
																{
																	if (terms_arr_per_div[thread_ID].size() == arr0){
																		cout << "Creating new! " << thread_ID << endl;
																		terms_arr_per_div[thread_ID].push_back( new int_type_t[term_increment] );
																	}
																}

																number_terms += term_increment;

															}

															//															if (arr1 != last_term_position % term_increment){
															//																cout << "arr1 and term mismatch .... " << arr1 << ", " << last_term_position << endl;
															//															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}

			//cout << "Old terms size " << number_terms << endl;
			//#pragma omp critical
			{
				cout << "thread id " << thread_ID << endl;
				cout << "current terms size " << number_terms << " and what should be the current terms size " << last_term_position << endl;
				start_indices_terms_per_div[thread_ID][counter] = last_term_position;
			}
			// compute the last slot ...

			//start_index_for_terms_arr[number_voxels_grid] = last_term_position;
		}
	}
}



void NewArrayForParallelOps(vector< int_type_t* >& terms, int_type_t term_increment){
#pragma omp critical
	{
		terms.push_back(new int_type_t[term_increment]);
	}

}



void ReconstructionStructure::ComputeProjectionInformationSpheresArrRepresentationParallelForOnlyThoseMarked(){

	char ch;

	// allocate all of the vectors ahead of time.
	vector< vector<double> > Xhm(number_divisions, vector<double>(4, 1));
	vector< vector<double> > x(number_divisions, vector<double>(3, 1));
	vector< vector<double> > Xhmp(number_divisions, vector<double>(4, 1));
	vector< vector<double> > xp(number_divisions, vector<double>(3, 1));

	// make these private
	int_type_t thread_ID;
	int r_index, c_index, r_max, c_max;
	double row_squared_distance_from_center;
	int_type_t last_term_position = 0;
	int_type_t number_terms = term_increment;
	int_type_t voxel_index = 0;
	double r_squared, r_squared_test, r;
	int_type_t interacting_p, interacting_p_within_image;
	bool in_image, original_in_image;
	int_type_t* arr_ptr;
	int_type_t counter;

	int_type_t arr0, arr1;
	// global to all threads
	int_type_t max_per_thread  = number_voxels_grid/number_divisions;
	int_type_t arr0_established_counter;
	int_type_t arr1_established_counter;



	vector< int_type_t* > start_indices_terms_per_div_temp;
	vector< vector< int_type_t* > > terms_arr_per_div_temp;
	start_indices_terms_per_div_temp.resize(number_divisions, 0);
	terms_arr_per_div_temp.resize(number_divisions, vector<int_type_t*>());

	for (int i = 0; i < number_divisions; i++){
		start_indices_terms_per_div_temp[i] = new int_type_t[max_per_thread + 2]; // to take care fo the ends ...
		terms_arr_per_div_temp[i].push_back(new int_type_t[term_increment]);
	}
	//vector<int_type_t*> arr_ptrs(number_divisions, 0);




	vector< vector<double> > V(number_cameras, vector<double>(3, 0));
	for (int_type_t inner_cam = 0;
			inner_cam < number_cameras; inner_cam++){
		if (plane_equations[inner_cam][0] != 0 || plane_equations[inner_cam][1] != 0){
			V[inner_cam][1] = plane_equations[inner_cam][0];
			V[inner_cam][0] = -plane_equations[inner_cam][1];
			V[inner_cam][2] = 0;
		}	else {
			V[inner_cam][0] = 0; //plane_equations[inner_cam][0];
			V[inner_cam][1] = -plane_equations[inner_cam][2];
			V[inner_cam][2] = plane_equations[inner_cam][1];
		}


		NormalizeVector(V[inner_cam]);
		//MultiplyVectorByScalar(V, inter_voxel_distance/2.0, V, 3);
		MultiplyVectorByScalar(V[inner_cam], inter_voxel_distance/2.0, V[inner_cam], 3);
	}

	//	//int some_found
	//
	//	// change the indices so that start index for terms incorporates the +1 so we don't have to test all of the time.
	//
	// this can be in parallel if we allocate terms ahead of time ... and vey nice b/c then each has its own slot.  Do this later for speed.
	//cout << "After preliminary alloc ... " << endl; cin >> ch;
	omp_set_num_threads(number_divisions);

	int_type_t starting_term_position, ending_term_position;

	// non-parallel for now .....
//	for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
#pragma omp parallel private(arr0, arr1, arr_ptr, thread_ID, counter, r_index, c_index, r_max, c_max, row_squared_distance_from_center, last_term_position, number_terms, voxel_index, r_squared, r_squared_test, r, interacting_p, interacting_p_within_image, in_image, original_in_image, arr0_established_counter, arr1_established_counter, starting_term_position, ending_term_position)

	{
		// preliminaries
		thread_ID = omp_get_thread_num();
		last_term_position = 0;
		number_terms = term_increment;
		voxel_index = 0;
		//		arr_ptr = terms_arr_per_div[thread_ID][0];
		counter = 0;
		arr0 = 0;
		arr1 = 0;
		arr0_established_counter = 0;
		arr1_established_counter = 0;

		//	for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
		{

			//			last_term_position = 0;
			//			number_terms = term_increment;
			//			voxel_index = 0;
			//			arr_ptr = terms_arr_per_div[thread_ID][0];
			//			counter = 0;
			//			arr0 = 0;
			//			arr1 = 0;

			// these probably made this fail ... r_index switching valuesall of the time ...
			for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

#pragma omp critical
				{
					if (x_index % 20 == 0){
						cout << "Computing projection information " << x_index << " out of " << number_voxels_per_dim[0] << " and thread id " << thread_ID << endl;
					}
				}
				for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
					for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){

						if (voxel_index % number_divisions == thread_ID){
							// can do a ptr here ...
							start_indices_terms_per_div_temp[thread_ID][counter] = last_term_position;
							counter++;

							// test here is if we place no restrictions on which ones get b-projected, do we get rid of the problem?
							// this is just copying over the current data ..

							//cout << "VI " << voxel_index << endl;

							// chjanged this here ... probably it ....
							starting_term_position = start_indices_terms_per_div[thread_ID][counter - 1];
							ending_term_position = start_indices_terms_per_div[thread_ID][counter];


							// this loop done if the projections have already been calculated.
							if (starting_term_position != ending_term_position){
								//if (set_by_first_d[voxel_index] == false && multi_dim_voxels_to_process[thread_ID][counter - 1] == false){
								// walk through

								//starting_term_position = start_indices_terms_per_div[thread_ID][counter - 1];
								//last_term_position = start_indices_terms_per_div[thread_ID][counter];

								//cout << "Start, end " << starting_term_position << ", " << ending_term_position << endl;
								for (int_type_t interact_index = starting_term_position;
										interact_index < ending_term_position; interact_index++){

									terms_arr_per_div_temp[thread_ID][arr0][arr1] = terms_arr_per_div[thread_ID][arr0_established_counter][arr1_established_counter];

									// increment the new AND old counters .....
									last_term_position++;
									arr1++;
									arr1_established_counter++;

									//
									//									cout << "interact index ... " << interact_index << " last term index " << ending_term_position << endl;
									//									//									if (arr0 >= 1 || arr1 > 16770000){
									//									cout << "new arr0 and arr1 " << arr0 << ", " << arr1 << "  old arr0 and arr1 " << arr0_established_counter << ", " << arr1_established_counter
									//											<< " term increment is " << term_increment << endl;
									//									//									}

									//									if (interact_index > 2000){
									//										cin >> ch;
									//									}
									// this is when we're changing over on the NEW terms ....
									if (arr1 == term_increment){

										//										cout << " arr0 and arr1 " << arr0 << ", " << arr1 << "  old arr0 and arr1 " << arr0_established_counter << ", " << arr1_established_counter
										//												<< " term increment is " << term_increment << endl;
										//										cin >> ch;

										arr1 = 0;
										arr0++;

										//	cout << "Creating new .... " << endl;
										//	cout << "new arr0 and arr1 " << arr0 << ", " << arr1 << "  old arr0 and arr1 " << arr0_established_counter << ", " << arr1_established_counter
										//		<< " term increment is " << term_increment << endl;
										// this will be critical in the future.
										//#pragma omp critical
										{
											if (terms_arr_per_div_temp[thread_ID].size() == arr0){
												cout << "Creating new, version 0! " << thread_ID << " at voxel " << voxel_index << endl;
												NewArrayForParallelOps(terms_arr_per_div_temp[thread_ID], term_increment);
												//terms_arr_per_div_temp[thread_ID].push_back( new int_type_t[term_increment] );

												//												cout << "test " << endl;
												//												cin >> ch;
												//												terms_arr_per_div_temp[thread_ID][arr0][0] = 1;
											}
										}

										number_terms += term_increment;

									}

									// this is when we're changing over on the OLD terms .... just need to update the counters ....
									// dealloc the old data ....
									if (arr1_established_counter == term_increment){

										cout << "Deleting ... " << thread_ID << endl;
										delete [] terms_arr_per_div[thread_ID][arr0_established_counter];
										terms_arr_per_div[thread_ID][arr0_established_counter] = 0;


										arr1_established_counter = 0;
										arr0_established_counter++;

										//										cout << "Test old " << terms_arr_per_div[thread_ID][arr0_established_counter - 1][arr1_established_counter - 1] << endl;
										//
										//										cout << "After deleting  .... " << endl;
										//										cout << "new arr0 and arr1 " << arr0 << ", " << arr1 << "  old arr0 and arr1 " << arr0_established_counter << ", " << arr1_established_counter
										//																													<< " term increment is " << term_increment << endl;
										//										cin >> ch;
										//										cout << "Test old " << terms_arr_per_div[thread_ID][arr0_established_counter][arr1_established_counter] << endl;
										//										cout << "Test old " << terms_arr_per_div[thread_ID][arr0_established_counter][arr1_established_counter + 1] << endl;
										//																				cin >> ch;


									}
								}
							}	else {

								if (multi_dim_voxels_to_process[thread_ID][counter - 1] == true){
									// now for this one, we're concerned with arr1 and arr0 -- don't touch arr0_establisted and arr1_established ....


									Xhm[thread_ID][0] = initial_offset[0] + x_index*inter_voxel_distance;
									Xhm[thread_ID][1] = initial_offset[1] + y_index*inter_voxel_distance;
									Xhm[thread_ID][2] = initial_offset[2] + z_index*inter_voxel_distance;


									for (int_type_t inner_cam = 0;
											inner_cam < number_cameras; inner_cam++){

										original_in_image = ProjectPointAndReturnIndex(projection_matrices[inner_cam],
												Xhm[thread_ID], x[thread_ID], rows, cols, interacting_p_within_image);

										// now, project a point -- this orthogonal try is not quite right ...
										// should be multiplying by an orthogonal vector will give us zero.


										for (int h = 0; h < 3; h++){
											Xhmp[thread_ID][h] = Xhm[thread_ID][h] + V[inner_cam][h];
										}

										in_image = ProjectPointAndReturnIndex(projection_matrices[inner_cam],
												Xhmp[thread_ID], xp[thread_ID], rows, cols, interacting_p_within_image);

										r_squared = SquaredDistance(x[thread_ID], xp[thread_ID], 2);
										// temp for space considerations ....
										r = sqrt(r_squared); ///2.0;
										//some_found = 0;

										if (r <= 0.5){
											x[thread_ID][0] = round(x[thread_ID][0]);
											x[thread_ID][1] = round(x[thread_ID][1]);
											r = 0.5;
										}
										// now walk through ....
										for (r_index = max(0.0, ceil(x[thread_ID][1] - r)), r_max = min(floor(x[thread_ID][1] + r), rows-1.0); r_index <= r_max; r_index++ ){
											// set up the bounds ...
											row_squared_distance_from_center = pow(r_index - x[thread_ID][1], 2);

											if (row_squared_distance_from_center <= r_squared){
												//col_squared_distance_from_center = r_squared - row_squared_distance_from_center;
												//temp_r = sqrt(col_squared_distance_from_center) + 0.25;

												//temp_r = r;
												for (c_index = max(0.0, ceil(x[thread_ID][0] - r)), c_max = min(floor(x[thread_ID][0] + r), cols-1.0); c_index <= c_max; c_index++ ){
													// make a vector ...
													xp[thread_ID][0] = c_index;  xp[thread_ID][1] = r_index;
													//r_squared_test = SquaredDistance(x, xp, 2);
													r_squared_test = pow(x[thread_ID][0] - xp[thread_ID][0], 2) + pow(x[thread_ID][1] - xp[thread_ID][1], 2);

													if (r_squared_test <= r_squared){
														//some_found++;

														interacting_p_within_image = r_index*cols + c_index;
														interacting_p = rows*cols*inner_cam + interacting_p_within_image;

														if (is_in_undistorted_image[interacting_p]){
															terms_arr_per_div_temp[thread_ID][arr0][arr1] = interacting_p;

															arr1++;
															last_term_position++;

															if (arr1 == term_increment){
																arr1 = 0;
																arr0++;

																cout << "Creating new, version2! " << thread_ID << " at voxel " << voxel_index << endl;
																NewArrayForParallelOps(terms_arr_per_div_temp[thread_ID], term_increment);
																cout << "After creating new " << thread_ID << endl;


																number_terms += term_increment;

															}

															//															if (arr1 != last_term_position % term_increment){
															//																cout << "arr1 and term mismatch .... " << arr1 << ", " << last_term_position << endl;
															//															}
														}
													}
												}
											}
										}
									}

								}
							}
						}
					}
				}
			}

			//cout << "Old terms size " << number_terms << endl;
			//#pragma omp critical
			{
				// we neverget to the last one ...

				if (arr1_established_counter != term_increment){
					delete [] terms_arr_per_div[thread_ID][arr0_established_counter];
				}
				cout << "thread id " << thread_ID << endl;
				cout << "current terms size " << number_terms << " and what should be the current terms size " << last_term_position << endl;
				//cin >> ch;
				start_indices_terms_per_div_temp[thread_ID][counter] = last_term_position;

				//delete [] start_indices_terms_per_div[thread_ID];

			}
			//			// compute the last slot ...

			//start_index_for_terms_arr[number_voxels_grid] = last_term_position;
		}
	}

	start_indices_terms_per_div.swap(start_indices_terms_per_div_temp);

	terms_arr_per_div.swap(terms_arr_per_div_temp);

}
