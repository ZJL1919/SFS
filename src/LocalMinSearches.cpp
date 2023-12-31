/*
 * LocalMinSearches.cpp
 *
 *  Created on: May 1, 2020
 *      Author: Jinglong
 */


#include "ReconstructionStructure.hpp"
#include "Helper0.hpp"
#include <omp.h>
#include <algorithm>


int_type_t ReconstructionStructure::LocalMinSearchParallel6Arr(bool* terms_changed_this_round, bool first){


	// in multi voxels to process, the values are voxel_index/number_divsions + thread_ID


	//  make another local search parallel -- independent threads to deal with all of this stuff?
	char ch;
	int_type_t number_changed = 0;
	int_type_t voxel_index;
	int_type_t last_term_position = 0;
	int_type_t arr_index_0, arr_index_1;
	int_type_t starting_term_position;
	int_type_t term_index;
	int first_d;
	bool current_value;
	bool will_change;


	//	if (terms_changed_this_round == 0){
	//		terms_changed_this_round.resize(coefficients.size(), false);
	//	}	else {
	for (int i = 0, n = coefficients.size(); i < n; i++){
		terms_changed_this_round[i] = false;
	}
	//	}



	//__gnu_parallel::random_shuffle(voxels_to_process.begin(), voxels_to_process.end());
	//__gnu_parallel::random_shuffle(voxels_to_process.begin(), voxels_to_process.end());


	int_type_t start_index, end_index, increment;

	increment = 100;

	vector<int> change_found(number_divisions, false);
	int_type_t thread_ID;
	int_type_t number_to_process_size = number_voxels_grid/number_divisions + 1;
	vector<int_type_t> change_indices(number_divisions, 0);
	bool master_still_test = true;
	int_type_t local_counter_index;
	vector<int_type_t> last_tested(number_divisions, 0);
	int_type_t thread_ID_to_update;
	int_type_t max_number_to_process;
	int_type_t thread_counter;
	int_type_t interact_index;
	int_type_t i;
	vector<bool> local_still_test(number_divisions, true);
	uint64_t old_error = sfis_error;

	/*__gnu_parallel::*/random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
	/*__gnu_parallel::*/random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
	//for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
	//	__gnu_parallel::random_shuffle(multi_voxels_to_process[thread_counter].begin(), multi_voxels_to_process[thread_counter].end());
	//	__gnu_parallel::random_shuffle(multi_voxels_to_process[thread_counter].begin(), multi_voxels_to_process[thread_counter].end());
	//}
	//start_index = 0;  end_index = min(voxels_to_process.size(), increment);

	cout << "Max size " << number_to_process_size << endl;
#pragma omp parallel private(local_counter_index, i, interact_index, current_value, first_d, voxel_index, thread_ID,start_index, end_index, starting_term_position, last_term_position, term_index, arr_index_0, arr_index_1)
	//for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
	{
		// set thread id, start and stops, thread id
		thread_ID = omp_get_thread_num();
		start_index = 0; end_index = min(number_to_process_size, increment);

		//while (local_still_test[thread_ID] == true)
		//#pragma omp barrier
		while (master_still_test == true)
		{

			change_found[thread_ID] = 0;

			for (i = start_index; i < end_index && change_found[thread_ID] == 0 && i < number_to_process_size; i++){
				if (i % 500000 == thread_ID){
					//if (i % 1 == 0){
					//#pragma omp critical
					{
						cout << "Local search " << i << " out of " << number_to_process_size << " in local min search, thread id " << thread_ID << endl;
					}
				}


				if (multi_dim_voxels_to_process[thread_ID][voxel_ordering[i]] == false){
					end_index++;
					last_tested[thread_ID] = i;

				}	else {
					local_counter_index = voxel_ordering[i];
					//local_counter_index = voxel_ordering.at(i);
					voxel_index = local_counter_index*number_divisions + thread_ID;

					if (voxel_index >= number_voxels_grid){
						cout << "Error -- voxel index too high ! " << voxel_index << " out of " << number_voxels_grid << endl;
						exit(1);
					}
					last_tested[thread_ID] = i;

					starting_term_position = start_indices_terms_per_div[thread_ID][local_counter_index];
					last_term_position = start_indices_terms_per_div[thread_ID][local_counter_index + 1];

					current_value = configuration_grid[voxel_index];

					first_d = 0;

					if (last_term_position != starting_term_position){

						arr_index_0 = starting_term_position/term_increment;
						arr_index_1 = starting_term_position % term_increment;


						// if there's fewer than 27 lookups, we don't do the look up.
//						if (last_term_position - starting_term_position > 27){
//
//						}



						if (current_value == true){
							for (interact_index = starting_term_position;
									interact_index < last_term_position; interact_index++){

								term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];

								//cout << "term index " << term_index << " out of " << number_pixels << endl;
								if (number_zeros_grid[term_index] == 0){
									first_d += coefficients[term_index];
								}


								arr_index_1++;
								if (arr_index_1 == term_increment){
									arr_index_0++;
									arr_index_1 = 0;
								}

							}
						}	else {
							//cout << "Impossible!" << endl;
							for (interact_index = starting_term_position;
									interact_index < last_term_position; interact_index++){

								term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
								if (number_zeros_grid[term_index] == 1){
									first_d += coefficients[term_index];
								}

								arr_index_1++;
								if (arr_index_1 == term_increment){
									arr_index_0++;
									arr_index_1 = 0;
								}
							}
						}


						//if ( ((first_d > 0 && current_value) ||  (first_d < 0 && !current_value))){
						if ( ((first_d >= 0 && current_value) ||  (first_d < 0 && !current_value))){
							// is there a different way to indicate that we've found a change, so short circuit the remainder?
							if (current_value){
								change_found[thread_ID] = -1;
								if (first_d == 0){
									change_found[thread_ID] = -2;
								}
							}	else {
								change_found[thread_ID] = 1;
							}



						}
					}
				}
			}
			// done the loop walking through increment number of items in the to-process list.

#pragma omp barrier // all threads stop here.

#pragma omp master
			{
				// we need to pick which one we're going to update... first, are there any changes?
				max_number_to_process = 0;
				thread_ID_to_update = number_divisions;
				master_still_test = false;

				//				cout << "Before! " << endl;
				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << " change " << change_found[thread_counter] << endl;
				//				}


				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
					if (change_found[thread_counter] != 0){

						if (number_to_process_size - last_tested[thread_counter] > max_number_to_process ){
							max_number_to_process = number_to_process_size - last_tested[thread_counter];
							thread_ID_to_update = thread_counter;
						}
						//	cout << "ID : " << thread_counter << " Sizes " <<  number_to_process_size[thread_counter] - last_tested[thread_counter] << endl;
					}	else {
						if (last_tested[thread_counter] < number_to_process_size){
							last_tested[thread_counter]++;
						}
					}
				}

				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
					if (change_found[thread_counter] == -2 || thread_counter == thread_ID_to_update){
						//if (thread_ID_to_update < number_divisions){

						local_counter_index = voxel_ordering[last_tested[thread_counter]];
						voxel_index = local_counter_index* number_divisions + thread_counter;

						//	cout << "Updating voxel index " << voxel_index << endl;

						number_changed++;

						current_value = configuration_grid[voxel_index];

						starting_term_position = start_indices_terms_per_div[thread_counter][local_counter_index];
						last_term_position = start_indices_terms_per_div[thread_counter][local_counter_index + 1];


						arr_index_0 = starting_term_position/term_increment;
						arr_index_1 = starting_term_position % term_increment;

						if (current_value){
							configuration_grid[voxel_index] = false;

							for (interact_index = starting_term_position;
									interact_index < last_term_position; interact_index++){
								term_index = terms_arr_per_div[thread_counter][arr_index_0][arr_index_1];
								number_zeros_grid[term_index]++;
								terms_changed_this_round[term_index] = true;

								arr_index_1++;
								if (arr_index_1 == term_increment){
									arr_index_0++;
									arr_index_1 = 0;
								}
								//cout << "incrementing number zeros grid"
							}

						}	else {
							configuration_grid[voxel_index] = true;

							for (interact_index = starting_term_position;
									interact_index < last_term_position; interact_index++){
								term_index = terms_arr_per_div[thread_counter][arr_index_0][arr_index_1];
								if (number_zeros_grid[term_index] == 0){
									cout << "About to go negative on number of zeros!" << endl;
									cout << "Voxel index, term index " << voxel_index << ", " << term_index << endl;
									exit(1);
								}

								number_zeros_grid[term_index]--;
								terms_changed_this_round[term_index] = true;

								arr_index_1++;
								if (arr_index_1 == term_increment){
									arr_index_0++;
									arr_index_1 = 0;
								}
							}
						}

						last_tested[thread_counter]++;

						//cout << "Number changed " << number_changed << endl;
					}
				}


				// just have master thread do this to avoid start and stop with the threads.
//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//					if (change_found[thread_counter] == -2 || thread_counter == thread_ID_to_update){
//					//if (thread_ID_to_update < number_divisions){
//
//						local_counter_index = voxel_ordering[last_tested[thread_ID_to_update]];
//						voxel_index = local_counter_index* number_divisions + thread_ID_to_update;
//
//						//	cout << "Updating voxel index " << voxel_index << endl;
//
//						number_changed++;
//
//						current_value = configuration_grid[voxel_index];
//
//						starting_term_position = start_indices_terms_per_div[thread_ID_to_update][local_counter_index];
//						last_term_position = start_indices_terms_per_div[thread_ID_to_update][local_counter_index + 1];
//
//
//						arr_index_0 = starting_term_position/term_increment;
//						arr_index_1 = starting_term_position % term_increment;
//
//						if (current_value){
//							configuration_grid[voxel_index] = false;
//
//							for (interact_index = starting_term_position;
//									interact_index < last_term_position; interact_index++){
//								term_index = terms_arr_per_div[thread_ID_to_update][arr_index_0][arr_index_1];
//								number_zeros_grid[term_index]++;
//								terms_changed_this_round[term_index] = true;
//
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//								//cout << "incrementing number zeros grid"
//							}
//
//						}	else {
//							configuration_grid[voxel_index] = true;
//
//							for (interact_index = starting_term_position;
//									interact_index < last_term_position; interact_index++){
//								term_index = terms_arr_per_div[thread_ID_to_update][arr_index_0][arr_index_1];
//								if (number_zeros_grid[term_index] == 0){
//									cout << "About to go negative on number of zeros!" << endl;
//									cout << "Voxel index, term index " << voxel_index << ", " << term_index << endl;
//									exit(1);
//								}
//
//								number_zeros_grid[term_index]--;
//								terms_changed_this_round[term_index] = true;
//
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//							}
//						}
//
//						last_tested[thread_ID_to_update]++;
//
//						//cout << "Number changed " << number_changed << endl;
//					}
//				}

				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
					//cout << "Thread counter " << thread_counter << " last tested " << last_tested[thread_counter] << endl;
					//local_still_test[thread_counter] = false;
					change_found[thread_ID] = 0;
					if (number_to_process_size > last_tested[thread_counter] ){
						master_still_test = true;
					}
				}
				//cout << "Mast still test? " << master_still_test << endl;

				//master_still_test = false;


				//				cout<< "After update!" << endl;
				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << endl;
				//				}
				//				cin >> ch;


				// what is error?
				//								sfis_error = function_constant;
				//								for (int_type_t p = 0; p < number_pixels; p++){
				//									if (number_zeros_grid[p] == 0){
				//										current_pixel_labeling[p] = true;
				//										sfis_error += coefficients[p];
				//									}	else {
				//										current_pixel_labeling[p] = false;
				//									}
				//								}
				//
				//								//cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
				//
				//								if (sfis_error > old_error){
				//									cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
				//									cout<< "Error problem! " << endl;
				//									cin >> ch;
				//								}	else{
				//									old_error = sfis_error;
				//								}
				//cout << "After one round " << endl;
				//cin >> ch;
				//usleep(100);

			} // end master section

#pragma omp barrier

			start_index = last_tested[thread_ID];  end_index = min(number_to_process_size, start_index + increment);
			//#pragma omp barrier
			//#pragma omp barriernumber_to_process_size[thread_ID]

		} // end while
	}

	sfis_error = function_constant;
	for (int_type_t p = 0; p < number_pixels; p++){
		if (number_zeros_grid[p] == 0){
			current_pixel_labeling[p] = true;
			sfis_error += coefficients[p];
		}	else {
			current_pixel_labeling[p] = false;
		}
	}

	cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;


	//cout << "After while loop " <<endl;
	//cin >> ch;

	//cout << "Number changed, " << number_changed << endl;
	//cin >> ch;

	//	// now go through and produce another list ....
	bool should_test_again;

	//vector< vector<int_type_t> >multi_swap(number_divisions, vector<int_type_t>());
#pragma omp parallel for
	for (int_type_t i = 0; i < number_divisions; i++){
		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
			multi_dim_voxels_to_process[i][j] = false;
		}
	}


	int_type_t number_to_test = 0;
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			number_to_test += multi_dim_voxels_to_process[i][j];
//		}
//	}

	cout << "Number to test after zeroing out ... " << number_to_test << endl;


	if (number_changed > 0){
#pragma omp parallel private(thread_ID, local_counter_index, voxel_index, interact_index, should_test_again, last_term_position, starting_term_position, arr_index_0, arr_index_1, term_index)
		{
			thread_ID = omp_get_thread_num();

			for (local_counter_index = 0, voxel_index = thread_ID; voxel_index < number_voxels_grid ; voxel_index = voxel_index + number_divisions, local_counter_index++){
				if (voxel_index % 100000000 == 0){
					cout << "Updating to search list " << voxel_index << " out of " << number_voxels_grid << " in local min search " << endl;
				}

				if (!set_by_first_d[voxel_index]){

					starting_term_position = start_indices_terms_per_div[thread_ID][local_counter_index];
					last_term_position = start_indices_terms_per_div[thread_ID][local_counter_index + 1];

					should_test_again = false;

					arr_index_0 = starting_term_position/term_increment;
					arr_index_1 = starting_term_position % term_increment;
					for (interact_index = starting_term_position;
							interact_index < last_term_position && should_test_again == false; interact_index++){

						term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
						if (terms_changed_this_round[term_index]){
							should_test_again = true;
						}
						arr_index_1++;
						if (arr_index_1 == term_increment){
							arr_index_0++;
							arr_index_1 = 0;
						}
					}
					if (should_test_again){
						//#pragma omp critical
						{
							multi_dim_voxels_to_process[thread_ID][local_counter_index] = true;

						}
					}
				}
			}
		}
	}

	number_to_test = 0;
	for (int_type_t i = 0; i < number_divisions; i++){
		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
			number_to_test += multi_dim_voxels_to_process[i][j];
		}
	}

	cout << "Number to test after going through terms " << number_to_test << endl;
	//	 cin >> ch;

	//	for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
	//
	//	}
	//
	//multi_voxels_to_process.swap(multi_swap);



	sfis_error = function_constant;
	for (int_type_t p = 0; p < number_pixels; p++){
		if (number_zeros_grid[p] == 0){
			current_pixel_labeling[p] = true;
			sfis_error += coefficients[p];
		}	else {
			current_pixel_labeling[p] = false;
		}
	}

	// now, determine which ones are candidates for changing next round ....
	//vector<>


	//return 0;
	return number_changed;
}


//int_type_t ReconstructionStructure::LocalMinSearchParallel6ArrWOHueristic(bool* terms_changed_this_round, bool first){
//
//
//	// in multi voxels to process, the values are voxel_index/number_divsions + thread_ID
//
//
//	//  make another local search parallel -- independent threads to deal with all of this stuff?
//	char ch;
//	int_type_t number_changed = 0;
//	int_type_t voxel_index;
//	int_type_t last_term_position = 0;
//	int_type_t arr_index_0, arr_index_1;
//	int_type_t starting_term_position;
//	int_type_t term_index;
//	int first_d;
//	bool current_value;
//	bool will_change;
//
//
//	//	if (terms_changed_this_round == 0){
//	//		terms_changed_this_round.resize(coefficients.size(), false);
//	//	}	else {
//	for (int i = 0, n = coefficients.size(); i < n; i++){
//		terms_changed_this_round[i] = false;
//	}
//	//	}
//
//
//
//	//__gnu_parallel::random_shuffle(voxels_to_process.begin(), voxels_to_process.end());
//	//__gnu_parallel::random_shuffle(voxels_to_process.begin(), voxels_to_process.end());
//
//
//	int_type_t start_index, end_index, increment;
//
//	increment = 100;
//
//	vector<int> change_found(number_divisions, false);
//	int_type_t thread_ID;
//	int_type_t number_to_process_size = number_voxels_grid/number_divisions + 1;
//	vector<int_type_t> change_indices(number_divisions, 0);
//	bool master_still_test = true;
//	int_type_t local_counter_index;
//	vector<int_type_t> last_tested(number_divisions, 0);
//	int_type_t thread_ID_to_update;
//	int_type_t max_number_to_process;
//	int_type_t thread_counter;
//	int_type_t interact_index;
//	int_type_t i;
//	vector<bool> local_still_test(number_divisions, true);
//	uint64_t old_error = sfis_error;
//
//	__gnu_parallel::random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
//	__gnu_parallel::random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
//	//for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//	//	__gnu_parallel::random_shuffle(multi_voxels_to_process[thread_counter].begin(), multi_voxels_to_process[thread_counter].end());
//	//	__gnu_parallel::random_shuffle(multi_voxels_to_process[thread_counter].begin(), multi_voxels_to_process[thread_counter].end());
//	//}
//	//start_index = 0;  end_index = min(voxels_to_process.size(), increment);
//
//	cout << "Max size " << number_to_process_size << endl;
//#pragma omp parallel private(local_counter_index, i, interact_index, current_value, first_d, voxel_index, thread_ID,start_index, end_index, starting_term_position, last_term_position, term_index, arr_index_0, arr_index_1)
//	//for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
//	{
//		// set thread id, start and stops, thread id
//		thread_ID = omp_get_thread_num();
//		start_index = 0; end_index = min(number_to_process_size, increment);
//
//		//while (local_still_test[thread_ID] == true)
//		//#pragma omp barrier
//		while (master_still_test == true)
//		{
//
//			change_found[thread_ID] = 0;
//
//			for (i = start_index; i < end_index && change_found[thread_ID] == 0 && i < number_to_process_size; i++){
//				if (i % 500000 == thread_ID){
//					//if (i % 1 == 0){
//					//#pragma omp critical
//					{
//						cout << "Local search " << i << " out of " << number_to_process_size << " in local min search, thread id " << thread_ID << endl;
//					}
//				}
//
//
//				if (multi_dim_voxels_to_process[thread_ID][voxel_ordering[i]] == false){
//					end_index++;
//					last_tested[thread_ID] = i;
//
//				}	else {
//					local_counter_index = voxel_ordering[i];
//					//local_counter_index = voxel_ordering.at(i);
//					voxel_index = local_counter_index*number_divisions + thread_ID;
//
//					if (voxel_index >= number_voxels_grid){
//						cout << "Error -- voxel index too high ! " << voxel_index << " out of " << number_voxels_grid << endl;
//						exit(1);
//					}
//					last_tested[thread_ID] = i;
//
//					starting_term_position = start_indices_terms_per_div[thread_ID][local_counter_index];
//					last_term_position = start_indices_terms_per_div[thread_ID][local_counter_index + 1];
//
//					current_value = configuration_grid[voxel_index];
//
//					first_d = 0;
//
//					if (last_term_position != starting_term_position){
//
//						arr_index_0 = starting_term_position/term_increment;
//						arr_index_1 = starting_term_position % term_increment;
//
//
//						// if there's fewer than 27 lookups, we don't do the look up.
////						if (last_term_position - starting_term_position > 27){
////
////						}
//
//
//
//						if (current_value == true){
//							for (interact_index = starting_term_position;
//									interact_index < last_term_position; interact_index++){
//
//								term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//
//								//cout << "term index " << term_index << " out of " << number_pixels << endl;
//								if (number_zeros_grid[term_index] == 0){
//									first_d += coefficients[term_index];
//								}
//
//
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//
//							}
//						}	else {
//							//cout << "Impossible!" << endl;
//							for (interact_index = starting_term_position;
//									interact_index < last_term_position; interact_index++){
//
//								term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//								if (number_zeros_grid[term_index] == 1){
//									first_d += coefficients[term_index];
//								}
//
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//							}
//						}
//
//
//						//if ( ((first_d > 0 && current_value) ||  (first_d < 0 && !current_value))){
//						if ( ((first_d > 0 && current_value) ||  (first_d < 0 && !current_value))){
//							// is there a different way to indicate that we've found a change, so short circuit the remainder?
//							if (current_value){
//								change_found[thread_ID] = -1;
//								if (first_d == 0){
//									change_found[thread_ID] = -2;
//								}
//							}	else {
//								change_found[thread_ID] = 1;
//							}
//
//
//
//						}
//					}
//				}
//			}
//			// done the loop walking through increment number of items in the to-process list.
//
//#pragma omp barrier // all threads stop here.
//
//#pragma omp master
//			{
//				// we need to pick which one we're going to update... first, are there any changes?
//				max_number_to_process = 0;
//				thread_ID_to_update = number_divisions;
//				master_still_test = false;
//
//				//				cout << "Before! " << endl;
//				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << " change " << change_found[thread_counter] << endl;
//				//				}
//
//
//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//					if (change_found[thread_counter] != 0){
//
//						if (number_to_process_size - last_tested[thread_counter] > max_number_to_process ){
//							max_number_to_process = number_to_process_size - last_tested[thread_counter];
//							thread_ID_to_update = thread_counter;
//						}
//						//	cout << "ID : " << thread_counter << " Sizes " <<  number_to_process_size[thread_counter] - last_tested[thread_counter] << endl;
//					}	else {
//						if (last_tested[thread_counter] < number_to_process_size){
//							last_tested[thread_counter]++;
//						}
//					}
//
//					// work on change found ..
//					//change_found[thread_counter] = 0;
//				}
//
//				//				cout<< "After !" << endl;
//				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << endl;
//				//				}
//				//
//				//cout << "thread to update " << thread_ID_to_update << " and dstance " << max_number_to_process << endl;
//				//				cin >> ch;
//
//
//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//					if (change_found[thread_counter] == -2 || thread_counter == thread_ID_to_update){
//						//if (thread_ID_to_update < number_divisions){
//
//						local_counter_index = voxel_ordering[last_tested[thread_counter]];
//						voxel_index = local_counter_index* number_divisions + thread_counter;
//
//						//	cout << "Updating voxel index " << voxel_index << endl;
//
//						number_changed++;
//
//						current_value = configuration_grid[voxel_index];
//
//						starting_term_position = start_indices_terms_per_div[thread_counter][local_counter_index];
//						last_term_position = start_indices_terms_per_div[thread_counter][local_counter_index + 1];
//
//
//						arr_index_0 = starting_term_position/term_increment;
//						arr_index_1 = starting_term_position % term_increment;
//
//						if (current_value){
//							configuration_grid[voxel_index] = false;
//
//							for (interact_index = starting_term_position;
//									interact_index < last_term_position; interact_index++){
//								term_index = terms_arr_per_div[thread_counter][arr_index_0][arr_index_1];
//								number_zeros_grid[term_index]++;
//								terms_changed_this_round[term_index] = true;
//
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//								//cout << "incrementing number zeros grid"
//							}
//
//						}	else {
//							configuration_grid[voxel_index] = true;
//
//							for (interact_index = starting_term_position;
//									interact_index < last_term_position; interact_index++){
//								term_index = terms_arr_per_div[thread_counter][arr_index_0][arr_index_1];
//								if (number_zeros_grid[term_index] == 0){
//									cout << "About to go negative on number of zeros!" << endl;
//									cout << "Voxel index, term index " << voxel_index << ", " << term_index << endl;
//									exit(1);
//								}
//
//								number_zeros_grid[term_index]--;
//								terms_changed_this_round[term_index] = true;
//
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//							}
//						}
//
//						last_tested[thread_counter]++;
//
//						//cout << "Number changed " << number_changed << endl;
//					}
//				}
//
//
//				// just have master thread do this to avoid start and stop with the threads.
////				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
////					if (change_found[thread_counter] == -2 || thread_counter == thread_ID_to_update){
////					//if (thread_ID_to_update < number_divisions){
////
////						local_counter_index = voxel_ordering[last_tested[thread_ID_to_update]];
////						voxel_index = local_counter_index* number_divisions + thread_ID_to_update;
////
////						//	cout << "Updating voxel index " << voxel_index << endl;
////
////						number_changed++;
////
////						current_value = configuration_grid[voxel_index];
////
////						starting_term_position = start_indices_terms_per_div[thread_ID_to_update][local_counter_index];
////						last_term_position = start_indices_terms_per_div[thread_ID_to_update][local_counter_index + 1];
////
////
////						arr_index_0 = starting_term_position/term_increment;
////						arr_index_1 = starting_term_position % term_increment;
////
////						if (current_value){
////							configuration_grid[voxel_index] = false;
////
////							for (interact_index = starting_term_position;
////									interact_index < last_term_position; interact_index++){
////								term_index = terms_arr_per_div[thread_ID_to_update][arr_index_0][arr_index_1];
////								number_zeros_grid[term_index]++;
////								terms_changed_this_round[term_index] = true;
////
////								arr_index_1++;
////								if (arr_index_1 == term_increment){
////									arr_index_0++;
////									arr_index_1 = 0;
////								}
////								//cout << "incrementing number zeros grid"
////							}
////
////						}	else {
////							configuration_grid[voxel_index] = true;
////
////							for (interact_index = starting_term_position;
////									interact_index < last_term_position; interact_index++){
////								term_index = terms_arr_per_div[thread_ID_to_update][arr_index_0][arr_index_1];
////								if (number_zeros_grid[term_index] == 0){
////									cout << "About to go negative on number of zeros!" << endl;
////									cout << "Voxel index, term index " << voxel_index << ", " << term_index << endl;
////									exit(1);
////								}
////
////								number_zeros_grid[term_index]--;
////								terms_changed_this_round[term_index] = true;
////
////								arr_index_1++;
////								if (arr_index_1 == term_increment){
////									arr_index_0++;
////									arr_index_1 = 0;
////								}
////							}
////						}
////
////						last_tested[thread_ID_to_update]++;
////
////						//cout << "Number changed " << number_changed << endl;
////					}
////				}
//
//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//					//cout << "Thread counter " << thread_counter << " last tested " << last_tested[thread_counter] << endl;
//					//local_still_test[thread_counter] = false;
//					change_found[thread_ID] = 0;
//					if (number_to_process_size > last_tested[thread_counter] ){
//						master_still_test = true;
//					}
//				}
//				//cout << "Mast still test? " << master_still_test << endl;
//
//				//master_still_test = false;
//
//
//				//				cout<< "After update!" << endl;
//				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << endl;
//				//				}
//				//				cin >> ch;
//
//
//				// what is error?
//				//								sfis_error = function_constant;
//				//								for (int_type_t p = 0; p < number_pixels; p++){
//				//									if (number_zeros_grid[p] == 0){
//				//										current_pixel_labeling[p] = true;
//				//										sfis_error += coefficients[p];
//				//									}	else {
//				//										current_pixel_labeling[p] = false;
//				//									}
//				//								}
//				//
//				//								//cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
//				//
//				//								if (sfis_error > old_error){
//				//									cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
//				//									cout<< "Error problem! " << endl;
//				//									cin >> ch;
//				//								}	else{
//				//									old_error = sfis_error;
//				//								}
//				//cout << "After one round " << endl;
//				//cin >> ch;
//				//usleep(100);
//
//			} // end master section
//
//#pragma omp barrier
//
//			start_index = last_tested[thread_ID];  end_index = min(number_to_process_size, start_index + increment);
//			//#pragma omp barrier
//			//#pragma omp barriernumber_to_process_size[thread_ID]
//
//		} // end while
//	}
//
//	sfis_error = function_constant;
//	for (int_type_t p = 0; p < number_pixels; p++){
//		if (number_zeros_grid[p] == 0){
//			current_pixel_labeling[p] = true;
//			sfis_error += coefficients[p];
//		}	else {
//			current_pixel_labeling[p] = false;
//		}
//	}
//
//	cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
//
//
//	//cout << "After while loop " <<endl;
//	//cin >> ch;
//
//	//cout << "Number changed, " << number_changed << endl;
//	//cin >> ch;
//
//	//	// now go through and produce another list ....
//	bool should_test_again;
//
//	//vector< vector<int_type_t> >multi_swap(number_divisions, vector<int_type_t>());
//#pragma omp parallel for
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			multi_dim_voxels_to_process[i][j] = false;
//		}
//	}
//
//
//	int_type_t number_to_test = 0;
////	for (int_type_t i = 0; i < number_divisions; i++){
////		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
////			number_to_test += multi_dim_voxels_to_process[i][j];
////		}
////	}
//
//	cout << "Number to test after zeroing out ... " << number_to_test << endl;
//
//
//	if (number_changed > 0){
//#pragma omp parallel private(thread_ID, local_counter_index, voxel_index, interact_index, should_test_again, last_term_position, starting_term_position, arr_index_0, arr_index_1, term_index)
//		{
//			thread_ID = omp_get_thread_num();
//
//			for (local_counter_index = 0, voxel_index = thread_ID; voxel_index < number_voxels_grid ; voxel_index = voxel_index + number_divisions, local_counter_index++){
//				if (voxel_index % 100000000 == 0){
//					cout << "Updating to search list " << voxel_index << " out of " << number_voxels_grid << " in local min search " << endl;
//				}
//
//				if (!set_by_first_d[voxel_index]){
//
//					starting_term_position = start_indices_terms_per_div[thread_ID][local_counter_index];
//					last_term_position = start_indices_terms_per_div[thread_ID][local_counter_index + 1];
//
//					should_test_again = false;
//
//					arr_index_0 = starting_term_position/term_increment;
//					arr_index_1 = starting_term_position % term_increment;
//					for (interact_index = starting_term_position;
//							interact_index < last_term_position && should_test_again == false; interact_index++){
//
//						term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//						if (terms_changed_this_round[term_index]){
//							should_test_again = true;
//						}
//						arr_index_1++;
//						if (arr_index_1 == term_increment){
//							arr_index_0++;
//							arr_index_1 = 0;
//						}
//					}
//					if (should_test_again){
//						//#pragma omp critical
//						{
//							multi_dim_voxels_to_process[thread_ID][local_counter_index] = true;
//
//						}
//					}
//				}
//			}
//		}
//	}
//
//	number_to_test = 0;
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			number_to_test += multi_dim_voxels_to_process[i][j];
//		}
//	}
//
//	cout << "Number to test after going through terms " << number_to_test << endl;
//	//	 cin >> ch;
//
//	//	for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//	//
//	//	}
//	//
//	//multi_voxels_to_process.swap(multi_swap);
//
//
//
//	sfis_error = function_constant;
//	for (int_type_t p = 0; p < number_pixels; p++){
//		if (number_zeros_grid[p] == 0){
//			current_pixel_labeling[p] = true;
//			sfis_error += coefficients[p];
//		}	else {
//			current_pixel_labeling[p] = false;
//		}
//	}
//
//	// now, determine which ones are candidates for changing next round ....
//	//vector<>
//
//
//	//return 0;
//	return number_changed;
//}

//int_type_t ReconstructionStructure::QuicklyAlterZeros(bool* terms_changed_this_round, bool first){
//
//
//	// in multi voxels to process, the values are voxel_index/number_divsions + thread_ID
//
//
//	//  make another local search parallel -- independent threads to deal with all of this stuff?
//	char ch;
//	int_type_t number_changed = 0;
//	int_type_t voxel_index;
//	int_type_t last_term_position = 0;
//	int_type_t arr_index_0, arr_index_1;
//	int_type_t starting_term_position;
//	int_type_t term_index;
//	int first_d;
//	bool current_value;
//	bool will_change;
//
//
//	//	if (terms_changed_this_round == 0){
//	//		terms_changed_this_round.resize(coefficients.size(), false);
//	//	}	else {
//	for (int i = 0, n = coefficients.size(); i < n; i++){
//		terms_changed_this_round[i] = false;
//	}
//	//	}
//
//
//
//	//__gnu_parallel::random_shuffle(voxels_to_process.begin(), voxels_to_process.end());
//	//__gnu_parallel::random_shuffle(voxels_to_process.begin(), voxels_to_process.end());
//
//
//	int_type_t start_index, end_index, increment;
//
//	increment = 100;
//
//	vector<int> change_found(number_divisions, false);
//	int_type_t thread_ID;
//	int_type_t number_to_process_size = number_voxels_grid/number_divisions + 1;
//	vector<int_type_t> change_indices(number_divisions, 0);
//	bool master_still_test = true;
//	int_type_t local_counter_index;
//	vector<int_type_t> last_tested(number_divisions, 0);
//	int_type_t thread_ID_to_update;
//	int_type_t max_number_to_process;
//	int_type_t thread_counter;
//	int_type_t interact_index;
//	int_type_t i;
//	vector<bool> local_still_test(number_divisions, true);
//	uint64_t old_error = sfis_error;
//
//	__gnu_parallel::random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
//	__gnu_parallel::random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
//	//for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//	//	__gnu_parallel::random_shuffle(multi_voxels_to_process[thread_counter].begin(), multi_voxels_to_process[thread_counter].end());
//	//	__gnu_parallel::random_shuffle(multi_voxels_to_process[thread_counter].begin(), multi_voxels_to_process[thread_counter].end());
//	//}
//	//start_index = 0;  end_index = min(voxels_to_process.size(), increment);
//
//	cout << "Max size " << number_to_process_size << endl;
//#pragma omp parallel private(local_counter_index, i, interact_index, current_value, first_d, voxel_index, thread_ID,start_index, end_index, starting_term_position, last_term_position, term_index, arr_index_0, arr_index_1)
//	//for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
//	{
//		// set thread id, start and stops, thread id
//		thread_ID = omp_get_thread_num();
//		start_index = 0; end_index = min(number_to_process_size, increment);
//
//		//while (local_still_test[thread_ID] == true)
//		//#pragma omp barrier
//		while (master_still_test == true)
//		{
//
//			change_found[thread_ID] = 0;
//
//			for (i = start_index; i < end_index && change_found[thread_ID] == 0 && i < number_to_process_size; i++){
//				if (i % 500000 == thread_ID){
//					//if (i % 1 == 0){
//					//#pragma omp critical
//					{
//						cout << "Local search " << i << " out of " << number_to_process_size << " in local min search, thread id " << thread_ID << endl;
//					}
//				}
//
//
//				if (multi_dim_voxels_to_process[thread_ID][voxel_ordering[i]] == false){
//					end_index++;
//					last_tested[thread_ID] = i;
//
//				}	else {
//					local_counter_index = voxel_ordering[i];
//					//local_counter_index = voxel_ordering.at(i);
//					voxel_index = local_counter_index*number_divisions + thread_ID;
//
//					if (voxel_index >= number_voxels_grid){
//						cout << "Error -- voxel index too high ! " << voxel_index << " out of " << number_voxels_grid << endl;
//						exit(1);
//					}
//					last_tested[thread_ID] = i;
//
//					starting_term_position = start_indices_terms_per_div[thread_ID][local_counter_index];
//					last_term_position = start_indices_terms_per_div[thread_ID][local_counter_index + 1];
//
//					current_value = configuration_grid[voxel_index];
//
//					first_d = 0;
//
//					if (last_term_position != starting_term_position){
//
//						arr_index_0 = starting_term_position/term_increment;
//						arr_index_1 = starting_term_position % term_increment;
//
//						if (current_value == true){
//							for (interact_index = starting_term_position;
//									interact_index < last_term_position; interact_index++){
//
//								term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//
//								//cout << "term index " << term_index << " out of " << number_pixels << endl;
//								if (number_zeros_grid[term_index] == 0){
//									first_d += coefficients[term_index];
//								}
//
//
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//
//							}
//						}	else {
//							//cout << "Impossible!" << endl;
//							for (interact_index = starting_term_position;
//									interact_index < last_term_position; interact_index++){
//
//								term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//								if (number_zeros_grid[term_index] == 1){
//									first_d += coefficients[term_index];
//								}
//
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//							}
//						}
//
//
//						//if ( ((first_d > 0 && current_value) ||  (first_d < 0 && !current_value))){
//						if ( ((first_d >= 0 && current_value) ||  (first_d < 0 && !current_value))){
//							// is there a different way to indicate that we've found a change, so short circuit the remainder?
//							if (current_value){
//								change_found[thread_ID] = -1;
//							}	else {
//								change_found[thread_ID] = 1;
//							}
//
//						}
//					}
//				}
//			}
//			// done the loop walking through increment number of items in the to-process list.
//
//#pragma omp barrier // all threads stop here.
//
//#pragma omp master
//			{
//				// we need to pick which one we're going to update... first, are there any changes?
//				max_number_to_process = 0;
//				thread_ID_to_update = number_divisions;
//				master_still_test = false;
//
//				//				cout << "Before! " << endl;
//				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << " change " << change_found[thread_counter] << endl;
//				//				}
//
//
//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//					if (change_found[thread_counter] != 0){
//
//						if (number_to_process_size - last_tested[thread_counter] > max_number_to_process ){
//							max_number_to_process = number_to_process_size - last_tested[thread_counter];
//							thread_ID_to_update = thread_counter;
//						}
//						//	cout << "ID : " << thread_counter << " Sizes " <<  number_to_process_size[thread_counter] - last_tested[thread_counter] << endl;
//					}	else {
//						if (last_tested[thread_counter] < number_to_process_size){
//							last_tested[thread_counter]++;
//						}
//					}
//
//
//
//
//					// work on change found ..
//					change_found[thread_counter] = 0;
//				}
//
//				//				cout<< "After !" << endl;
//				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << endl;
//				//				}
//				//
//				//cout << "thread to update " << thread_ID_to_update << " and dstance " << max_number_to_process << endl;
//				//				cin >> ch;
//
//				// just have master thread do this to avoid start and stop with the threads.
//				if (thread_ID_to_update < number_divisions){
//
//					local_counter_index = voxel_ordering[last_tested[thread_ID_to_update]];
//					voxel_index = local_counter_index* number_divisions + thread_ID_to_update;
//
//					//	cout << "Updating voxel index " << voxel_index << endl;
//
//					number_changed++;
//
//					current_value = configuration_grid[voxel_index];
//
//					starting_term_position = start_indices_terms_per_div[thread_ID_to_update][local_counter_index];
//					last_term_position = start_indices_terms_per_div[thread_ID_to_update][local_counter_index + 1];
//
//
//					arr_index_0 = starting_term_position/term_increment;
//					arr_index_1 = starting_term_position % term_increment;
//
//					if (current_value){
//						configuration_grid[voxel_index] = false;
//
//						for (interact_index = starting_term_position;
//								interact_index < last_term_position; interact_index++){
//							term_index = terms_arr_per_div[thread_ID_to_update][arr_index_0][arr_index_1];
//							number_zeros_grid[term_index]++;
//							terms_changed_this_round[term_index] = true;
//
//							arr_index_1++;
//							if (arr_index_1 == term_increment){
//								arr_index_0++;
//								arr_index_1 = 0;
//							}
//							//cout << "incrementing number zeros grid"
//						}
//
//					}	else {
//						configuration_grid[voxel_index] = true;
//
//						for (interact_index = starting_term_position;
//								interact_index < last_term_position; interact_index++){
//							term_index = terms_arr_per_div[thread_ID_to_update][arr_index_0][arr_index_1];
//							if (number_zeros_grid[term_index] == 0){
//								cout << "About to go negative on number of zeros!" << endl;
//								cout << "Voxel index, term index " << voxel_index << ", " << term_index << endl;
//								exit(1);
//							}
//
//							number_zeros_grid[term_index]--;
//							terms_changed_this_round[term_index] = true;
//
//							arr_index_1++;
//							if (arr_index_1 == term_increment){
//								arr_index_0++;
//								arr_index_1 = 0;
//							}
//						}
//					}
//
//					last_tested[thread_ID_to_update]++;
//
//					//cout << "Number changed " << number_changed << endl;
//				}
//
//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//					//cout << "Thread counter " << thread_counter << " last tested " << last_tested[thread_counter] << endl;
//					//local_still_test[thread_counter] = false;
//					if (number_to_process_size > last_tested[thread_counter] ){
//						master_still_test = true;
//					}
//				}
//				//cout << "Mast still test? " << master_still_test << endl;
//
//				//master_still_test = false;
//
//
//				//				cout<< "After update!" << endl;
//				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << endl;
//				//				}
//				//				cin >> ch;
//
//
//				// what is error?
//				//								sfis_error = function_constant;
//				//								for (int_type_t p = 0; p < number_pixels; p++){
//				//									if (number_zeros_grid[p] == 0){
//				//										current_pixel_labeling[p] = true;
//				//										sfis_error += coefficients[p];
//				//									}	else {
//				//										current_pixel_labeling[p] = false;
//				//									}
//				//								}
//				//
//				//								//cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
//				//
//				//								if (sfis_error > old_error){
//				//									cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
//				//									cout<< "Error problem! " << endl;
//				//									cin >> ch;
//				//								}	else{
//				//									old_error = sfis_error;
//				//								}
//				//cout << "After one round " << endl;
//				//cin >> ch;
//				//usleep(100);
//
//			} // end master section
//
//#pragma omp barrier
//
//			start_index = last_tested[thread_ID];  end_index = min(number_to_process_size, start_index + increment);
//			//#pragma omp barrier
//			//#pragma omp barriernumber_to_process_size[thread_ID]
//
//		} // end while
//	}
//
//	sfis_error = function_constant;
//	for (int_type_t p = 0; p < number_pixels; p++){
//		if (number_zeros_grid[p] == 0){
//			current_pixel_labeling[p] = true;
//			sfis_error += coefficients[p];
//		}	else {
//			current_pixel_labeling[p] = false;
//		}
//	}
//
//	cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
//
//
//	//cout << "After while loop " <<endl;
//	//cin >> ch;
//
//	//cout << "Number changed, " << number_changed << endl;
//	//cin >> ch;
//
//	//	// now go through and produce another list ....
//	bool should_test_again;
//
//	//vector< vector<int_type_t> >multi_swap(number_divisions, vector<int_type_t>());
//#pragma omp parallel for
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			multi_dim_voxels_to_process[i][j] = false;
//		}
//	}
//
//
//	int_type_t number_to_test = 0;
////	for (int_type_t i = 0; i < number_divisions; i++){
////		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
////			number_to_test += multi_dim_voxels_to_process[i][j];
////		}
////	}
//
//	cout << "Number to test after zeroing out ... " << number_to_test << endl;
//
//
//	if (number_changed > 0){
//#pragma omp parallel private(thread_ID, local_counter_index, voxel_index, interact_index, should_test_again, last_term_position, starting_term_position, arr_index_0, arr_index_1, term_index)
//		{
//			thread_ID = omp_get_thread_num();
//
//			for (local_counter_index = 0, voxel_index = thread_ID; voxel_index < number_voxels_grid ; voxel_index = voxel_index + number_divisions, local_counter_index++){
//				if (voxel_index % 100000000 == 0){
//					cout << "Updating to search list " << voxel_index << " out of " << number_voxels_grid << " in local min search " << endl;
//				}
//
//				if (!set_by_first_d[voxel_index]){
//
//					starting_term_position = start_indices_terms_per_div[thread_ID][local_counter_index];
//					last_term_position = start_indices_terms_per_div[thread_ID][local_counter_index + 1];
//
//					should_test_again = false;
//
//					arr_index_0 = starting_term_position/term_increment;
//					arr_index_1 = starting_term_position % term_increment;
//					for (interact_index = starting_term_position;
//							interact_index < last_term_position && should_test_again == false; interact_index++){
//
//						term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//						if (terms_changed_this_round[term_index]){
//							should_test_again = true;
//						}
//						arr_index_1++;
//						if (arr_index_1 == term_increment){
//							arr_index_0++;
//							arr_index_1 = 0;
//						}
//					}
//					if (should_test_again){
//						//#pragma omp critical
//						{
//							multi_dim_voxels_to_process[thread_ID][local_counter_index] = true;
//
//						}
//					}
//				}
//			}
//		}
//	}
//
//	number_to_test = 0;
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			number_to_test += multi_dim_voxels_to_process[i][j];
//		}
//	}
//
//	cout << "Number to test after going through terms " << number_to_test << endl;
//	//	 cin >> ch;
//
//	//	for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//	//
//	//	}
//	//
//	//multi_voxels_to_process.swap(multi_swap);
//
//
//
//	sfis_error = function_constant;
//	for (int_type_t p = 0; p < number_pixels; p++){
//		if (number_zeros_grid[p] == 0){
//			current_pixel_labeling[p] = true;
//			sfis_error += coefficients[p];
//		}	else {
//			current_pixel_labeling[p] = false;
//		}
//	}
//
//	// now, determine which ones are candidates for changing next round ....
//	//vector<>
//
//
//	//return 0;
//	return number_changed;
//}

int_type_t ReconstructionStructure::QuicklyAlterZeros(bool* terms_changed_this_round, bool first){


	// in multi voxels to process, the values are voxel_index/number_divsions + thread_ID


	//  make another local search parallel -- independent threads to deal with all of this stuff?
	char ch;
	int_type_t number_changed = 0;
	int_type_t voxel_index;
	int_type_t last_term_position = 0;
	int_type_t arr_index_0, arr_index_1;
	int_type_t starting_term_position;
	int_type_t term_index;
	int first_d;
	bool current_value;
	bool will_change;


	int_type_t start_index, end_index, increment;

	increment = 100;

	vector<int> change_found(number_divisions, false);
	int_type_t thread_ID;
	int_type_t number_to_process_size = number_voxels_grid/number_divisions + 1;
	vector<int_type_t> change_indices(number_divisions, 0);
	bool master_still_test = true;
	int_type_t local_counter_index;
	vector<int_type_t> last_tested(number_divisions, 0);
	int_type_t thread_ID_to_update;
	int_type_t max_number_to_process;
	int_type_t thread_counter;
	int_type_t interact_index;
	int_type_t i;
	vector<bool> local_still_test(number_divisions, true);
	uint64_t old_error = sfis_error;

	/*__gnu_parallel::*/random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
	/*__gnu_parallel::*/random_shuffle(voxel_ordering.begin(), voxel_ordering.end());

	cout << "Max size " << number_to_process_size << endl;
#pragma omp parallel private(local_counter_index, i, interact_index, current_value, first_d, voxel_index, thread_ID,start_index, end_index, starting_term_position, last_term_position, term_index, arr_index_0, arr_index_1)
	//for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
	{
		// set thread id, start and stops, thread id
		thread_ID = omp_get_thread_num();
		start_index = 0; end_index = min(number_to_process_size, increment);

		//while (local_still_test[thread_ID] == true)
		//#pragma omp barrier
		while (master_still_test == true)
		{

			change_found[thread_ID] = 0;

			for (i = start_index; i < end_index && change_found[thread_ID] == 0 && i < number_to_process_size; i++){
				if (i % 500000 == thread_ID){
					//if (i % 1 == 0){
					//#pragma omp critical
					{
						cout << "Local search " << i << " out of " << number_to_process_size << " in local min search, thread id " << thread_ID << endl;
					}
				}


				if (multi_dim_voxels_to_process[thread_ID][voxel_ordering[i]] == false){
					end_index++;
					last_tested[thread_ID] = i;

				}	else {

						local_counter_index = voxel_ordering[i];
						//local_counter_index = voxel_ordering.at(i);
						voxel_index = local_counter_index*number_divisions + thread_ID;

						if (configuration_grid[voxel_index] == true){

						if (voxel_index >= number_voxels_grid){
							cout << "Error -- voxel index too high ! " << voxel_index << " out of " << number_voxels_grid << endl;
							exit(1);
						}
						last_tested[thread_ID] = i;

						starting_term_position = start_indices_terms_per_div[thread_ID][local_counter_index];
						last_term_position = start_indices_terms_per_div[thread_ID][local_counter_index + 1];

						current_value = true;
						//current_value = configuration_grid[voxel_index];

						first_d = 0;

						if (last_term_position != starting_term_position){

							arr_index_0 = starting_term_position/term_increment;
							arr_index_1 = starting_term_position % term_increment;

							if (current_value == true){
								for (interact_index = starting_term_position;
										interact_index < last_term_position; interact_index++){

									term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];

									//cout << "term index " << term_index << " out of " << number_pixels << endl;
									if (number_zeros_grid[term_index] == 0){
										first_d += coefficients[term_index];
									}


									arr_index_1++;
									if (arr_index_1 == term_increment){
										arr_index_0++;
										arr_index_1 = 0;
									}

								}
							}	else {
								//cout << "Impossible!" << endl;
								for (interact_index = starting_term_position;
										interact_index < last_term_position; interact_index++){

									term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
									if (number_zeros_grid[term_index] == 1){
										first_d += coefficients[term_index];
									}

									arr_index_1++;
									if (arr_index_1 == term_increment){
										arr_index_0++;
										arr_index_1 = 0;
									}
								}
							}


							//if ( ((first_d > 0 && current_value) ||  (first_d < 0 && !current_value))){
							if ( first_d == 0 && current_value ){
								// is there a different way to indicate that we've found a change, so short circuit the remainder?
								if (current_value){
									change_found[thread_ID] = -1;
								}	else {
									change_found[thread_ID] = 1;
								}

							}
						}
					}
				}
			}
			// done the loop walking through increment number of items in the to-process list.

#pragma omp barrier // all threads stop here.

#pragma omp master
			{
				// we need to pick which one we're going to update... first, are there any changes?
				max_number_to_process = 0;
				thread_ID_to_update = number_divisions;
				master_still_test = false;

				//				cout << "Before! " << endl;
				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << " change " << change_found[thread_counter] << endl;
				//				}


				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
					if (change_found[thread_counter] != 0){

						if (number_to_process_size - last_tested[thread_counter] > max_number_to_process ){
							max_number_to_process = number_to_process_size - last_tested[thread_counter];
							thread_ID_to_update = thread_counter;
						}
						//	cout << "ID : " << thread_counter << " Sizes " <<  number_to_process_size[thread_counter] - last_tested[thread_counter] << endl;
					}	else {
						if (last_tested[thread_counter] < number_to_process_size){
							last_tested[thread_counter]++;
						}
					}




					// work on change found ..
					change_found[thread_counter] = 0;
				}

				//				cout<< "After !" << endl;
				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << endl;
				//				}
				//
				//cout << "thread to update " << thread_ID_to_update << " and dstance " << max_number_to_process << endl;
				//				cin >> ch;

				// just have master thread do this to avoid start and stop with the threads.
				if (thread_ID_to_update < number_divisions){

					local_counter_index = voxel_ordering[last_tested[thread_ID_to_update]];
					voxel_index = local_counter_index* number_divisions + thread_ID_to_update;

					//	cout << "Updating voxel index " << voxel_index << endl;

					number_changed++;

					current_value = configuration_grid[voxel_index];

					starting_term_position = start_indices_terms_per_div[thread_ID_to_update][local_counter_index];
					last_term_position = start_indices_terms_per_div[thread_ID_to_update][local_counter_index + 1];


					arr_index_0 = starting_term_position/term_increment;
					arr_index_1 = starting_term_position % term_increment;

					if (current_value){
						configuration_grid[voxel_index] = false;

						for (interact_index = starting_term_position;
								interact_index < last_term_position; interact_index++){
							term_index = terms_arr_per_div[thread_ID_to_update][arr_index_0][arr_index_1];
							number_zeros_grid[term_index]++;
							terms_changed_this_round[term_index] = true;

							arr_index_1++;
							if (arr_index_1 == term_increment){
								arr_index_0++;
								arr_index_1 = 0;
							}
							//cout << "incrementing number zeros grid"
						}

					}	else {
						configuration_grid[voxel_index] = true;

						for (interact_index = starting_term_position;
								interact_index < last_term_position; interact_index++){
							term_index = terms_arr_per_div[thread_ID_to_update][arr_index_0][arr_index_1];
							if (number_zeros_grid[term_index] == 0){
								cout << "About to go negative on number of zeros!" << endl;
								cout << "Voxel index, term index " << voxel_index << ", " << term_index << endl;
								exit(1);
							}

							number_zeros_grid[term_index]--;
							terms_changed_this_round[term_index] = true;

							arr_index_1++;
							if (arr_index_1 == term_increment){
								arr_index_0++;
								arr_index_1 = 0;
							}
						}
					}

					last_tested[thread_ID_to_update]++;

					//cout << "Number changed " << number_changed << endl;
				}

				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
					//cout << "Thread counter " << thread_counter << " last tested " << last_tested[thread_counter] << endl;
					//local_still_test[thread_counter] = false;
					if (number_to_process_size > last_tested[thread_counter] ){
						master_still_test = true;
					}
				}
				//cout << "Mast still test? " << master_still_test << endl;

				//master_still_test = false;


				//				cout<< "After update!" << endl;
				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << endl;
				//				}
				//				cin >> ch;


				// what is error?
				//								sfis_error = function_constant;
				//								for (int_type_t p = 0; p < number_pixels; p++){
				//									if (number_zeros_grid[p] == 0){
				//										current_pixel_labeling[p] = true;
				//										sfis_error += coefficients[p];
				//									}	else {
				//										current_pixel_labeling[p] = false;
				//									}
				//								}
				//
				//								//cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
				//
				//								if (sfis_error > old_error){
				//									cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
				//									cout<< "Error problem! " << endl;
				//									cin >> ch;
				//								}	else{
				//									old_error = sfis_error;
				//								}
				//cout << "After one round " << endl;
				//cin >> ch;
				//usleep(100);

			} // end master section

#pragma omp barrier

			start_index = last_tested[thread_ID];  end_index = min(number_to_process_size, start_index + increment);
			//#pragma omp barrier
			//#pragma omp barriernumber_to_process_size[thread_ID]

		} // end while
	}

	sfis_error = function_constant;
	for (int_type_t p = 0; p < number_pixels; p++){
		if (number_zeros_grid[p] == 0){
			current_pixel_labeling[p] = true;
			sfis_error += coefficients[p];
		}	else {
			current_pixel_labeling[p] = false;
		}
	}

	cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;


	//cout << "After while loop " <<endl;
	//cin >> ch;

	//cout << "Number changed, " << number_changed << endl;
	//cin >> ch;

	//	// now go through and produce another list ....
	bool should_test_again;

	//vector< vector<int_type_t> >multi_swap(number_divisions, vector<int_type_t>());



	//	 cin >> ch;

	//	for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
	//
	//	}
	//
	//multi_voxels_to_process.swap(multi_swap);


// error does not change ....
//	sfis_error = function_constant;
//	for (int_type_t p = 0; p < number_pixels; p++){
//		if (number_zeros_grid[p] == 0){
//			current_pixel_labeling[p] = true;
//			sfis_error += coefficients[p];
//		}	else {
//			current_pixel_labeling[p] = false;
//		}
//	}

	// now, determine which ones are candidates for changing next round ....
	//vector<>


	//return 0;
	return number_changed;
}

//int_type_t ReconstructionStructure::LocalMinSearchParallel6ArrExperimental(bool* terms_changed_this_round, bool first, int_type_t* number_zeros_write){
//
//
//	// in multi voxels to process, the values are voxel_index/number_divsions + thread_ID
//
//
//	//  make another local search parallel -- independent threads to deal with all of this stuff?
//	char ch;
//	int_type_t number_changed = 0;
//	int_type_t voxel_index;
//	int_type_t last_term_position = 0;
//	int_type_t arr_index_0, arr_index_1;
//	int_type_t starting_term_position;
//	int_type_t term_index;
//	int first_d;
//	bool current_value;
//	bool will_change;
//
//
//	//	if (terms_changed_this_round == 0){
//	//		terms_changed_this_round.resize(coefficients.size(), false);
//	//	}	else {
//	for (int i = 0, n = coefficients.size(); i < n; i++){
//		terms_changed_this_round[i] = false;
//	}
//	//	}
//
//
//
//	//__gnu_parallel::random_shuffle(voxels_to_process.begin(), voxels_to_process.end());
//	//__gnu_parallel::random_shuffle(voxels_to_process.begin(), voxels_to_process.end());
//
//
//	int_type_t start_index, end_index, increment;
//
//	increment = 1000;
//
//	vector<int> change_found(number_divisions, false);
//	int_type_t thread_ID;
//	int_type_t number_to_process_size = number_voxels_grid/number_divisions + 1;
//	vector<int_type_t> change_indices(number_divisions, 0);
//	bool master_still_test = true;
//	int_type_t local_counter_index;
//	vector<int_type_t> last_tested(number_divisions, 0);
//	int_type_t thread_ID_to_update;
//	int_type_t max_number_to_process;
//	int_type_t thread_counter;
//	int_type_t interact_index;
//	int_type_t i;
//	vector<bool> local_still_test(number_divisions, true);
//	uint64_t old_error = sfis_error;
//	bool some_change = false;
//
//	__gnu_parallel::random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
//	__gnu_parallel::random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
//	//for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//	//	__gnu_parallel::random_shuffle(multi_voxels_to_process[thread_counter].begin(), multi_voxels_to_process[thread_counter].end());
//	//	__gnu_parallel::random_shuffle(multi_voxels_to_process[thread_counter].begin(), multi_voxels_to_process[thread_counter].end());
//	//}
//	//start_index = 0;  end_index = min(voxels_to_process.size(), increment);
//
//	cout << "Max size " << number_to_process_size << endl;
//	//#pragma omp parallel private(local_counter_index, i, interact_index, current_value, first_d, voxel_index, thread_ID,start_index, end_index, starting_term_position, last_term_position, term_index, arr_index_0, arr_index_1)
//	for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
//	{
//		// set thread id, start and stops, thread id
//		//thread_ID = omp_get_thread_num();
//		start_index = 0; end_index = min(number_to_process_size, increment);
//
//		//while (local_still_test[thread_ID] == true)
//		//#pragma omp barrier
//		while (master_still_test == true)
//		{
//
//			change_found[thread_ID] = 0;
//
//			//for (i = start_index; i < end_index && change_found[thread_ID] == 0 && i < number_to_process_size; i++){
//			//	for (i = start_index; i < end_index && change_found[thread_ID] == 0 && i < number_to_process_size; i++){
//			for (i = start_index; i < end_index && i < number_to_process_size; i++){
//				if (i % 1000 == thread_ID){
//					//if (i % 1 == 0){
//					//#pragma omp critical
//					{
//						cout << "Local search " << i << " out of " << number_to_process_size << " in local min search, thread id " << thread_ID << endl;
//					}
//				}
//
//
//				if (multi_dim_voxels_to_process[thread_ID][voxel_ordering[i]] == false){
//					end_index++;
//					last_tested[thread_ID] = i;
//
//				}	else {
//					local_counter_index = voxel_ordering[i];
//					//local_counter_index = voxel_ordering.at(i);
//					voxel_index = local_counter_index*number_divisions + thread_ID;
//
//					if (voxel_index >= number_voxels_grid){
//						cout << "Error -- voxel index too high ! " << voxel_index << " out of " << number_voxels_grid << endl;
//						exit(1);
//					}
//					last_tested[thread_ID] = i;
//
//					starting_term_position = start_indices_terms_per_div[thread_ID][local_counter_index];
//					last_term_position = start_indices_terms_per_div[thread_ID][local_counter_index + 1];
//
//					current_value = configuration_grid[voxel_index];
//
//					first_d = 0;
//
//					if (last_term_position != starting_term_position){
//
//						arr_index_0 = starting_term_position/term_increment;
//						arr_index_1 = starting_term_position % term_increment;
//
//						if (current_value == true){
//							for (int_type_t interact_index = starting_term_position;
//									interact_index < last_term_position; interact_index++){
//
//								term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//
//								//cout << "term index " << term_index << " out of " << number_pixels << endl;
//								if (number_zeros_grid[term_index] == 0){
//									first_d += coefficients[term_index];
//								}
//
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//
//							}
//						}	else {
//							for (interact_index = starting_term_position;
//									interact_index < last_term_position; interact_index++){
//
//								term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//								if (number_zeros_grid[term_index] == 1){
//									first_d += coefficients[term_index];
//								}
//
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//							}
//						}
//
//
//						if ( ((first_d >= 0 && current_value) ||  (first_d < 0 && !current_value))){
//							// is there a different way to indicate that we've found a change, so short circuit the remainder?
//							// go ahead and change .....
//#pragma omp critical
//							{
//								some_change = true;
//
//								arr_index_0 = starting_term_position/term_increment;
//								arr_index_1 = starting_term_position % term_increment;
//
//								if (current_value){
//									configuration_grid[voxel_index] = false;
//
//									for (interact_index = starting_term_position;
//											interact_index < last_term_position; interact_index++){
//										term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//										number_zeros_write[term_index]++;
//										// do this later?
//										terms_changed_this_round[term_index] = true;
//
//										arr_index_1++;
//										if (arr_index_1 == term_increment){
//											arr_index_0++;
//											arr_index_1 = 0;
//										}
//										//cout << "incrementing number zeros grid"
//									}
//
//								}	else {
//									configuration_grid[voxel_index] = true;
//
//									for (interact_index = starting_term_position;
//											interact_index < last_term_position; interact_index++){
//										term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//										if (number_zeros_write[term_index] == 0){
//											cout << "About to go negative on number of zeros!" << endl;
//											cout << "Voxel index, term index " << voxel_index << ", " << term_index << endl;
//											exit(1);
//										}
//
//										number_zeros_write[term_index]--;
//										terms_changed_this_round[term_index] = true;
//
//										arr_index_1++;
//										if (arr_index_1 == term_increment){
//											arr_index_0++;
//											arr_index_1 = 0;
//										}
//									}
//								}
//								//cout << "Changed " << thread_ID <<  " " << voxel_index << " s, e " << start_index << "   i " << i << " end index " << end_index << endl;
//
//							}
//
//							// we'll stop after one found for each .... ? if good strategy?
//							if (current_value){
//								change_found[thread_ID] = -1;
//							}	else {
//								change_found[thread_ID] = 1;
//							}
//
//							number_changed++;
//
//						}
//					}
//				}
//			}
//			// done the loop walking through increment number of items in the to-process list.
//
//#pragma omp barrier // all threads stop here.
//
//			//			if (some_change == true){
//			//				// flash over the current number zeros ... maybe do this in parallel with a parallel section?
//			//				for (int_type_t j = 0; j < number_pixels; j++){
//			//					if (j % number_divisions == thread_ID){
//			//						number_zeros_grid[j] = number_zeros_write[j];
//			//					}
//			//				}
//			//			}
//
//
//#pragma omp master
//			{
//				//				cout << "grid " << endl;
//				if (some_change == true){
//					// flash over the current number zeros ... maybe do this in parallel with a parallel section?
//					for (int_type_t j = 0; j < number_pixels; j++){
//						number_zeros_grid[j] = number_zeros_write[j];
//					}
//				}
//
//				some_change = false;
//				//			cout << "Updated .... " << endl;
//				//				char ch;
//				//				cin >> ch;
//
//				master_still_test = false;
//
//			} // end master section
//
//#pragma omp barrier
//
//#pragma omp critical
//			{
//				start_index = last_tested[thread_ID] + 1;  end_index = min(number_to_process_size, start_index + increment);
//				if (start_index < end_index){
//					master_still_test = true;
//#pragma omp flush(master_still_test)
//				}
//				//cout << "Master still test? " << master_still_test << endl;
//			}
//
//			//#pragma omp critical
//			//			{
//			//				cout << "start and end next round  " << thread_ID << " s, e " << start_index << " end index " << end_index << endl;
//			//				cout << "Last tested " << last_tested[thread_ID] << endl;
//			//			}
//#pragma omp barrier
//			//#pragma omp barriernumber_to_process_size[thread_ID]
//
//		} // end while
//	}
//
//	sfis_error = function_constant;
//	for (int_type_t p = 0; p < number_pixels; p++){
//		if (number_zeros_grid[p] == 0){
//			current_pixel_labeling[p] = true;
//			sfis_error += coefficients[p];
//		}	else {
//			current_pixel_labeling[p] = false;
//		}
//	}
//
//	cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
//
//
//	//cout << "After while loop " <<endl;
//	//cin >> ch;
//
//	//cout << "Number changed, " << number_changed << endl;
//	//cin >> ch;
//
//	//	// now go through and produce another list ....
//	bool should_test_again;
//
//	//vector< vector<int_type_t> >multi_swap(number_divisions, vector<int_type_t>());
//#pragma omp parallel for
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			multi_dim_voxels_to_process[i][j] = false;
//		}
//	}
//
//
//	int_type_t number_to_test = 0;
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			number_to_test += multi_dim_voxels_to_process[i][j];
//		}
//	}
//
//	cout << "Number to test after zeroing out ... " << number_to_test << endl;
//
//
//	if (number_changed > 0){
//#pragma omp parallel private(thread_ID, local_counter_index, voxel_index, interact_index, should_test_again, last_term_position, starting_term_position, arr_index_0, arr_index_1, term_index)
//		{
//			thread_ID = omp_get_thread_num();
//
//			for (local_counter_index = 0, voxel_index = thread_ID; voxel_index < number_voxels_grid ; voxel_index = voxel_index + number_divisions, local_counter_index++){
//				if (voxel_index % 100000000 == 0){
//					cout << "Updating to search list " << voxel_index << " out of " << number_voxels_grid << " in local min search " << endl;
//				}
//
//				if (!set_by_first_d[voxel_index]){
//
//					starting_term_position = start_indices_terms_per_div[thread_ID][local_counter_index];
//					last_term_position = start_indices_terms_per_div[thread_ID][local_counter_index + 1];
//
//					should_test_again = false;
//
//					arr_index_0 = starting_term_position/term_increment;
//					arr_index_1 = starting_term_position % term_increment;
//					for (interact_index = starting_term_position;
//							interact_index < last_term_position && should_test_again == false; interact_index++){
//
//						term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
//						if (terms_changed_this_round[term_index]){
//							should_test_again = true;
//						}
//						arr_index_1++;
//						if (arr_index_1 == term_increment){
//							arr_index_0++;
//							arr_index_1 = 0;
//						}
//					}
//					if (should_test_again){
//						//#pragma omp critical
//						{
//							multi_dim_voxels_to_process[thread_ID][local_counter_index] = true;
//
//						}
//					}
//				}
//			}
//		}
//	}
//
//	number_to_test = 0;
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			number_to_test += multi_dim_voxels_to_process[i][j];
//		}
//	}
//
//	cout << "Number to test after going through terms " << number_to_test << endl;
//	//	 cin >> ch;
//
//	//	for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//	//
//	//	}
//	//
//	//multi_voxels_to_process.swap(multi_swap);
//
//
//
//	sfis_error = function_constant;
//	for (int_type_t p = 0; p < number_pixels; p++){
//		if (number_zeros_grid[p] == 0){
//			current_pixel_labeling[p] = true;
//			sfis_error += coefficients[p];
//		}	else {
//			current_pixel_labeling[p] = false;
//		}
//	}
//
//	// now, determine which ones are candidates for changing next round ....
//	//vector<>
//
//
//	//return 0;
//	return number_changed;
//}

//int_type_t ReconstructionStructure::LocalMinSearchParallel6ArrTiles(bool* terms_changed_this_round, bool first){
//
//
//	// in multi voxels to process, the values are voxel_index/number_divsions + thread_ID
//
//
//	//  make another local search parallel -- independent threads to deal with all of this stuff?
//	char ch;
//	int_type_t number_changed = 0;
//	int_type_t voxel_index;
//	int_type_t last_term_position = 0;
//	int_type_t arr_index_0, arr_index_1;
//	int_type_t starting_term_position;
//	int_type_t term_index;
//	int first_d;
//	bool current_value;
//	bool will_change;
//
//
//	//	if (terms_changed_this_round == 0){
//	//		terms_changed_this_round.resize(coefficients.size(), false);
//	//	}	else {
//	for (int i = 0, n = coefficients.size(); i < n; i++){
//		terms_changed_this_round[i] = false;
//	}
//
//
//	int_type_t start_index, end_index, increment;
//
//	increment = 100;
//
//	vector<int> change_found(number_divisions, false);
//	int_type_t thread_ID;
//	int_type_t number_to_process_size = number_voxels_grid/number_divisions + 1;
//	vector<int_type_t> change_indices(number_divisions, 0);
//	bool master_still_test = true;
//	int_type_t local_counter_index;
//	vector<int_type_t> last_tested(number_divisions, 0);
//	int_type_t thread_ID_to_update;
//	int_type_t max_number_to_process;
//	int_type_t thread_counter;
//	int_type_t interact_index;
//	int_type_t i;
//	vector<bool> local_still_test(number_divisions, true);
//	uint64_t old_error = sfis_error;
//
//
//	int_type_t number_tiles = start_indices_terms_per_div_tiles.size();
//	cout << "Local search, number of tiles is " << number_tiles << endl;
//	bool tile_found = false;
//	int tile_counter = 0;
//
//	__gnu_parallel::random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
//	__gnu_parallel::random_shuffle(voxel_ordering.begin(), voxel_ordering.end());
//
//
//	cout << "Max size " << number_to_process_size << endl;
//	cout << "Number divisions " << number_divisions << endl;
//
//	// add in the the tile counter, bool tile var in parallel local vars
//#pragma omp parallel private(tile_counter, local_counter_index, i, interact_index, current_value, first_d, voxel_index, thread_ID,start_index, end_index, starting_term_position, last_term_position, term_index, arr_index_0, arr_index_1)
//	//for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
//	{
//		//cout << "Thread id " << thread_ID << endl;
//		// set thread id, start and stops, thread id
//		thread_ID = omp_get_thread_num();
//		start_index = 0; end_index = min(number_to_process_size, increment);
//
//		//while (local_still_test[thread_ID] == true)
//		//#pragma omp barrier
//		while (master_still_test == true)
//		{
//
//			change_found[thread_ID] = 0;
//
//			for (i = start_index; i < end_index && change_found[thread_ID] == 0 && i < number_to_process_size; i++){
//				if (i % 100 == thread_ID){
//					//if (i % 1 == 0){
//					//#pragma omp critical
//					{
//						cout << "Local search " << i << " out of " << number_to_process_size << " in local min search, thread id " << thread_ID << endl;
//					}
//				}
//
//
//				if (multi_dim_voxels_to_process[thread_ID][voxel_ordering[i]] == false){
//					end_index++;
//					last_tested[thread_ID] = i;
//
//				}	else {
//					local_counter_index = voxel_ordering[i];
//					//local_counter_index = voxel_ordering.at(i);
//					voxel_index = local_counter_index*number_divisions + thread_ID;
//
//					if (voxel_index >= number_voxels_grid){
//						cout << "Error -- voxel index too high ! " << voxel_index << " out of " << number_voxels_grid << endl;
//						exit(1);
//					}
//					last_tested[thread_ID] = i;
//
//					// now we need to find it in the tiles ...
//					//tile_found = false;
//					//tile_counter = tile_number[voxel_index];
//					tile_counter = tile_number[thread_ID][local_counter_index];
//					if (tile_counter != -1)
//						//					if (tile_number[voxel_index])
//						//					for (tile_counter = 0; tile_counter < number_tiles && !tile_found; tile_counter++)
//					{
//
//						starting_term_position = start_indices_terms_per_div_tiles[tile_counter][thread_ID][local_counter_index];
//						last_term_position = start_indices_terms_per_div_tiles[tile_counter][thread_ID][local_counter_index + 1];
//
//						//if (last_term_position != starting_term_position)
//						{
//
//							//tile_found = true;
//							current_value = configuration_grid[voxel_index];
//							first_d = 0;
//
//							arr_index_0 = starting_term_position/term_increment;
//							arr_index_1 = starting_term_position % term_increment;
//
//							if (current_value == true){
//								for (int_type_t interact_index = starting_term_position;
//										interact_index < last_term_position; interact_index++){
//
//									term_index = terms_arr_per_div_tiles[tile_counter][thread_ID][arr_index_0][arr_index_1];
//
//									//cout << "term index " << term_index << " out of " << number_pixels << endl;
//									if (number_zeros_grid[term_index] == 0){
//										first_d += coefficients[term_index];
//									}
//
//									arr_index_1++;
//									if (arr_index_1 == term_increment){
//										arr_index_0++;
//										arr_index_1 = 0;
//									}
//
//								}
//							}	else {
//								for (interact_index = starting_term_position;
//										interact_index < last_term_position; interact_index++){
//
//									term_index = terms_arr_per_div_tiles[tile_counter][thread_ID][arr_index_0][arr_index_1];
//									if (number_zeros_grid[term_index] == 1){
//										first_d += coefficients[term_index];
//									}
//
//									arr_index_1++;
//									if (arr_index_1 == term_increment){
//										arr_index_0++;
//										arr_index_1 = 0;
//									}
//								}
//							}
//
//
//							//							if (!current_value){
//							//								first_d  = -1;
//							//							}	else {
//							//								first_d = -1;
//							//							}
//
//							if ( ((first_d >= 0 && current_value) ||  (first_d < 0 && !current_value))){
//								// is there a different way to indicate that we've found a change, so short circuit the remainder?
//								if (current_value){
//									change_found[thread_ID] = -1;
//								}	else {
//									change_found[thread_ID] = 1;
//								}
//
//							}
//						}
//					}
//				}
//			}
//			// done the loop walking through increment number of items in the to-process list.
//
//#pragma omp barrier // all threads stop here.
//
//#pragma omp master
//			{
//				// we need to pick which one we're going to update... first, are there any changes?
//				max_number_to_process = 0;
//				thread_ID_to_update = number_divisions;
//				master_still_test = false;
//
//				//				cout << "Before! " << endl;
//				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << " change " << change_found[thread_counter] << endl;
//				//				}
//
//
//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//					if (change_found[thread_counter] != 0){
//
//						if (number_to_process_size - last_tested[thread_counter] > max_number_to_process ){
//							max_number_to_process = number_to_process_size - last_tested[thread_counter];
//							thread_ID_to_update = thread_counter;
//						}
//						//	cout << "ID : " << thread_counter << " Sizes " <<  number_to_process_size[thread_counter] - last_tested[thread_counter] << endl;
//					}	else {
//						if (last_tested[thread_counter] < number_to_process_size){
//							last_tested[thread_counter]++;
//						}
//					}
//
//
//
//
//					// work on change found ..
//					change_found[thread_counter] = 0;
//				}
//
//				//				cout<< "After !" << endl;
//				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << endl;
//				//				}
//				//
//				//cout << "thread to update " << thread_ID_to_update << " and dstance " << max_number_to_process << endl;
//				//				cin >> ch;
//
//				// just have master thread do this to avoid start and stop with the threads.
//				if (thread_ID_to_update < number_divisions){
//
//					local_counter_index = voxel_ordering[last_tested[thread_ID_to_update]];
//					voxel_index = local_counter_index* number_divisions + thread_ID_to_update;
//
//					//	cout << "Updating voxel index " << voxel_index << endl;
//
//					number_changed++;
//
//
//					//tile_found = false;
//
//					//tile_counter = tile_number[voxel_index];
//					tile_counter = tile_number[thread_ID][local_counter_index];
//
//					if (tile_counter != -1)
//						//for (tile_counter = 0; tile_counter < number_tiles && !tile_found; tile_counter++)
//					{
//
//						starting_term_position = start_indices_terms_per_div_tiles[tile_counter][thread_ID_to_update][local_counter_index];
//						last_term_position = start_indices_terms_per_div_tiles[tile_counter][thread_ID_to_update][local_counter_index + 1];
//
//						//if (starting_term_position != last_term_position)
//						{
//							//tile_found = true;
//
//							current_value = configuration_grid[voxel_index];
//
//							arr_index_0 = starting_term_position/term_increment;
//							arr_index_1 = starting_term_position % term_increment;
//
//							if (current_value){
//								configuration_grid[voxel_index] = false;
//
//								for (interact_index = starting_term_position;
//										interact_index < last_term_position; interact_index++){
//									term_index = terms_arr_per_div_tiles[tile_counter][thread_ID_to_update][arr_index_0][arr_index_1];
//									number_zeros_grid[term_index]++;
//									terms_changed_this_round[term_index] = true;
//
//									arr_index_1++;
//									if (arr_index_1 == term_increment){
//										arr_index_0++;
//										arr_index_1 = 0;
//									}
//									//cout << "incrementing number zeros grid"
//								}
//
//							}	else {
//								configuration_grid[voxel_index] = true;
//
//								for (interact_index = starting_term_position;
//										interact_index < last_term_position; interact_index++){
//									term_index = terms_arr_per_div_tiles[tile_counter][thread_ID_to_update][arr_index_0][arr_index_1];
//									if (number_zeros_grid[term_index] == 0){
//										cout << "About to go negative on number of zeros!" << endl;
//										cout << "Voxel index, term index " << voxel_index << ", " << term_index << endl;
//										exit(1);
//									}
//
//									number_zeros_grid[term_index]--;
//									terms_changed_this_round[term_index] = true;
//
//									arr_index_1++;
//									if (arr_index_1 == term_increment){
//										arr_index_0++;
//										arr_index_1 = 0;
//									}
//								}
//							}
//
//							last_tested[thread_ID_to_update]++;
//						}
//					}
//
//					//cout << "Number changed " << number_changed << endl;
//				}
//
//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//					//cout << "Thread counter " << thread_counter << " last tested " << last_tested[thread_counter] << endl;
//					//local_still_test[thread_counter] = false;
//					if (number_to_process_size > last_tested[thread_counter] ){
//						master_still_test = true;
//					}
//				}
//				//cout << "Mast still test? " << master_still_test << endl;
//
//				//master_still_test = false;
//
//
//				//				cout<< "After update!" << endl;
//				//				for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//				//					cout << "thread number " << thread_counter << " and last tested " << last_tested[thread_counter] << endl;
//				//				}
//				//				cin >> ch;
//
//
//				// what is error?
//				//								sfis_error = function_constant;
//				//								for (int_type_t p = 0; p < number_pixels; p++){
//				//									if (number_zeros_grid[p] == 0){
//				//										current_pixel_labeling[p] = true;
//				//										sfis_error += coefficients[p];
//				//									}	else {
//				//										current_pixel_labeling[p] = false;
//				//									}
//				//								}
//				//
//				//								//cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
//				//
//				//								if (sfis_error > old_error){
//				//									cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
//				//									cout<< "Error problem! " << endl;
//				//									cin >> ch;
//				//								}	else{
//				//									old_error = sfis_error;
//				//								}
//				//cout << "After one round " << endl;
//				//cin >> ch;
//				//usleep(100);
//
//			} // end master section
//
//#pragma omp barrier
//
//			start_index = last_tested[thread_ID];  end_index = min(number_to_process_size, start_index + increment);
//			//#pragma omp barrier
//			//#pragma omp barriernumber_to_process_size[thread_ID]
//
//		} // end while
//	}
//
//	sfis_error = function_constant;
//	for (int_type_t p = 0; p < number_pixels; p++){
//		if (number_zeros_grid[p] == 0){
//			current_pixel_labeling[p] = true;
//			sfis_error += coefficients[p];
//		}	else {
//			current_pixel_labeling[p] = false;
//		}
//	}
//
//	cout << " Error  " << FormatWithCommas<int>(sfis_error/127) << "   number changed " << number_changed << endl;
//
//
//	//cout << "After while loop " <<endl;
//	//cin >> ch;
//
//	//cout << "Number changed, " << number_changed << endl;
//	//cin >> ch;
//
//	//	// now go through and produce another list ....
//	bool should_test_again;
//
//	//vector< vector<int_type_t> >multi_swap(number_divisions, vector<int_type_t>());
//#pragma omp parallel for
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			multi_dim_voxels_to_process[i][j] = false;
//		}
//	}
//
//
//	int_type_t number_to_test = 0;
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			number_to_test += multi_dim_voxels_to_process[i][j];
//		}
//	}
//
//	cout << "Number to test after zeroing out ... " << number_to_test << endl;
//
//
//	if (number_changed > 0){
//
//		// add in the the tile counter, bool tile var in parallel local vars
//
//#pragma omp parallel private(tile_counter, thread_ID, local_counter_index, voxel_index, interact_index, should_test_again, last_term_position, starting_term_position, arr_index_0, arr_index_1, term_index)
//		//		for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
//		{
//			thread_ID = omp_get_thread_num();
//
//			for (local_counter_index = 0, voxel_index = thread_ID; voxel_index < number_voxels_grid ; voxel_index = voxel_index + number_divisions, local_counter_index++){
//				if (voxel_index % 100000000 == 0){
//					cout << "Updating to search list " << voxel_index << " out of " << number_voxels_grid << " in local min search " << endl;
//				}
//
//				if (!set_by_first_d[voxel_index]){
//
//					tile_found = false;
//
//					//tile_counter = tile_number[voxel_index];
//					tile_counter = tile_number[thread_ID][local_counter_index];
//
//					if (tile_counter != -1)
//						//for (tile_counter = 0; tile_counter < number_tiles && !tile_found; tile_counter++)
//					{
//
//						starting_term_position = start_indices_terms_per_div_tiles[tile_counter][thread_ID][local_counter_index];
//						last_term_position = start_indices_terms_per_div_tiles[tile_counter][thread_ID][local_counter_index + 1];
//
//						//if (starting_term_position != last_term_position)
//						{
//							//tile_found = true;
//							should_test_again = false;
//
//							arr_index_0 = starting_term_position/term_increment;
//							arr_index_1 = starting_term_position % term_increment;
//							for (interact_index = starting_term_position;
//									interact_index < last_term_position && should_test_again == false; interact_index++){
//
//								term_index = terms_arr_per_div_tiles[tile_counter][thread_ID][arr_index_0][arr_index_1];
//								if (terms_changed_this_round[term_index]){
//									should_test_again = true;
//								}
//								arr_index_1++;
//								if (arr_index_1 == term_increment){
//									arr_index_0++;
//									arr_index_1 = 0;
//								}
//							}
//							if (should_test_again){
//								//#pragma omp critical
//								{
//									multi_dim_voxels_to_process[thread_ID][local_counter_index] = true;
//
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//
//	number_to_test = 0;
//	for (int_type_t i = 0; i < number_divisions; i++){
//		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
//			number_to_test += multi_dim_voxels_to_process[i][j];
//		}
//	}
//
//	cout << "Number to test after going through terms " << number_to_test << endl;
//	//	 cin >> ch;
//
//	//	for (thread_counter = 0; thread_counter < number_divisions; thread_counter++ ){
//	//
//	//	}
//	//
//	//multi_voxels_to_process.swap(multi_swap);
//
//
//
//	sfis_error = function_constant;
//	for (int_type_t p = 0; p < number_pixels; p++){
//		if (number_zeros_grid[p] == 0){
//			current_pixel_labeling[p] = true;
//			sfis_error += coefficients[p];
//		}	else {
//			current_pixel_labeling[p] = false;
//		}
//	}
//
//	// now, determine which ones are candidates for changing next round ....
//	//vector<>
//
//
//	//return 0;
//	return number_changed;
//}
