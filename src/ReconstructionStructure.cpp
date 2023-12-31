/*
 * ReconstructionStructure.cpp
 *
 *  Created on: July 22, 2020
 *      Author: Jinglong
 */


#include "ReconstructionStructure.hpp"
#include "Helper0.hpp"
#include <omp.h>

ReconstructionStructure::ReconstructionStructure(){
	number_voxels_grid = 0;
	inter_plane_space = -1;
}

ReconstructionStructure::~ReconstructionStructure(){

}

void ReconstructionStructure::DeallocBigArrays(){


	if (start_indices_terms_per_div.size() > 0){
		for (int_type_t i = 0; i < start_indices_terms_per_div.size(); i++){
			delete [] start_indices_terms_per_div[i];
		}

		for (int j = 0; j < number_divisions; j++){
			for (int_type_t i = 0; i < terms_arr_per_div[j].size(); i++){
				delete [] terms_arr_per_div[j][i];
			}
		}
	}


	if (number_terms_per_div.size() > 0){
		cout << "Dealloc! " << endl;
	}

	if (tile_number.size() != 0){
		for (int_type_t i = 0; i < number_divisions; i++){
			delete [] tile_number[i];
		}
	}
}


void ReconstructionStructure::CreateGridMinimal(vector< vector<double> >& boundingbox, double division){
	inter_voxel_distance = division;
	BB = boundingbox;
	initial_offset.resize(3, 0);
	number_voxels_per_dim.resize(3, 0);

	initial_offset[0] = BB[0][0] + inter_voxel_distance/2.0;
	initial_offset[1] = BB[0][1] + inter_voxel_distance/2.0;
	initial_offset[2] = BB[0][2] + inter_voxel_distance/2.0;

	for (int i = 0; i < 3; i++){
		number_voxels_per_dim[i] = floor((BB[1][i] - BB[0][i])/inter_voxel_distance);
		cout << "Number voxels per " << number_voxels_per_dim[i] << endl;
	}
	//number_voxels_per_dim[1] = floor((BB.max(1) - BB.min(1))/inter_voxel_distance);
	//number_voxels_per_dim[2] = floor((BB.max(2) - BB.min(2))/inter_voxel_distance);

	//	number_voxels_per_dim[0] = floor((BB.max(0) - BB.min(0))/inter_voxel_distance) - 1;
	//	number_voxels_per_dim[1] = floor((BB.max(1) - BB.min(1))/inter_voxel_distance);
	//	number_voxels_per_dim[2] = floor((BB.max(2) - BB.min(2))/inter_voxel_distance);

	// configuration grid is only for the voxels.
	number_voxels_grid = number_voxels_per_dim[0]*number_voxels_per_dim[1]*number_voxels_per_dim[2];

	// later - make faster, make an array for ALL of these things ...
	configuration_grid.resize(number_voxels_grid, true);

	//number_of_positive_terms.resize(number_voxels_grid, 0);
	set_by_first_d.resize(number_voxels_grid, false);

	// can make all of these arrays ....
	//	start_index_for_terms_arr = new int_type_t[number_voxels_grid + 1];
}

void ReconstructionStructure::AddProjectionMatrices(vector<ImageData*>& cameras){
	number_cameras = cameras.size();
	rows = cameras[0]->imcv.rows;
	cols = cameras[0]->imcv.cols;

			Matrix P(3, 4);
	Matrix Pinv;

	vector<vector<double> > Pvv(3, vector<double>(4, 0));
	vector<vector<double> > Pivv(3, vector<double>(3, 0));
	vector<double> Cvv(4, 1);
	vector<double> Dvv(3);
	vector<double> lambDvv(3);
	vector<double> Xvv(3);
	vector<double> xvv(3);
	vector<double> Rvv(3);
	vector<double> plane_equation(4);


	// do projection matrices....
	for (int cam = 0; cam < number_cameras; cam++){

		for (int i = 0; i < 3; i++){
			P.Row(i+1) << cameras[cam]->P.at<double>(i, 0) << cameras[cam]->P.at<double>(i, 1) << cameras[cam]->P.at<double>(i, 2) << cameras[cam]->P.at<double>(i, 3);
		}

//		// only first two rows ....
//		for (int i = 0; i < 2; i++){
//			for (int j = 0; j < 4; j++){
//				P(i+1, j+1) /= double(downsample);
//			}
//		}

		//cout << " P " << endl << P << endl;
		//cout << "Original representation " << endl;
		//MatrixPrint(cameras[cam]->P);
		Pinv = P.t() * (P*P.t()).i();

		Pinv = Pinv /Pinv (4, 3);

		//cout << "Pinv " << endl << Pinv << endl;


		// do loop to read in ....
		for (int i = 0; i < 3; i++){
			Cvv[i] = (cameras[cam]->C.at<double>(i, 0));
			for (int c = 0; c < 3; c++){
				Pivv[i][c] = Pinv(i + 1, c + 1);
			}

			for (int c = 0; c < 4; c++){
				Pvv[i][c] = P(i + 1, c + 1);
			}
		}

		//cout << "New mat projection matrices " << endl << P << endl;


				cout << "Projection matrices ... " << endl;
				for (int i = 0; i < 3; i++){

					for (int c = 0; c < 4; c++){
						cout << Pvv[i][c] << " "; //= P(i + 1, c + 1);
					}
					cout << endl;
				}

		projection_matrices.push_back(Pvv);
		inverse_projection_matrices.push_back(Pivv);
		camera_centers.push_back(Cvv);

		plane_equation[0] = cameras[cam]->P.at<double>(2, 0);
		plane_equation[1] = cameras[cam]->P.at<double>(2, 1);
		plane_equation[2] = cameras[cam]->P.at<double>(2, 2);
		plane_equation[3] = cameras[cam]->P.at<double>(2, 3);

		NormalizePlane(plane_equation);
		plane_equations.push_back(plane_equation);


		//		plane_equation[0] = cameras[cam]->PrincipalPlane.a();
		//		plane_equation[1] = cameras[cam]->PrincipalPlane.b();
		//		plane_equation[2] = cameras[cam]->PrincipalPlane.c();
		//		plane_equation[3] = cameras[cam]->PrincipalPlane.d();
		//		NormalizePlane(plane_equation);
		//		plane_equations.push_back(plane_equation);

	}
}

void ReconstructionStructure::CreateGrid(vector< vector<double> >& boundingbox, vector<ImageData*>& cameras, double division, int downsample){

	//	vector<bool> configuration_grid;
	//		vector<bool> part_of_visual_hull;
	//
	//		vector<double> initial_offset;

	inter_voxel_distance = division;
	BB = boundingbox;
	initial_offset.resize(3, 0);
	number_voxels_per_dim.resize(3, 0);

	initial_offset[0] = BB[0][0] + inter_voxel_distance/2.0;
	initial_offset[1] = BB[0][1] + inter_voxel_distance/2.0;
	initial_offset[2] = BB[0][2] + inter_voxel_distance/2.0;

	for (int i = 0; i < 3; i++){
		number_voxels_per_dim[i] = floor((BB[1][i] - BB[0][i])/inter_voxel_distance);
		cout << "Number voxels per " << number_voxels_per_dim[i] << endl;
	}


	// configuration grid is only for the voxels.
	number_voxels_grid = number_voxels_per_dim[0]*number_voxels_per_dim[1]*number_voxels_per_dim[2];

	// later - make faster, make an array for ALL of these things ...
	configuration_grid.resize(number_voxels_grid, true);//bool ͵     

	set_by_first_d.resize(number_voxels_grid, false);//bool      


	term_increment = pow(2, 24);//    
	//term_increment = pow(2, 16);

	number_cameras = cameras.size();// м   ͼ  

	rows = cameras[0]->imcv.rows/downsample;//  0  ͼ             
	cols = cameras[0]->imcv.cols/downsample;

	number_divisions = omp_get_max_threads();//   Դ         ߳     

	number_pixels = number_cameras*rows*cols;//    ͼ  һ   ж      ص 

	current_pixel_labeling.resize(number_pixels, true);
	is_in_undistorted_image.resize(number_pixels, true);
	number_zeros_grid.resize(number_pixels, 0);

	for (int_type_t i = 0; i < number_divisions; i++){
		multi_dim_voxels_to_process.push_back(new bool[number_voxels_grid/number_divisions  + 1]);
		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
			multi_dim_voxels_to_process[i][j] = false;
		}
	}

	voxel_ordering.resize(number_voxels_grid/number_divisions  + 1, 0);
	for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
		voxel_ordering[j] = j;
	}

	//cout << "Midway through ... " << endl;
	//char ch; cin >> ch;


	Matrix P(3, 4);
	Matrix Pinv;

	vector<vector<double> > Pvv(3, vector<double>(4, 0));
	vector<vector<double> > Pivv(3, vector<double>(3, 0));
	vector<double> Cvv(4, 1);
	vector<double> Dvv(3);
	vector<double> lambDvv(3);
	vector<double> Xvv(3);
	vector<double> xvv(3);
	vector<double> Rvv(3);
	vector<double> plane_equation(4);


	// do projection matrices....
	for (int cam = 0; cam < number_cameras; cam++){

		for (int i = 0; i < 3; i++){
			P.Row(i+1) << cameras[cam]->P.at<double>(i, 0) << cameras[cam]->P.at<double>(i, 1) << cameras[cam]->P.at<double>(i, 2) << cameras[cam]->P.at<double>(i, 3);
		}
		// only first two rows ....
		for (int i = 0; i < 2; i++){
			for (int j = 0; j < 4; j++){
				P(i+1, j+1) /= double(downsample);
			}
		}//downsample

		//cout << " P " << endl << P << endl;
		//cout << "Original representation " << endl;
		//MatrixPrint(cameras[cam]->P);
		Pinv = P.t() * (P*P.t()).i();
		Pinv = Pinv /Pinv (4, 3);

		//cout << "Pinv " << endl << Pinv << endl;


		// do loop to read in ....
		for (int i = 0; i < 3; i++){
			Cvv[i] = (cameras[cam]->C.at<double>(i, 0));
			for (int c = 0; c < 3; c++){
				Pivv[i][c] = Pinv(i + 1, c + 1);
			}

			for (int c = 0; c < 4; c++){
				Pvv[i][c] = P(i + 1, c + 1);
			}
		}

		//cout << "New mat projection matrices " << endl << P << endl;


		//		cout << "Projection matrices ... " << endl;
		//		for (int i = 0; i < 3; i++){
		//
		//			for (int c = 0; c < 4; c++){
		//				cout << Pvv[i][c] << " "; //= P(i + 1, c + 1);
		//			}
		//			cout << endl;
		//		}

		projection_matrices.push_back(Pvv);
		inverse_projection_matrices.push_back(Pivv);
		camera_centers.push_back(Cvv);//    

		plane_equation[0] = cameras[cam]->P.at<double>(2, 0);
		plane_equation[1] = cameras[cam]->P.at<double>(2, 1);
		plane_equation[2] = cameras[cam]->P.at<double>(2, 2);
		plane_equation[3] = cameras[cam]->P.at<double>(2, 3);//    

		NormalizePlane(plane_equation);
		plane_equations.push_back(plane_equation);

	}


	int index;
	// for now, assume that all images are in color.  Technically, we can deallocate the undistortion map after the next loop ....
	bool color = true;
	bool to_test;//    

	//cout << "rows, cols " << rows << ", " <<


	cv::Mat undistort_map(rows, cols, CV_8UC1, cv::Scalar(0,0,0));
	for (int_type_t cam = 0; cam < number_cameras; cam++){


		cout << "original map " << cameras[cam]->undistorted_map_original.rows << " " << cameras[cam]->undistorted_map_original.cols << endl;
		cout << "Resize map " << undistort_map.rows << ", " << undistort_map.cols << endl;

		color = cameras[cam]->undistorted_map_original.channels() == 3;
		uchar* p = 0;
		// resize for each time ...
		if (downsample > 1){
			cv::resize(cameras[cam]->undistorted_map_original, undistort_map, undistort_map.size());
			p = &undistort_map.data[0];
		}	else {
			p = &cameras[cam]->undistorted_map_original.data[0];
		}

		if (color){
			for (int r = 0; r < rows; r++){
				for (int c = 0; c < cols; c++){

					to_test = *p > 150;

					if (!to_test){
						// we don't want to look this up for the majority of pixels that are "inside"
						// can also use an iterator -- speed?
						is_in_undistorted_image[cam*(rows*cols) + r*cols + c] = to_test;
					}

					// increment by three ....
					p++;
					p++;
					p++;
				}
			}


		}	else {
			//cout << "Implement for grayscale! " << endl;
			for (int r = 0; r < rows; r++){
				for (int c = 0; c < cols; c++){

					to_test = *p > 150;
					
					//    ֵ>150  is_in_undistorted_imageΪtrue,    Ϊfalse.

					if (!to_test){
						// we don't want to look this up for the majority of pixels that are "inside"
						// can also use an iterator -- speed?
						is_in_undistorted_image[cam*(rows*cols) + r*cols + c] = to_test;
					}

					// increment by one b/c grayscale ....
					p++;
				}
			}

		}
		
	}

}

void ReconstructionStructure::CreateGridUsingAncestor(ReconstructionStructure* RS, vector< vector<double> >& boundingbox, vector<ImageData*>& cameras, double division, int downsample){

	//	vector<bool> configuration_grid;
	//		vector<bool> part_of_visual_hull;
	//
	//		vector<double> initial_offset;

	inter_voxel_distance = division;//  ʾ   ؿ ߳ 
	BB = boundingbox; 
	initial_offset.resize(3, 0);// 洢   ½    ؿ    ĵ      
	number_voxels_per_dim.resize(3, 0);

	initial_offset[0] = BB[0][0] + inter_voxel_distance/2.0;
	initial_offset[1] = BB[0][1] + inter_voxel_distance/2.0;
	initial_offset[2] = BB[0][2] + inter_voxel_distance/2.0;//       ½    ؿ    ĵ      

	for (int i = 0; i < 3; i++){
		number_voxels_per_dim[i] = floor((BB[1][i] - BB[0][i])/inter_voxel_distance);
		cout << "Number voxels per " << number_voxels_per_dim[i] << endl;
	}//    ÿ     м      ؿ 
	//number_voxels_per_dim[1] = floor((BB.max(1) - BB.min(1))/inter_voxel_distance);
	//number_voxels_per_dim[2] = floor((BB.max(2) - BB.min(2))/inter_voxel_distance);

	//	number_voxels_per_dim[0] = floor((BB.max(0) - BB.min(0))/inter_voxel_distance) - 1;
	//	number_voxels_per_dim[1] = floor((BB.max(1) - BB.min(1))/inter_voxel_distance);
	//	number_voxels_per_dim[2] = floor((BB.max(2) - BB.min(2))/inter_voxel_distance);

	// configuration grid is only for the voxels.
	number_voxels_grid = number_voxels_per_dim[0]*number_voxels_per_dim[1]*number_voxels_per_dim[2];//һ   м      ؿ 



	// later - make faster, make an array for ALL of these things ...
	configuration_grid.resize(number_voxels_grid, true);

	//number_of_positive_terms.resize(number_voxels_grid, 0);
	set_by_first_d.resize(number_voxels_grid, false);

	// can make all of these arrays ....
	//	start_index_for_terms_arr = new int_type_t[number_voxels_grid + 1];
	term_increment = pow(2, 24);

	number_cameras = cameras.size();

	//	rows = cameras[0]->im_original.rows/downsample;
	//	cols = cameras[0]->im_original.cols/downsample;

	rows = cameras[0]->imcv.rows/downsample;
	cols = cameras[0]->imcv.cols/downsample;

	number_divisions = omp_get_max_threads();

	number_pixels = number_cameras*rows*cols;

	current_pixel_labeling.resize(number_pixels, true);
	// this can be an array instead of a vector for speed .... especially b/c number of pixels can be large ...
	//is_in_undistorted_image.resize(number_pixels, true);
	is_in_undistorted_image = RS->is_in_undistorted_image;
	number_zeros_grid.resize(number_pixels, 0);//Ӧ ñ ʾÿ     ض Ӧ      ǩΪ0      

	for (int_type_t i = 0; i < number_divisions; i++){
		multi_dim_voxels_to_process.push_back(new bool[number_voxels_grid/number_divisions  + 1]);
		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
			multi_dim_voxels_to_process[i][j] = false;
		}
	}

	cout << "Before .... " << endl;

		tile_number.resize(number_divisions, 0);
		for (int_type_t i = 0; i < number_divisions; i++){
			tile_number[i] = new int16_t[number_voxels_grid/number_divisions + 2];
		}


		//tile_number = new int16_t[number_voxels_grid];
		for (int_type_t i = 0; i < number_divisions; i++){
			for (int_type_t j = 0; j < number_voxels_grid/number_divisions + 2; j++){
				tile_number[i][j] = -1;
			}
		}
		cout << "After " << endl;

	voxel_ordering.resize(number_voxels_grid/number_divisions  + 1, 0);
	for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
		voxel_ordering[j] = j;
	}

	//cout << "Midway through ... " << endl;
	//char ch; cin >> ch;


	Matrix P(3, 4);
	Matrix Pinv;

	vector<vector<double> > Pvv(3, vector<double>(4, 0));
	vector<vector<double> > Pivv(3, vector<double>(3, 0));
	vector<double> Cvv(4, 1);
	vector<double> Dvv(3);
	vector<double> lambDvv(3);
	vector<double> Xvv(3);
	vector<double> xvv(3);
	vector<double> Rvv(3);
	vector<double> plane_equation(4);//ʲô  ˼    


	// do projection matrices....
	for (int cam = 0; cam < number_cameras; cam++){

		for (int i = 0; i < 3; i++){
			P.Row(i+1) << cameras[cam]->P.at<double>(i, 0) << cameras[cam]->P.at<double>(i, 1) << cameras[cam]->P.at<double>(i, 2) << cameras[cam]->P.at<double>(i, 3);
		}

		// only first two rows ....
		for (int i = 0; i < 2; i++){
			for (int j = 0; j < 4; j++){
				P(i+1, j+1) /= double(downsample);
			}
		}

		//cout << " P " << endl << P << endl;
		//cout << "Original representation " << endl;
		//MatrixPrint(cameras[cam]->P);
		Pinv = P.t() * (P*P.t()).i();

		Pinv = Pinv /Pinv (4, 3);

		//cout << "Pinv " << endl << Pinv << endl;


		// do loop to read in ....
		for (int i = 0; i < 3; i++){
			Cvv[i] = (cameras[cam]->C.at<double>(i, 0));
			for (int c = 0; c < 3; c++){
				Pivv[i][c] = Pinv(i + 1, c + 1);
			}

			for (int c = 0; c < 4; c++){
				Pvv[i][c] = P(i + 1, c + 1);
			}
		}

		//cout << "New mat projection matrices " << endl << P << endl;


		//		cout << "Projection matrices ... " << endl;
		//		for (int i = 0; i < 3; i++){
		//
		//			for (int c = 0; c < 4; c++){
		//				cout << Pvv[i][c] << " "; //= P(i + 1, c + 1);
		//			}
		//			cout << endl;
		//		}

		projection_matrices.push_back(Pvv);
		inverse_projection_matrices.push_back(Pivv);
		camera_centers.push_back(Cvv);

		plane_equation[0] = cameras[cam]->P.at<double>(2, 0);
		plane_equation[1] = cameras[cam]->P.at<double>(2, 1);
		plane_equation[2] = cameras[cam]->P.at<double>(2, 2);
		plane_equation[3] = cameras[cam]->P.at<double>(2, 3);

		NormalizePlane(plane_equation);
		plane_equations.push_back(plane_equation);


		//		plane_equation[0] = cameras[cam]->PrincipalPlane.a();
		//		plane_equation[1] = cameras[cam]->PrincipalPlane.b();
		//		plane_equation[2] = cameras[cam]->PrincipalPlane.c();
		//		plane_equation[3] = cameras[cam]->PrincipalPlane.d();
		//		NormalizePlane(plane_equation);
		//		plane_equations.push_back(plane_equation);

	}


	int index;
	// for now, assume that all images are in color.  Technically, we can deallocate the undistortion map after the next loop ....
	bool color = true;
	bool to_test;

	//cout << "rows, cols " << rows << ", " <<


	//	cv::Mat undistort_map(rows, cols, CV_8UC1, cv::Scalar(0,0,0));
	//	// can do this faster with a ptr .... not relevant, b/c so fast already?
	//	for (int_type_t cam = 0; cam < number_cameras; cam++){
	//
	//
	//		cout << "original map " << cameras[cam]->undistorted_map_original.rows << " " << cameras[cam]->undistorted_map_original.cols << endl;
	//		cout << "Resize map " << undistort_map.rows << ", " << undistort_map.cols << endl;
	//
	//		color = cameras[cam]->undistorted_map_original.channels() == 3;
	//		uchar* p = 0;
	//		// resize for each time ...
	//		if (downsample > 1){
	//			cv::resize(cameras[cam]->undistorted_map_original, undistort_map, undistort_map.size());
	//			p = &undistort_map.data[0];
	//		}	else {
	//			p = &cameras[cam]->undistorted_map_original.data[0];
	//		}
	//
	//		if (color){
	//			for (int r = 0; r < rows; r++){
	//				for (int c = 0; c < cols; c++){
	//
	//					to_test = *p > 150;
	//
	//					if (!to_test){
	//						// we don't want to look this up for the majority of pixels that are "inside"
	//						// can also use an iterator -- speed?
	//						is_in_undistorted_image[cam*(rows*cols) + r*cols + c] = to_test;
	//					}
	//
	//					// increment by three ....
	//					p++;
	//					p++;
	//					p++;
	//				}
	//			}
	//
	//
	//		}	else {
	//			//cout << "Implement for grayscale! " << endl;
	//			for (int r = 0; r < rows; r++){
	//				for (int c = 0; c < cols; c++){
	//
	//					to_test = *p > 150;
	//
	//					if (!to_test){
	//						// we don't want to look this up for the majority of pixels that are "inside"
	//						// can also use an iterator -- speed?
	//						is_in_undistorted_image[cam*(rows*cols) + r*cols + c] = to_test;
	//					}
	//
	//					// increment by one b/c grayscale ....
	//					p++;
	//				}
	//			}
	//
	//		}
	//	}

	//usleep(100000);


}

void ReconstructionStructure::CreateGridUsingAncestorWithProjection(ReconstructionStructure* RS, vector< vector<double> >& boundingbox, double division, int downsample){

	//	vector<bool> configuration_grid;
	//		vector<bool> part_of_visual_hull;
	//
	//		vector<double> initial_offset;

	inter_voxel_distance = division;
	BB = boundingbox;
	initial_offset.resize(3, 0);
	number_voxels_per_dim.resize(3, 0);

	initial_offset[0] = BB[0][0] + inter_voxel_distance/2.0;
	initial_offset[1] = BB[0][1] + inter_voxel_distance/2.0;
	initial_offset[2] = BB[0][2] + inter_voxel_distance/2.0;

	for (int i = 0; i < 3; i++){
		number_voxels_per_dim[i] = floor((BB[1][i] - BB[0][i])/inter_voxel_distance);
		cout << "Number voxels per " << number_voxels_per_dim[i] << endl;
	}
	//number_voxels_per_dim[1] = floor((BB.max(1) - BB.min(1))/inter_voxel_distance);
	//number_voxels_per_dim[2] = floor((BB.max(2) - BB.min(2))/inter_voxel_distance);

	//	number_voxels_per_dim[0] = floor((BB.max(0) - BB.min(0))/inter_voxel_distance) - 1;
	//	number_voxels_per_dim[1] = floor((BB.max(1) - BB.min(1))/inter_voxel_distance);
	//	number_voxels_per_dim[2] = floor((BB.max(2) - BB.min(2))/inter_voxel_distance);

	// configuration grid is only for the voxels.
	number_voxels_grid = number_voxels_per_dim[0]*number_voxels_per_dim[1]*number_voxels_per_dim[2];



	// later - make faster, make an array for ALL of these things ...
	configuration_grid.resize(number_voxels_grid, true);

	//number_of_positive_terms.resize(number_voxels_grid, 0);
	set_by_first_d.resize(number_voxels_grid, false);

	// can make all of these arrays ....
	//	start_index_for_terms_arr = new int_type_t[number_voxels_grid + 1];
	term_increment = pow(2, 24);

	number_cameras = RS->number_cameras;

	//	rows = cameras[0]->im_original.rows/downsample;
	//	cols = cameras[0]->im_original.cols/downsample;

	rows = RS->rows;
	cols = RS->cols;

	number_divisions = omp_get_max_threads();

	number_pixels = number_cameras*rows*cols;

	current_pixel_labeling.resize(number_pixels, true);
	// this can be an array instead of a vector for speed .... especially b/c number of pixels can be large ...
	//is_in_undistorted_image.resize(number_pixels, true);
	is_in_undistorted_image = RS->is_in_undistorted_image;
	number_zeros_grid.resize(number_pixels, 0);

	for (int_type_t i = 0; i < number_divisions; i++){
		multi_dim_voxels_to_process.push_back(new bool[number_voxels_grid/number_divisions  + 1]);
		for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
			multi_dim_voxels_to_process[i][j] = false;
		}
	}

	cout << "Before .... " << endl;

		tile_number.resize(number_divisions, 0);
		for (int_type_t i = 0; i < number_divisions; i++){
			tile_number[i] = new int16_t[number_voxels_grid/number_divisions + 2];
		}


		//tile_number = new int16_t[number_voxels_grid];
		for (int_type_t i = 0; i < number_divisions; i++){
			for (int_type_t j = 0; j < number_voxels_grid/number_divisions + 2; j++){
				tile_number[i][j] = -1;
			}
		}
		cout << "After " << endl;

	voxel_ordering.resize(number_voxels_grid/number_divisions  + 1, 0);
	for (int_type_t j = 0; j < number_voxels_grid/number_divisions  + 1; j++){
		voxel_ordering[j] = j;
	}

	//cout << "Midway through ... " << endl;
	//char ch; cin >> ch;


	Matrix P(3, 4);
	Matrix Pinv;

	vector<vector<double> > Pvv(3, vector<double>(4, 0));
	vector<vector<double> > Pivv(3, vector<double>(3, 0));
	vector<double> Cvv(4, 1);
	vector<double> Dvv(3);
	vector<double> lambDvv(3);
	vector<double> Xvv(3);
	vector<double> xvv(3);
	vector<double> Rvv(3);
	vector<double> plane_equation(4);

	//projection_matrices = RS->projection_matrices;
	//inverse_projection_matrices = RS->inverse_projection_matrices;
	//camera_centers = RS->camera_centers;
	//plane_equations = RS->plane_equations;


	for (int cam = 0; cam < number_cameras; cam++){
		// projection matrices ....
		//out << "cam " << cam << endl;
		for (int i = 0; i < 3; i++){
			Cvv[i] = RS->camera_centers[cam][i];
			for (int c = 0; c < 4; c++){
				Pvv[i][c] = RS->projection_matrices[cam][i][c];
				//out << RS_best->projection_matrices[cam][i][c] << " ";
			}

			for (int c = 0; c < 3; c++){
				Pivv[i][c] = RS->inverse_projection_matrices[cam][i][c];
			}
			//out << endl;

		}

		for (int i = 0; i < 4; i++){
			plane_equation[i] = RS->plane_equations[cam][i];
		}


		projection_matrices.push_back(Pvv);
		inverse_projection_matrices.push_back(Pivv);
		camera_centers.push_back(Cvv);
		plane_equations.push_back(plane_equation);

		//out << endl;
	}


//	// do projection matrices....
//	for (int cam = 0; cam < number_cameras; cam++){
//		projection_matrices.push_back(RS->projection_matrices[cam]);
//		inverse_projection_matrices.push_back(RS->inverse_projection_matrices[cam]);
//		camera_centers.push_back(RS->camera_centers[cam]);
//		plane_equations.push_back(RS->plane_equations[cam]);
//	}


	int index;
	// for now, assume that all images are in color.  Technically, we can deallocate the undistortion map after the next loop ....
	bool color = true;
	bool to_test;






}

void ReconstructionStructure::ResetProjectionMatricesAndNumberZeros(vector<ImageData*>& cameras, int downsample){

	char ch;
	cout << "Downsample " << downsample << endl;
	cout << "Before number zeros grid " << endl; //cin >> ch;
	for (int_type_t i = 0; i < number_voxels_grid; i++){
		number_zeros_grid[i] = 0;
	}
	cout << "After number zeros grid " << endl; //cin >> ch;

	Matrix P(3, 4);
	Matrix Pinv;

	vector<vector<double> > Pvv(3, vector<double>(4, 0));
	vector<vector<double> > Pivv(3, vector<double>(3, 0));
	vector<double> Cvv(4, 1);
	vector<double> Dvv(3);
	vector<double> lambDvv(3);
	vector<double> Xvv(3);
	vector<double> xvv(3);
	vector<double> Rvv(3);
	vector<double> plane_equation(4);


	// do projection matrices....
	for (int cam = 0; cam < number_cameras; cam++){
		cout << "cam  " << cam << endl; //cin >> ch;

		for (int i = 0; i < 3; i++){
			P.Row(i+1) << cameras[cam]->P.at<double>(i, 0) << cameras[cam]->P.at<double>(i, 1) << cameras[cam]->P.at<double>(i, 2) << cameras[cam]->P.at<double>(i, 3);
		}

		// only first two rows ....
		for (int i = 0; i < 2; i++){
			for (int j = 0; j < 4; j++){
				P(i+1, j+1) /= double(downsample);
			}
		}

		//cout << "Original representation " << endl;
		//MatrixPrint(cameras[cam]->P);
		Pinv = P.t() * (P*P.t()).i();

		Pinv = Pinv /Pinv (4, 3);

		//cout << "Pinv " << endl << Pinv << endl;


		// do loop to read in ....
		for (int i = 0; i < 3; i++){
			Cvv[i] = (cameras[cam]->C.at<double>(i, 0));
			for (int c = 0; c < 3; c++){
				Pivv[i][c] = Pinv(i + 1, c + 1);
			}

			for (int c = 0; c < 4; c++){
				Pvv[i][c] = P(i + 1, c + 1);
			}
		}


		projection_matrices[cam] = (Pvv);
		inverse_projection_matrices[cam] = (Pivv);
		camera_centers[cam] = (Cvv);

		//		plane_equation[0] = cameras[cam]->PrincipalPlane.a();
		//		plane_equation[1] = cameras[cam]->PrincipalPlane.b();
		//		plane_equation[2] = cameras[cam]->PrincipalPlane.c();
		//		plane_equation[3] = cameras[cam]->PrincipalPlane.d();
		//		NormalizePlane(plane_equation);
		//		plane_equations.push_back(plane_equation);

	}

}





void ReconstructionStructure::MarkToTestAsFalseFirstDsOnly(){
	int_type_t voxel_index = 0;

	for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){
		for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
			for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){

				if (set_by_first_d[voxel_index] == true && configuration_grid[voxel_index] == false){
					multi_dim_voxels_to_process[voxel_index % number_divisions][voxel_index/number_divisions] = false;
					//cout << "Set some ... " << endl;

				}	else {
					//configuration_grid[voxel_index] = true;

					//					if (configuration_grid[voxel_index] == true){
					//						multi_dim_voxels_to_process[voxel_index % number_divisions][voxel_index/number_divisions] = false;
					//					}
				}

//				if (multi_dim_voxels_to_process[voxel_index % number_divisions][voxel_index/number_divisions] == true){
//					configuration_grid[voxel_index] = false;
//				}

			}
		}
	}
	//char ch; cin >> ch;
}

int_type_t ReconstructionStructure::MarkShellAsToProcess(){

	int_type_t voxel_index = 0;
	int_type_t neighbor_index;
	bool has_empty_neighbor;

	int_type_t voxels_to_process = 0;

	for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

		if (x_index % 1000 == 0){
			cout << "Marking shell... " << x_index << " out of " << number_voxels_per_dim[0] << endl;
		}
		for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
			for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){

				if (set_by_first_d[voxel_index] == true && configuration_grid[voxel_index] == false){
					has_empty_neighbor = false;
					// search for an empty neighbor ...
					// do 4-conntcted neighbors for now ...
					if (x_index > 0){
						neighbor_index =  (x_index - 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;
						if (set_by_first_d[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && x_index < number_voxels_per_dim[0] - 1){
						neighbor_index =  (x_index + 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;
						if (set_by_first_d[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && y_index > 0){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y_index - 1)*number_voxels_per_dim[2] + z_index;
						if (set_by_first_d[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && y_index < number_voxels_per_dim[1] - 1){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y_index + 1)*number_voxels_per_dim[2] + z_index;
						if (set_by_first_d[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && z_index > 0){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index - 1;
						if (set_by_first_d[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && z_index < number_voxels_per_dim[2] - 1){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index + 1;
						if (set_by_first_d[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (has_empty_neighbor){
						multi_dim_voxels_to_process[voxel_index % number_divisions][voxel_index/number_divisions] = true;
						voxels_to_process++;
					}
				}
			}
		}
	}

	return voxels_to_process;
}


int_type_t ReconstructionStructure::PropogateLabelsAccordingToChild(ReconstructionStructure& RS_child, int division, int inner_distance, int distance){

	//bool verbose = division == 1;
	bool verbose = false;
	// there's somethign wonky going on here with the neighbor indices .... check it out.
	// assume that parent is a div 2 ...  so instead of the 2* and 2* + 1, we do 2* and 2*  - 1

	//voxel_index =  x_index*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;

	// bad use of omp?
	//#pragma omp parallel for
	char ch;
	cout << "Distance for propogation is " << distance << endl;
	cout << "Inner distance is " << inner_distance << endl;
	cout << "Line 506 " << endl;
	if (verbose){
		cin >> ch;
	}
	//cin >> ch;
	if (distance > 0){
		for (int_type_t i = 0; i < RS_child.number_voxels_grid; i++){
			RS_child.set_by_first_d[i] = true;//Ϊʲô  ô      
		}
	}
	cout << "Number of voxels in the grid " << number_voxels_grid << endl;
	cout << "Line 613 " << endl;
	if (verbose){
		cin >> ch;
	}
	//cin >> ch;

	// ?? too big?
	bool* live_voxels = new bool[number_voxels_grid];
	//bool live_voxels[number_voxels_grid]; // gets destrpyed at the end of the function
	for (int_type_t i = 0; i < number_voxels_grid; i++){
		live_voxels[i] = false; //ʲô  ˼    
	}
	cout << "Line 620 " << endl;
	//cin >> ch;
	vector<int_type_t> testing_list;


	int_type_t voxels_to_test = 0;

	int_type_t voxel_index = 0;
	int_type_t small_index, neighbor_index;
	bool first_d_value;
	bool config_value;
	bool has_empty_neighbor, has_occupied_neighbor;
	int_type_t x0, y0, z0, startx0, starty0, startz0;
	int_type_t kid_distance = division*distance;

	int number_falses = 0;

	// TODO later parallel if appropriate ...
	for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

		if (x_index % 100 == 0){
			cout << "Propogating ... " << x_index << " out of " << number_voxels_per_dim[0] << endl;
		}
		for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
			for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){
				//
				//first_d_value = set_by_first_d[voxel_index];
				config_value = configuration_grid[voxel_index];


				if (config_value == false)
				{
					number_falses++;
					for (x0 = division*x_index; x0 < division*x_index + division; x0++){
						for (y0 = division*y_index; y0 < division*y_index + division; y0++ ){
							for (z0 = division*z_index; z0 < division*z_index + division; z0++ ){//    division==2,  ʾһ     ؿ ָ ɰ˸            ؿ ı߳ 
								small_index =  x0*RS_child.number_voxels_per_dim[1]*RS_child.number_voxels_per_dim[2] + y0*RS_child.number_voxels_per_dim[2] + z0;
								RS_child.configuration_grid[small_index] = false;
								//RS_child.set_by_first_d[small_index] = false;
							}
						}
					}//       ر ǩ          

					if (verbose){
						cout << "voxel index " << voxel_index << " " << small_index << endl;
					}


					has_empty_neighbor = false;
					// search for an empty neighbor ...
					// do 4-conntcted neighbors for now ...
					if (x_index > 0){
						neighbor_index =  (x_index - 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}
					//  Щ ж      ֱ   ʲô  ˼    ΪʲôҪ  ô    
					if (!has_empty_neighbor && x_index < number_voxels_per_dim[0] - 1){
						neighbor_index =  (x_index + 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && y_index > 0){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y_index - 1)*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && y_index < number_voxels_per_dim[1] - 1){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y_index + 1)*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && z_index > 0){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index - 1;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && z_index < number_voxels_per_dim[2] - 1){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index + 1;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}//ΪʲôҪ    ô   ж         


					//has_empty_neighbor = false;
					if (has_empty_neighbor){

						live_voxels[voxel_index] = true;//live_voxels    ʲô    
						testing_list.push_back(voxel_index);//testing_list    ʲô  ˼    

						//						// make the kids have a fasl set by first d. ...
						//						//cout << "limits " << max(0, division*(x_index - distance)) << ", " <<  min(RS_child.number_voxels_per_dim[0], division*(x_index + distance)) << endl;
						//						for (x0 = max(0, division*(x_index - distance)); x0 < min(RS_child.number_voxels_per_dim[0], division*(x_index + distance + 1)); x0++){
						//							for (y0 = max(0, division*(y_index - distance)); y0 < min(RS_child.number_voxels_per_dim[1], division*(y_index + distance + 1)); y0++){
						//								for (z0 = max(0, division*(z_index - distance)); z0 < min(RS_child.number_voxels_per_dim[2], division*(z_index + distance + 1)); z0++){
						//									small_index =  x0*RS_child.number_voxels_per_dim[1]*RS_child.number_voxels_per_dim[2] + y0*RS_child.number_voxels_per_dim[2] + z0;
						//									RS_child.set_by_first_d[small_index] = false;
						//								}
						//							}
						//						}
					}
				}
			}
		}
	}

	cout << "Line 730, testing list size .... " << testing_list.size() <<  endl;
	if (verbose){
		cin >> ch;
	}
	//cin >> ch;
	// expand the live voxels ....
	// for speed ... change to a C-style array
	cout << "Testing list " << testing_list.size() << endl;
	for (int i = 0; i < distance; i++){
		cout << "Distance " << distance << endl;
		vector<int_type_t> next_test;
		for (int_type_t j = 0, jn = testing_list.size(); j < jn; j++){

			voxel_index = testing_list[j];
			ReturnXYZIndicesFromIndex(voxel_index, x0, y0, z0);

			// test 6 connected directions ...
			if (x0 > 0){
				neighbor_index =  (x0 - 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y0*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (x0 < number_voxels_per_dim[0] - 1){
				neighbor_index =  (x0 + 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y0*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (y0 > 0){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0 - 1)*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (y0 < number_voxels_per_dim[1] - 1){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0 + 1)*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if ( z0 > 0){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0)*number_voxels_per_dim[2] + z0 - 1;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (z0 < number_voxels_per_dim[2] - 1){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0)*number_voxels_per_dim[2] + z0 + 1;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

		}

		next_test.swap(testing_list);
		cout << "Next round expand ... " << testing_list.size() << endl;
	}

	voxel_index = 0;
	for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

		if (x_index % 100 == 0){
			cout << "Filling in ... " << x_index << " out of " << number_voxels_per_dim[0] << endl;
		}
		for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
			for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){
				//

				if (live_voxels[voxel_index] == true)
				{

					for (x0 = division*x_index; x0 < division*x_index + division; x0++){
						for (y0 = division*y_index; y0 < division*y_index + division; y0++ ){
							for (z0 = division*z_index; z0 < division*z_index + division; z0++ ){
								small_index =  x0*RS_child.number_voxels_per_dim[1]*RS_child.number_voxels_per_dim[2] + y0*RS_child.number_voxels_per_dim[2] + z0;
								RS_child.set_by_first_d[small_index] = false;
							}
						}
					}
				}
			}
		}
	}


	// also make the "shell" part of the ones to process
	voxels_to_test = RS_child.MarkShellAsToProcess();

	cout << "Shell marked .. " << voxels_to_test << endl;

	//char ch;
	if (verbose){
		cin >> ch;
	}

	for (int_type_t i = 0; i < RS_child.number_voxels_grid; i++){
		if (RS_child.set_by_first_d[i] == false){
			voxels_to_test++;
			RS_child.multi_dim_voxels_to_process[i % RS_child.number_divisions][i / RS_child.number_divisions] = true;

			//RS_child.configuration_grid[i] = true;
		}
	}

	// so all of the false first ds are also in multi dim voxels to process

	// walk trhough and find border ones ....

	cout << "Voxels to test ... " << voxels_to_test << endl;
	cout << "Number falses " << number_falses << endl;
	if (verbose){
		cin >> ch;
	}
	//cin >> ch;
	delete [] live_voxels;

	return voxels_to_test;
}


int_type_t ReconstructionStructure::PropogateLabelsAndMark(int distance){

	// there's somethign wonky going on here with the neighbor indices .... check it out.
	// assume that parent is a div 2 ...  so instead of the 2* and 2* + 1, we do 2* and 2*  - 1

	//voxel_index =  x_index*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;

	// bad use of omp?
	//#pragma omp parallel for
	if (distance > 0){
		for (int_type_t i = 0; i < number_voxels_grid; i++){
			set_by_first_d[i] = true;
		}
	}

	bool* live_voxels = new bool[number_voxels_grid];
	//bool live_voxels[number_voxels_grid]; // gets destrpyed at the end of the function
	for (int_type_t i = 0; i < number_voxels_grid; i++){
		live_voxels[i] = false;
	}
	vector<int_type_t> testing_list;


	int_type_t voxels_to_test = 0;

	int_type_t voxel_index = 0;
	int_type_t small_index, neighbor_index;
	bool first_d_value;
	bool config_value;
	bool has_empty_neighbor, has_occupied_neighbor;
	int_type_t x0, y0, z0, startx0, starty0, startz0;


	// TODO later parallel if appropriate ...
	for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

		if (x_index % 100 == 0){
			cout << "Propogating ... " << x_index << " out of " << number_voxels_per_dim[0] << endl;
		}
		for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
			for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){
				//
				//first_d_value = set_by_first_d[voxel_index];
				config_value = configuration_grid[voxel_index];


				if (config_value == false)
				{

					has_empty_neighbor = false;
					// search for an empty neighbor ...
					// do 4-conntcted neighbors for now ...
					if (x_index > 0){
						neighbor_index =  (x_index - 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && x_index < number_voxels_per_dim[0] - 1){
						neighbor_index =  (x_index + 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && y_index > 0){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y_index - 1)*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && y_index < number_voxels_per_dim[1] - 1){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y_index + 1)*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && z_index > 0){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index - 1;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && z_index < number_voxels_per_dim[2] - 1){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index + 1;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}


					//has_empty_neighbor = false;
					if (has_empty_neighbor){

						live_voxels[voxel_index] = true;
						testing_list.push_back(voxel_index);
					}
				}
			}
		}
	}

	// expand the live voxels ....
	// for speed ... change to a C-style array
	cout << "Testing list " << testing_list.size() << endl;
	for (int i = 0; i < distance; i++){
		cout << "Distance " << distance << endl;
		vector<int_type_t> next_test;
		for (int_type_t j = 0, jn = testing_list.size(); j < jn; j++){

			voxel_index = testing_list[j];
			ReturnXYZIndicesFromIndex(voxel_index, x0, y0, z0);

			// test 6 connected directions ...
			if (x0 > 0){
				neighbor_index =  (x0 - 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y0*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (x0 < number_voxels_per_dim[0] - 1){
				neighbor_index =  (x0 + 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y0*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (y0 > 0){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0 - 1)*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (y0 < number_voxels_per_dim[1] - 1){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0 + 1)*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if ( z0 > 0){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0)*number_voxels_per_dim[2] + z0 - 1;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (z0 < number_voxels_per_dim[2] - 1){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0)*number_voxels_per_dim[2] + z0 + 1;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

		}

		next_test.swap(testing_list);
		cout << "Next round expand ... " << testing_list.size() << endl;
	}

	voxel_index = 0;
	for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

		if (x_index % 100 == 0){
			cout << "Filling in ... " << x_index << " out of " << number_voxels_per_dim[0] << endl;
		}
		for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
			for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){
				//

				if (live_voxels[voxel_index] == true)
				{
					set_by_first_d[voxel_index] = false;
				}
			}
		}
	}


	// also make the "shell" part of the ones to process
	voxels_to_test = MarkShellAsToProcess();

	cout << "Shell marked .. " << voxels_to_test << endl;

	char ch;

	for (int_type_t i = 0; i < number_voxels_grid; i++){


		if (set_by_first_d[i] == false){
			voxels_to_test++;
			multi_dim_voxels_to_process[i % number_divisions][i / number_divisions] = true;
		}
	}

	// walk trhough and find border ones ....

	cout << "Voxels to test ... " << voxels_to_test << endl;
	//cin >> ch;

	delete [] live_voxels;

	return voxels_to_test;
}


int_type_t ReconstructionStructure::SameLevelExpansion(int distance){

	bool* live_voxels = new bool[number_voxels_grid];

	for (int_type_t i = 0; i < number_voxels_grid; i++){
		live_voxels[i] = false;
	}
	vector<int_type_t> testing_list;


	int_type_t voxels_to_test = 0;

	int_type_t voxel_index = 0;
	int_type_t small_index, neighbor_index;
	bool first_d_value;
	bool config_value;
	bool has_empty_neighbor, has_occupied_neighbor;
	int_type_t x0, y0, z0, startx0, starty0, startz0;


	// TODO later parallel if appropriate ...
	for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

		if (x_index % 100 == 0){
			cout << "Propogating ... " << x_index << " out of " << number_voxels_per_dim[0] << endl;
		}
		for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
			for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){
				//
				first_d_value = set_by_first_d[voxel_index];
				config_value = configuration_grid[voxel_index];


				// in other words, this one was live for previous iterations .... want to mark all of the voxels within a certain distance from this voxel.
				// we're only doing the shell -- outer layer -- of the reconstruction to save on time and overlap.
				if (config_value == false && first_d_value == false)
				{

					has_empty_neighbor = false;
					// search for an empty neighbor ...
					// do 4-conntcted neighbors for now ...
					if (x_index > 0){
						neighbor_index =  (x_index - 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && x_index < number_voxels_per_dim[0] - 1){
						neighbor_index =  (x_index + 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && y_index > 0){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y_index - 1)*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && y_index < number_voxels_per_dim[1] - 1){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y_index + 1)*number_voxels_per_dim[2] + z_index;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && z_index > 0){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index - 1;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}

					if (!has_empty_neighbor && z_index < number_voxels_per_dim[2] - 1){
						neighbor_index =  (x_index)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index + 1;
						if (configuration_grid[neighbor_index] == true){
							has_empty_neighbor = true;
						}
					}


					//has_empty_neighbor = false;
					if (has_empty_neighbor){

						live_voxels[voxel_index] = true;
						testing_list.push_back(voxel_index);
					}
				}
			}
		}
	}

	// expand the live voxels ....
	// for speed ... change to a C-style array
	cout << "Testing list " << testing_list.size() << endl;
	for (int i = 0; i < distance; i++){
		cout << "Distance " << distance << endl;
		vector<int_type_t> next_test;
		for (int_type_t j = 0, jn = testing_list.size(); j < jn; j++){

			voxel_index = testing_list[j];
			ReturnXYZIndicesFromIndex(voxel_index, x0, y0, z0);

			// test 6 connected directions ...
			if (x0 > 0){
				neighbor_index =  (x0 - 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y0*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (x0 < number_voxels_per_dim[0] - 1){
				neighbor_index =  (x0 + 1)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y0*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (y0 > 0){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0 - 1)*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (y0 < number_voxels_per_dim[1] - 1){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0 + 1)*number_voxels_per_dim[2] + z0;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if ( z0 > 0){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0)*number_voxels_per_dim[2] + z0 - 1;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

			if (z0 < number_voxels_per_dim[2] - 1){
				neighbor_index =  (x0)*number_voxels_per_dim[1]*number_voxels_per_dim[2] + (y0)*number_voxels_per_dim[2] + z0 + 1;
				if (live_voxels[neighbor_index] == false){
					live_voxels[neighbor_index] = true;
					next_test.push_back(neighbor_index);
				}
			}

		}

		next_test.swap(testing_list);
		cout << "Next round expand ... " << testing_list.size() << endl;
	}

	voxel_index = 0;
	voxels_to_test = 0;
	for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

		if (x_index % 100 == 0){
			cout << "Filling in ... " << x_index << " out of " << number_voxels_per_dim[0] << endl;
		}
		for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
			for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){
				//

				if (live_voxels[voxel_index] == true && set_by_first_d[voxel_index] == true)
				{
					set_by_first_d[voxel_index] = false;
					voxels_to_test++;
					multi_dim_voxels_to_process[voxel_index % number_divisions][voxel_index / number_divisions] = true;
				}
			}
		}
	}


	// also make the "shell" part of the ones to process


	char ch;

	//	for (int_type_t i = 0; i < number_voxels_grid; i++){
	//
	//
	//		if (set_by_first_d[i] == false){
	//			voxels_to_test++;
	//			multi_dim_voxels_to_process[i % number_divisions][i / number_divisions] = true;
	//		}
	//	}

	// walk trhough and find border ones ....

	cout << "Voxels to test ... " << voxels_to_test << endl;
	//cin >> ch;

	delete [] live_voxels;

	return voxels_to_test;
}

//int_type_t ReconstructionStructure::PropogateLabelsAccordingToChild(ReconstructionStructure& RS_child, int distance){
//
//	// there's somethign wonky going on here with the neighbor indices .... check it out.
//	// assume that parent is a div 2 ...  so instead of the 2* and 2* + 1, we do 2* and 2*  - 1
//
//	//voxel_index =  x_index*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;
//
//	// bad use of omp?
//	//#pragma omp parallel for
//	for (int_type_t i = 0; i < RS_child.number_voxels_grid; i++){
//		RS_child.set_by_first_d[i] = true;
//	}
//
//	int_type_t voxels_to_test = 0;
//
//	int_type_t voxel_index = 0;
//	int_type_t small_index, neighbor_index;
//	bool first_d_value;
//	bool config_value;
//	bool has_empty_neighbor, has_occupied_neighbor;
//	int_type_t x0, y0, z0, startx0, starty0, startz0;
//	int_type_t kid_distance = 2*distance;
//
//	for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){
//
//		if (x_index % 100 == 0){
//			cout << "Propogating ... " << x_index << " out of " << number_voxels_per_dim[0] << endl;
//		}
//		for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
//			for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){
//				//
//				first_d_value = set_by_first_d[voxel_index];
//				config_value = configuration_grid[voxel_index];
//
//				// can we do has an occupied neighbor from here?
//				// copy over the occupied voxels ....
//				//if (first_d_value == true || config_value == false){
//				has_occupied_neighbor = false;
////				if (config_value == true){
////					// go through neighbors
////
////
////
////				}
//
//
//
//				//if (config_value == false || first_d_value == false)
//				{
//					//										x0 = 2*(x_index);
//					//										y0 = 2*(y_index);
//					//										z0 = 2*(z_index);
//					//					//					//x0 = x_index; y0 = y_index; z0 = z_index;
//					//										small_index =  x0*RS_child.number_voxels_per_dim[1]*RS_child.number_voxels_per_dim[2] + y0*RS_child.number_voxels_per_dim[2] + z0;
//					//										RS_child.configuration_grid[small_index] = config_value;
//					//										RS_child.set_by_first_d[small_index] = first_d_value;
//
//					//first, determine if this voxel is on the edge -- within distance of an empty voxel
//					has_empty_neighbor = false;
//					if (x_index < distance){
//						startx0 = 0;
//					}	else {
//						startx0 = x_index - distance;
//					}
//
//					if (y_index < distance){
//						starty0 = 0;
//					}	else {
//						starty0 = y_index - distance;
//					}
//
//					if (z_index < distance){
//						startz0 = 0;
//					}	else {
//						startz0 = z_index - distance;
//					}
//
//
//
//					for (x0 = startx0 ; x0 < min(x_index + distance + 1, number_voxels_per_dim[0]) && !has_empty_neighbor; x0++){
//						for (y0 = starty0 ; y0 < min(y_index + distance + 1, number_voxels_per_dim[1]) && !has_empty_neighbor; y0++ ){
//							for (z0 = startz0 ; z0 < min(z_index + distance + 1, number_voxels_per_dim[2]) && !has_empty_neighbor; z0++ ){
//								neighbor_index = x0*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y0*number_voxels_per_dim[2] + z0;
//								if (configuration_grid[neighbor_index] == true){
//									has_empty_neighbor = true;
//								}
//							}
//						}
//					}
//
//
//					for (x0 = 2*x_index; x0 < 2*x_index + 2; x0++){
//						for (y0 = 2*y_index; y0 < 2*y_index + 2; y0++ ){
//							for (z0 = 2*z_index; z0 < 2*z_index + 2; z0++ ){
//								small_index =  x0*RS_child.number_voxels_per_dim[1]*RS_child.number_voxels_per_dim[2] + y0*RS_child.number_voxels_per_dim[2] + z0;
//								RS_child.configuration_grid[small_index] = config_value;
//								//RS_child.configuration_grid[small_index] = has_empty_neighbor;
//								RS_child.set_by_first_d[small_index] = false;
//							}
//						}
//					}
//					//voxels_to_test += 8;
////
////					if (has_empty_neighbor){
////						//cout << "Has empty neighbor!" << endl;
////						// surely there's a better approach to this ....
////						if (2*x_index < kid_distance){
////							startx0 = 0;
////						}	else {
////							startx0 = 2*x_index - kid_distance;
////						}
////
////						if (2*y_index < kid_distance){
////							starty0 = 0;
////						}	else {
////							starty0 = 2*y_index - kid_distance;
////						}
////
////						if (2*z_index < kid_distance){
////							startz0 = 0;
////						}	else {
////							startz0 = 2*z_index - kid_distance;
////						}
////						for (x0 = startx0 ; x0 <  min(2*x_index + kid_distance + 1, RS_child.number_voxels_per_dim[0]); x0++){
////							for (y0 = starty0 ; y0 < min(2*y_index + kid_distance + 1, RS_child.number_voxels_per_dim[1]); y0++ ){
////								for (z0 = startz0 ; z0 < min(2*z_index + kid_distance + 1, RS_child.number_voxels_per_dim[2]); z0++ ){
////									small_index =  x0*RS_child.number_voxels_per_dim[1]*RS_child.number_voxels_per_dim[2] + y0*RS_child.number_voxels_per_dim[2] + z0;
////									RS_child.set_by_first_d[small_index] = false;
////									//RS_child.configuration_grid[small_index] = true;
////								}
////							}
////						}
////					}
//				}
//			}
//		}
//	}
//
//	// also make the "shell" part of the ones to process
//
//
//	for (int_type_t i = 0; i < RS_child.number_voxels_grid; i++){
//		RS_child.set_by_first_d[i] = false;
//		if (RS_child.set_by_first_d[i] == false){
//			voxels_to_test++;
//			RS_child.multi_dim_voxels_to_process[i % RS_child.number_divisions][i / RS_child.number_divisions] = true;
//		}
//	}
//
//	// walk trhough and find border ones ....
//
//	cout << "Voxels to test ... " << voxels_to_test << endl;
//	//char ch; cin >> ch;
//
//
//	return voxels_to_test;
//}

int_type_t ReconstructionStructure::NumberZeros(){

	int_type_t n = 0;

	for (int_type_t i = 0; i < number_voxels_grid; i++){
		if (configuration_grid[i] == false){
			n++;
		}

	}

	return n;
}


void ReconstructionStructure::UpdateNumberZeros(){

	char ch;
	int_type_t voxel_index = 0;

	int_type_t last_term_position = 0;
	int_type_t interact_index  = 0;
	int_type_t arr_index_0, arr_index_1;
	int_type_t term_index, starting_term_position;
	//int_type_t* term_index_ptr = 0;
	int_type_t thread_ID;
	int_type_t counter;

	cout << "Number of voxels " << number_voxels_grid << endl;

	// TODO set the start indices per term per div manually in the b-p part ...
	bool found = false;
	//	for (int_type_t p = 0; p < number_pixels; p++){
	//		if (number_zeros_grid[p] > 0){
	//			cout << "Found one .... this shouldn't happen " << p << endl;
	//			visual_hull_pixels[p] = true;
	//
	//			current_pixel_labeling[p] = false; // technically, can do this later ...
	//		}
	//	}

	int_type_t number_tiles = start_indices_terms_per_div_tiles.size();

	if (number_tiles == 0){
		cout << "No tiles! " << endl;
		//#pragma omp parallel private(counter, voxel_index, thread_ID, last_term_position,  interact_index, arr_index_0, arr_index_1, starting_term_position, term_index)
		{

			//#pragma omp critical
			for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
			{
				//thread_ID = omp_get_thread_num();
				voxel_index = 0;
				counter = 0;

				//#pragma omp critical
				{
					cout << "thread id " << thread_ID << endl;
					// now, update the terms .... so that the images are right
					for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

						if (x_index % 100 == 0){
							cout << "Incrementing .... " << x_index << " out of " << number_voxels_per_dim[0] << " thread ID " << thread_ID << endl;
						}
						for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
							for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){

								if (voxel_index % number_divisions == thread_ID){

									//cout << "voxel index " << endl;
									//starting_term_position = RS_child.start_indices_terms_per_div.at(thread_ID)[counter];
									//last_term_position = RS_child.start_indices_terms_per_div.at(thread_ID)[counter + 1];
									starting_term_position = start_indices_terms_per_div[thread_ID][counter];
									last_term_position = start_indices_terms_per_div[thread_ID][counter + 1];
									counter++;

									//								if (set_by_first_d[voxel_index] == false){
									//									configuration_grid[voxel_index] = false;
									//								}


									if (configuration_grid[voxel_index] == false){
										found = false;
										//cout << "Some false ... " << endl;
										arr_index_0 = starting_term_position/term_increment;
										arr_index_1 = starting_term_position % term_increment;
										//term_index_ptr = &terms_arr[arr_index_0][arr_index_1];

										//cout << "Number of terms " << last_term_position - starting_term_position << endl;
										for (interact_index = starting_term_position;
												interact_index < last_term_position; interact_index++){

											//term_index = RS_child.terms_arr_per_div.at(thread_ID).at(arr_index_0)[arr_index_1];
											term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];


											//#pragma omp critical
											{
												number_zeros_grid[term_index]++;
												//RS_child.number_zeros_grid[term_index]++;
											}
											//										if (term_index == 1074644){
											//											cout << "One occupied voxel projects to this term ... " << voxel_index << " " << number_zeros_grid[term_index] << endl;
											//											int_type_t x, y, z;
											//											ReturnXYZIndicesFromIndex(voxel_index, x, y, z);
											//											cout << "xyz " << x << ", " << y << ", " << z << endl;
											//											cout << "interact index , starting term " << interact_index << ", " << starting_term_position << endl;
											//											cout << "How many terms? " << last_term_position - starting_term_position << endl;
											//											cout << "arr0, arr1 " << arr_index_0 << ", " << arr_index_1 << endl;
											//											cout<< "term increment " << term_increment << endl;
											//											found = true;
											//										}
											//
											//										if (found){
											//											cout << "previous term ... " << terms_arr_per_div[thread_ID][arr_index_0][arr_index_1 - 1] << endl;
											//											cout << "Previous term image " << terms_arr_per_div[thread_ID][arr_index_0][arr_index_1 - 1]/(rows*cols)
											//													<< ", " << terms_arr_per_div[thread_ID][arr_index_0][arr_index_1 - 1] % (rows*cols) << endl;
											//											cout << term_index << ", " << term_index/(rows*cols) << ", " << term_index % (rows*cols) << endl;
											//											cout << "thread id " << thread_ID << endl;
											//										}

											arr_index_1++;
											//++term_index_ptr;
											if (arr_index_1 == term_increment){
												arr_index_0++;
												arr_index_1 = 0;
												//term_index_ptr = &terms_arr[arr_index_0][arr_index_1];
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
	}	else {
		cout << "Number of tiles for update number of zeros .... " << number_tiles << endl;

		if (number_tiles != 1){
			cout << "Number tiles is greater than 1 in update, not allowed " << endl;
			exit(1);
		}

		//#pragma omp parallel private(counter, voxel_index, thread_ID, last_term_position,  interact_index, arr_index_0, arr_index_1, starting_term_position, term_index)
		{

			//#pragma omp critical
			for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
			{
				//thread_ID = omp_get_thread_num();
				voxel_index = 0;
				counter = 0;

				//#pragma omp critical
				{
					cout << "thread id " << thread_ID << endl;
					// now, update the terms .... so that the images are right
					for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

						if (x_index % 100 == 0){
							cout << "Incrementing .... " << x_index << " out of " << number_voxels_per_dim[0] << " thread ID " << thread_ID << endl;
						}
						for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
							for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){

								if (voxel_index % number_divisions == thread_ID){

									//cout << "voxel index " << endl;
									//starting_term_position = RS_child.start_indices_terms_per_div.at(thread_ID)[counter];
									//last_term_position = RS_child.start_indices_terms_per_div.at(thread_ID)[counter + 1];
									starting_term_position = start_indices_terms_per_div_tiles[0][thread_ID][counter];
									last_term_position = start_indices_terms_per_div_tiles[0][thread_ID][counter + 1];
									counter++;

									//								if (set_by_first_d[voxel_index] == false){
									//									configuration_grid[voxel_index] = false;
									//								}


									if (configuration_grid[voxel_index] == false){
										found = false;
										//cout << "Some false ... " << endl;
										arr_index_0 = starting_term_position/term_increment;
										arr_index_1 = starting_term_position % term_increment;
										//term_index_ptr = &terms_arr[arr_index_0][arr_index_1];

										//cout << "Number of terms " << last_term_position - starting_term_position << endl;
										for (interact_index = starting_term_position;
												interact_index < last_term_position; interact_index++){

											//term_index = RS_child.terms_arr_per_div.at(thread_ID).at(arr_index_0)[arr_index_1];
											term_index = terms_arr_per_div_tiles[0][thread_ID][arr_index_0][arr_index_1];


											//#pragma omp critical
											{
												number_zeros_grid[term_index]++;
												//RS_child.number_zeros_grid[term_index]++;
											}
											//										if (term_index == 1074644){
											//											cout << "One occupied voxel projects to this term ... " << voxel_index << " " << number_zeros_grid[term_index] << endl;
											//											int_type_t x, y, z;
											//											ReturnXYZIndicesFromIndex(voxel_index, x, y, z);
											//											cout << "xyz " << x << ", " << y << ", " << z << endl;
											//											cout << "interact index , starting term " << interact_index << ", " << starting_term_position << endl;
											//											cout << "How many terms? " << last_term_position - starting_term_position << endl;
											//											cout << "arr0, arr1 " << arr_index_0 << ", " << arr_index_1 << endl;
											//											cout<< "term increment " << term_increment << endl;
											//											found = true;
											//										}
											//
											//										if (found){
											//											cout << "previous term ... " << terms_arr_per_div[thread_ID][arr_index_0][arr_index_1 - 1] << endl;
											//											cout << "Previous term image " << terms_arr_per_div[thread_ID][arr_index_0][arr_index_1 - 1]/(rows*cols)
											//													<< ", " << terms_arr_per_div[thread_ID][arr_index_0][arr_index_1 - 1] % (rows*cols) << endl;
											//											cout << term_index << ", " << term_index/(rows*cols) << ", " << term_index % (rows*cols) << endl;
											//											cout << "thread id " << thread_ID << endl;
											//										}

											arr_index_1++;
											//++term_index_ptr;
											if (arr_index_1 == term_increment){
												arr_index_0++;
												arr_index_1 = 0;
												//term_index_ptr = &terms_arr[arr_index_0][arr_index_1];
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

	sfis_error = function_constant;
	cout << "Completed the increment ... " << endl;
	//char fg; cin >> fg;
	for (int_type_t p = 0; p < number_pixels; p++){
		visual_hull_pixels[p] = false;
		if (number_zeros_grid[p] > 0){
			//cout << "Found one .... " << p << endl;

			current_pixel_labeling[p] = false; // technically, can do this later ...

			//			if (p == 839*cols + 724){
			//				cout << "yes, term " << 839*cols + 724 << " is occupied " << number_zeros_grid[p] << endl;
			//			}
			//			if (number_zeros_grid[p] == 1){
			//				cout << " a zeros count that is one ... " << p << endl;
			//			}
		}	else {
			current_pixel_labeling[p] = true;
			sfis_error += coefficients[p];
		}
	}

	//cout << "Done with update" << endl;


	//	for (int_type_t i = 0; i < number_cameras; i++){
	//
	//
	//
	//
	//
	//	}

	//cin >> ch;

}


void ReconstructionStructure::CheckProjectionForErrors(){

	char ch;
	int_type_t voxel_index = 0;

	int_type_t last_term_position = 0;
	int_type_t previous_term_position = 0;
	int_type_t interact_index  = 0;
	int_type_t arr_index_0, arr_index_1;
	int_type_t term_index, starting_term_position;
	//int_type_t* term_index_ptr = 0;
	int_type_t thread_ID;
	int_type_t counter;

	cout << "Number of voxels " << number_voxels_grid << endl;

	// TODO set the start indices per term per div manually in the b-p part ...
	bool found = false;
	//	for (int_type_t p = 0; p < number_pixels; p++){
	//		if (number_zeros_grid[p] > 0){
	//			cout << "Found one .... this shouldn't happen " << p << endl;
	//			visual_hull_pixels[p] = true;
	//
	//			current_pixel_labeling[p] = false; // technically, can do this later ...
	//		}
	//	}

	//#pragma omp parallel private(counter, voxel_index, thread_ID, last_term_position,  interact_index, arr_index_0, arr_index_1, starting_term_position, term_index)
	{

		//#pragma omp critical
		//for (thread_ID = 0; thread_ID < number_divisions; thread_ID++)
		for (thread_ID = 0; thread_ID < 1; thread_ID++)
		{
			//thread_ID = omp_get_thread_num();
			voxel_index = 0;
			counter = 0;

			//#pragma omp critical
			{
				cout << "thread id " << thread_ID << endl;
				// now, update the terms .... so that the images are right
				for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

					if (x_index % 100 == 0){
						cout << "Checking for errors  .... " << x_index << " out of " << number_voxels_per_dim[0] << " thread ID " << thread_ID << endl;
					}
					for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
						for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){

							if (voxel_index % number_divisions == thread_ID){

								//cout << "voxel index " << endl;

								starting_term_position = start_indices_terms_per_div[thread_ID][counter];
								last_term_position = start_indices_terms_per_div[thread_ID][counter + 1];
								counter++;

								//								if (set_by_first_d[voxel_index] == false){
								//									configuration_grid[voxel_index] = false;
								//								}


								//if (configuration_grid[voxel_index] == false)
								{
									found = false;
									//cout << "Some false ... " << endl;
									arr_index_0 = starting_term_position/term_increment;
									arr_index_1 = starting_term_position % term_increment;
									//term_index_ptr = &terms_arr[arr_index_0][arr_index_1];

									//cout << "Number of terms " << last_term_position - starting_term_position << endl;
									for (interact_index = starting_term_position;
											interact_index < last_term_position; interact_index++){

										//term_index = RS_child.terms_arr_per_div.at(thread_ID).at(arr_index_0)[arr_index_1];
										term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];

										if (interact_index == starting_term_position){
											previous_term_position = term_index;
										}	else {

											// compare term and previous term ...
											if (term_index < previous_term_position){
												cout << "ERROR thread " << thread_ID << "  " << previous_term_position << " " << term_index << endl;
											}
											previous_term_position = term_index;
										}

										if (voxel_index == 16845 && (interact_index - starting_term_position) > 4000){
											cout << "Term " << term_index << endl;
										}


										arr_index_1++;
										//++term_index_ptr;
										if (arr_index_1 == term_increment){
											arr_index_0++;
											arr_index_1 = 0;
											//term_index_ptr = &terms_arr[arr_index_0][arr_index_1];
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


	cout << "Completed the increment ... " << endl;
	//char fg; cin >> fg;
	for (int_type_t p = 0; p < number_pixels; p++){
		if (number_zeros_grid[p] > 0){
			//cout << "Found one .... " << p << endl;
			visual_hull_pixels[p] = true;

			current_pixel_labeling[p] = false; // technically, can do this later ...

			if (p == 839*cols + 724){
				cout << "yes, term " << 839*cols + 724 << " is occupied " << number_zeros_grid[p] << endl;
			}
			//			if (number_zeros_grid[p] == 1){
			//				cout << " a zeros count that is one ... " << p << endl;
			//			}
		}
	}

	//cout << "Done with update" << endl;


	//	for (int_type_t i = 0; i < number_cameras; i++){
	//
	//
	//
	//
	//
	//	}

	cout << "After checking for errors ... " << endl;
	cin >> ch;

}



void ReconstructionStructure::PrintFeaturesToFile(std::ofstream& out){

	out << "inter voxel distance " << inter_voxel_distance << endl;
	out << endl << "Number of voxels per dim " << endl;
	for (int i = 0; i < 3; i++){
		out << number_voxels_per_dim[i] << endl;
	}

	out << endl << "Bounding box " << endl;

	for (int i = 0; i < 3; i++){
		out << BB[0][i] << "  " << BB[1][i] << endl;
	}

	out << endl << "Number of voxels in the grid: " << FormatWithCommas<int_type_t>(number_voxels_grid) << endl;

	//out << FormatWithCommas<int_type_t>(number_voxels_grid) << endl;

	out << "Number of pixels " << FormatWithCommas<int>(number_pixels) << endl;
}


void ReconstructionStructure::ConstructSfISCostFunction(vector<ImageData*>& cameras, int downsample){

	// resize the coefficient vector....
	sfis_error = 0;
	coefficients.resize(number_pixels, 0);
	visual_hull_pixels.resize(number_pixels, false);
	function_constant = 0;

	int im_value;
	int prob_occ;
	int prob_empty;

	bool color;

	//sleep(1);
	// assume gray ....
	cv::Mat im(rows, cols, CV_8UC1, cv::Scalar(0,0,0));

	cout << "Number of cameras in construct sfis cost function " <<  number_cameras << endl;
	for (int_type_t cam = 0, pixel_index = 0; cam < number_cameras; cam++){

		uchar* p = 0;
		// resize for each time ...  get rid fo resize -- resize in the
		if (downsample > 1){
			cv::resize(cameras[cam]->imcv, im, im.size());
			//cameras[cam]->imcv = im.clone;
			//cv::imwrite("temp1.png", im);
			p = &im.data[0];
		}	else {
			p = &cameras[cam]->imcv.data[0];
			//cv::imwrite("temp2.png", cameras[cam]->imcv);
		}


		color = cameras[cam]->imcv.channels() == 3;


		//uchar* p = &cameras[cam]->imcv.data[0];
		cout << "camera " << cam << " color " << color << endl;

		if (color){
			for (int r = 0; r < rows; r++){
				for (int c = 0; c < cols; c++, pixel_index++){

					if (is_in_undistorted_image[cam*(rows*cols) + r*cols + c]){
						im_value = int(*p);

						prob_occ = im_value;
						prob_empty = 255 - prob_occ;

						prob_occ /= 2;
						prob_empty /=2;

						// we divide by two so that we can get away with using int8_t integers.
						function_constant +=  prob_empty;

						// ditto on dividing by two -- faster with an iterator?
						coefficients[cam*(rows*cols) + r*cols + c] = (prob_occ - prob_empty);
						//coefficients[pixel_index] = (prob_occ - prob_empty);

						sfis_error += prob_occ;
						// increment p by 3 ....

					}
					p++;
					p++;
					p++;
				}
			}
		}	else {

			for (int r = 0; r < rows; r++){
				//cout << "Rows " << r << endl;
				for (int c = 0; c < cols; c++, pixel_index++){

					if (is_in_undistorted_image[cam*(rows*cols) + r*cols + c]){
						im_value = int(*p);

						prob_occ = im_value;
						prob_empty = 255 - prob_occ;

						prob_occ /= 2;
						prob_empty /=2;

						// we divide by two so that we can get away with using int8_t integers.
						function_constant +=  prob_empty;//ʲô  ˼      ʼ  Ϊ   б      ص      ֵ

						// ditto on dividing by two -- faster with an iterator?
						coefficients[cam*(rows*cols) + r*cols + c] = (prob_occ - prob_empty);
						//coefficients[cam*(rows*cols) + r*cols + c] = min(127, (prob_occ - prob_empty) + 127/4);
						//coefficients[pixel_index] = (prob_occ - prob_empty);

						sfis_error += prob_occ;//  ʼ  Ϊ    Ŀ     ص      ֵ  һ  
						// increment p by 3 ....

					}
					
					p++;

				}
			}
		}
	}


}

void ReconstructionStructure::ConstructSfISCostFunctionUsingAncestor(ReconstructionStructure* ancestor){

	// resize the coefficient vector....
	sfis_error = 0;
	//coefficients.resize(number_pixels, 0);
	visual_hull_pixels.resize(number_pixels, false);
	function_constant = ancestor->function_constant;

	coefficients = ancestor->coefficients;

	sfis_error = function_constant;
}


void ReconstructionStructure::ConstructSfISCostFunction(vector<ImageData*>& cameras){

	// resize the coefficient vector....
	sfis_error = 0;
	coefficients.resize(number_pixels, 0);
	visual_hull_pixels.resize(number_pixels, false);
	function_constant = 0;

	int im_value;
	int prob_occ;
	int prob_empty;

	//sleep(1);

	cout << "Number of cameras in construct sfis cost function " <<  number_cameras << endl;
	for (int_type_t cam = 0, pixel_index = 0; cam < number_cameras; cam++){

		uchar* p = &cameras[cam]->imcv.data[0];
		cout << "camera " << cam << endl;

		for (int r = 0; r < rows; r++){
			for (int c = 0; c < cols; c++, pixel_index++){

				if (is_in_undistorted_image[cam*(rows*cols) + r*cols + c]){
					im_value = int(*p);

					prob_occ = im_value;
					prob_empty = 255 - prob_occ;

					prob_occ /= 2;
					prob_empty /=2;

					// we divide by two so that we can get away with using int8_t integers.
					function_constant +=  prob_empty;

					// ditto on dividing by two -- faster with an iterator?
					coefficients[cam*(rows*cols) + r*cols + c] = (prob_occ - prob_empty);
					//coefficients[pixel_index] = (prob_occ - prob_empty);

					sfis_error += prob_occ;
					// increment p by 3 ....

				}
				p++;
				p++;
				p++;
			}
		}
	}


}

//inline bool ComputeIntersectionOfRayWithPlane(const vector<double>& C, const vector<double>& V, int plane_axis, double plane_value, const Bbox_3& BB,
//		vector<double>& lambdaV, vector<double>& X){



//void ReconstructionStructure::VisualHull(bool verbose, ofstream& out){
//
//
//	// at this point, this is the visual hull
//	int number_of_flips = 0;
//	int_type_t p = 0;
//	int_type_t interacting_p_within_image;
//	int_type_t interacting_p;
//	double plane_value;
//	int64_t first_derivative = 0;
//	bool current_value;
//	bool current_value_interacting;
//	int camera_number;
//	int within_image_index;
//	vector<int_type_t> list_of_interacting_pixels;
//	vector<double> lambdaV(3, 0);
//	vector<double> V(3, 0);
//	vector<double> Xhm(4, 1);
//	vector<double> x(3, 1);
//	bool valid_point;
//	bool in_image;
//	bool flip;
//
//	vector<vector<double> > Pvv(3, vector<double>(4, 0));
//	vector<vector<double> > Pivv(3, vector<double>(3, 0));
//	vector<double> Cvv(3);
//	vector<double> Dvv(3);
//	vector<double> lambDvv(3);
//	vector<double> Xvv(3);
//	vector<double> xvv(3);
//	vector<double> Rvv(3);
//
//	int in_count = 0;
//
//
//	vector<bool> tested_this_round(number_pixels, false);
//
//	//__gnu_parallel::random_shuffle(random_pixels.begin(), random_pixels.end());
//
//
//	// randomly pick a plane ... scramble first ....
//	for (int i = 0; i < 3; i++){
//		//	__gnu_parallel::random_shuffle(all_planes[i].begin(), all_planes[i].end());
//	}
//
//	// so we're going to choose according to plane, and then walk through all of the pixels for that plane.
//
//
//	// TODo there's a erroneous access somewhere in here ....
//	// AND we need to set the min number of cameras parameter
//	// AND check on whether the 1st Ds are being computed correctly ....
//	// PLSU the function value is going wonky ....
//
//
//
//	// To try -- compare the occupied versus the number of zeros business to make sure that is correct.
//	vector<bool> list_of_current_values;
//
//
//	// a test ... project 0, 0, 0, 1 to each image ....
//	Xhm[0] =0;
//	Xhm[1] = 0;
//	Xhm[2] = 0;
//	Xhm[3] = 1;
//
//	vector< vector<double> > image_xs;
//
//	//	for (int cam_index = 0; cam_index < number_cameras; cam_index++){
//	//
//	//
//	//		in_image = ProjectPointAndReturnIndex(projection_matrices[cam_index],
//	//				Xhm, x, rows, cols, interacting_p_within_image);
//	//
//	//		cout << "Above was in the image? " << endl;
//	//		cout << "p within image " << interacting_p_within_image << endl;
//	//		cout << "coefficient value " << int(coefficients[cam_index*rows*cols + interacting_p_within_image]);
//	//		cout << "---------------------------" << endl;
//	//
//	//		image_xs.push_back(x);
//	//		list_of_interacting_pixels.push_back(interacting_p_within_image);
//	//
//	//		//			if (in_image){
//	//		//
//	//		//			}
//	//	}
//	//
//	//	for (int cam_index = 0; cam_index < number_cameras; cam_index++){
//	//
//	//		within_image_index =list_of_interacting_pixels[cam_index];
//	//
//	//		V[0] = pixel_vectors[cam_index][within_image_index*3];
//	//		V[1] = pixel_vectors[cam_index][within_image_index*3 + 1];
//	//		V[2] = pixel_vectors[cam_index][within_image_index*3 + 2];
//	//
//	//		valid_point = ComputeIntersectionOfRayWithPlane(camera_centers[cam_index],
//	//				V, 0, 0, BB, lambdaV, Xhm);
//	//		Xhm[3] = 1;
//	//
//	//		cout << "Reprojected for " << cam_index << endl;
//	//		cout << Xhm[0] << " " << Xhm[1] << " " << Xhm[2] << " " << Xhm[3] << endl;
//	//		cout << "V  " << V[0] << " " << V[1] << " " << V[2] << endl;
//	//		cout << "-------------------------------------" << endl;
//	//	}
//
//	// do the opposite ....
//
//	//exit(1);
//
//	bool some_negative;
//	// an east test -- flip all to occupied, the write.  See what happens.
//	for (int axis = 0; axis < 3; axis++){
//		// for beginning tests, choose the center plane...
//		//for (int plane_index = all_planes[axis].size()/2- 2 , np = all_planes[axis].size()/2 + 3; plane_index < np; plane_index++){
//		for (int plane_index = 0, np = all_planes[axis].size(); plane_index < np; plane_index++){
//			plane_value = all_planes[axis][plane_index];
//
//
//			cout << "axis is " << axis << endl;
//			cout << "value is " << plane_value << endl;
//
//			//char lk; cin >> lk;
//
//
//
//			// pick the pixel in the middle of the first image .... (r/2)*cols + cols/2
//			//int start_index = (rows/2)*cols + cols/2;
//			// there's still an indexing error somewhere ....
//			for (int pixel_index = 0; pixel_index < number_pixels; pixel_index++){
//				p = random_pixels[pixel_index];
//
//				//				if (coefficients[p] < 0){
//				//					cout << "LEss zero coefficient! " << int(coefficients[p]) <<endl; //char ch; cin >> ch;
//				//				}
//
//				if (pixel_index % 1000000 == 0 && verbose){
//					cout << "p: " << p << "  pixel index " << pixel_index << " out of " << number_pixels << " with number in " << in_count << endl;
//					cout << "Number of flips " << number_of_flips << endl;
//					//cout << "Number of zeros already " << number_of_zeros[p] << endl;
//					//cout << "V " << V[0] << " " << V[1] << " " << V[2] << endl;
//				}
//
//
//				//if (!tested_this_round[p])
//				{
//
//					// what's the originating camera?
//					camera_number = p / (rows*cols);
//
//					// within image index ...
//					within_image_index = p % (rows*cols);
//
//					current_value =
//							occupied_point_representation.at(camera_number*3 + axis).at(rows*cols*plane_index
//									+ within_image_index);
//
//					V[0] = pixel_vectors[camera_number][within_image_index*3];
//					V[1] = pixel_vectors[camera_number][within_image_index*3 + 1];
//					V[2] = pixel_vectors[camera_number][within_image_index*3 + 2];
//
//					valid_point = ComputeIntersectionOfRayWithPlane(camera_centers[camera_number],
//							V, axis, plane_value, BB, lambdaV, Xhm);
//					Xhm[3] = 1;
//
//					//					cout << "Valid point " << valid_point << endl;
//					//					cout << "seed " << endl;
//					//					cout << "pixel p " << p << "   projects to " << Xhm[0] << " " << Xhm[1]
//					//					                                                                  << " " << Xhm[2] << " " << Xhm[3] << endl;
//					//					cout << "V  " << V[0] << " " << V[1] << " " << V[2] << endl;
//					//					cout << "Camera number " << camera_number << endl;
//					//					cout << "winitn image index " << within_image_index << endl;
//					//					cout << "current value " << current_value << endl;
//					//					cout << "Number of zeros " << number_of_zeros[p] << endl;
//					//					cout << "-----------------------------------------------" << endl;
//					list_of_current_values.clear();
//					if (valid_point){
//						in_count++;
//
//						// make a list of all interacting terms with this one
//						list_of_interacting_pixels.push_back(p);
//
//						//						cout << "seed " << endl;
//						//						cout << "pixel p " << p << "projects to " << Xhm[0] << " " << Xhm[1]
//						//						                                                                  << " " << Xhm[2] << " " << Xhm[3] << endl;
//						//
//						//						cout << "Camera number " << camera_number << endl;
//						//						cout << "winitn image index " << within_image_index << endl;
//						//						cout << "current value " << current_value << endl;
//						//						cout << "Number of zeros " << number_of_zeros[p] << endl;
//						//						cout << "-----------------------------------------------" << endl;
//						// project to all other cameras .... record pixel.
//						for (int cam_index = 0; cam_index < number_cameras; cam_index++){
//
//							if (cam_index != camera_number){
//
//								in_image = ProjectPointAndReturnIndex(projection_matrices[cam_index],
//										Xhm, x, rows, cols, interacting_p_within_image);
//
//								if (in_image){
//
//
//									//char fg; cin >> fg;
//
//									interacting_p = rows*cols*cam_index + interacting_p_within_image;
//
//									list_of_interacting_pixels.push_back(interacting_p);
//									//									cout << "In image! for camera " << cam_index << endl;
//									//									cout << "within image index " << interacting_p_within_image << endl;
//									//									cout << "Global index " << interacting_p << endl;
//									//									cout << "---------------------------------------------" << endl;
//								}
//
//							}
//						}
//
//						//exit(1);
//
//
//						//.... sum up the first d values according to the number zeros vector
//						first_derivative = 0;
//
//						flip = false;
//
//						some_negative = false;
//
//						if (list_of_interacting_pixels.size() > 1){
//
//							list_of_current_values.push_back(current_value);
//
//							for (int pv = 0, pvn = list_of_interacting_pixels.size(); pv < pvn; pv++){
//
//								interacting_p = list_of_interacting_pixels[pv];
//
//								// we need to get the current value of this pixel
//								current_value_interacting = occupied_point_representation.at((interacting_p / (rows*cols))*3 +
//										axis).at(rows*cols*plane_index
//												+ interacting_p % (rows*cols));
//								list_of_current_values.push_back(current_value_interacting);
//
//								//								cout << "list of interacting .... " << endl;
//								//								cout << "index " << interacting_p << endl;
//								//								cout << "current value " << current_value_interacting << endl;
//								//								cout << "number of zeros " << number_of_zeros[interacting_p] << endl;
//								//								cout << "coefficients " << int(coefficients[interacting_p]) << endl;
//								if (coefficients.at(interacting_p) < 0){
//									some_negative = true;
//								}
//
//								// this was it for the visual hulll -- there's something wrong with the
//								// other formulation.
//								//if ( (number_of_zeros.at(interacting_p) == 0) ){
//								// okay, got it.
//								if ( (number_of_zeros.at(interacting_p) == 0) ||
//										(!current_value && (current_value_interacting == false)
//												&& (number_of_zeros.at(interacting_p) == 1))){
//
//									// this sometimes happens that the interacting p is actually (occupied), but that
//									// the seed value is empty, b/c of the projections to image planes.
//
//									// this is not quite right -- the false pixel has to be the CURRENT one
//
//
//									first_derivative += int(coefficients.at(interacting_p));
//									//cout << "updated first d! " << first_derivative << endl;
//								}
//
//								//								if ((current_value_interacting == true)	&& (number_of_zeros.at(interacting_p) > 0)){
//								//									cout << "Error!  Mismatch in current value and number zeros A" << endl;
//								//									exit(1);
//								//								}
//
//								if ((current_value_interacting == false) && (number_of_zeros.at(interacting_p) == 0)){
//									cout << "Error!  Mismatch in current value and number zeros B" << endl;
//									exit(1);
//								}
//							}
//
//
//							// are we going to flip?  Change the value for self and all others that interact with self, update flipped for this pixel vector
//							// update the error estimate.
//							flip = !some_negative && current_value;
//							//flip = (current_value && first_derivative >= 0); // || (!current_value && first_derivative < 0);
//						}
//
//						if (flip){
//
//							//if (sfis_error/127 < 10000){
//							//cout << "sfis error " << sfis_error/127 << endl;
//							//}
//							sfis_error -= fabs(first_derivative);
//							//							if (first_derivative > 0){
//							//								sfis_error -= first_derivative;
//							//							}	else {
//							//								sfis_error += first_derivative;
//							//							}
//
//							//							if (first_derivative != 0){
//							//								cout << "first d " << first_derivative << endl;
//							//							}
//
//							// this covers self as well.
//							for (int pv = 0, pvn = list_of_interacting_pixels.size(); pv < pvn; pv++){
//
//								interacting_p = list_of_interacting_pixels[pv];
//
//								// update the occupied points
//								occupied_point_representation.at((interacting_p / (rows*cols))*3 +
//										axis).at(rows*cols*plane_index
//												+ interacting_p % (rows*cols)) = !current_value;
//
//
//								// if the current value of the interacting pixel is already 0 or 1, we don't update.
//
//								// update the number of zeros vector.?
//								// if first d > 0, number fo zeros goes up
//								//number_of_zeros.at(interacting_p)++;
//
//
//								/// TODO -- before monkeying with this, check whether the point was already true or false,
//								// b/c of the possible mismatch via projection matrices
//								if (first_derivative >= 0){
//									number_of_zeros.at(interacting_p)++;
//								}	else {
//									number_of_zeros.at(interacting_p)--;
//								}
//
//								//cout << "updating number of zeros for " << interacting_p << " to " << number_of_zeros[interacting_p] << endl;
//								//cout << "new value is " << occupied_point_representation[(interacting_p / (rows*cols))*3 +
//								//                                                        axis][rows*cols*plane_index
//								//                                                             + interacting_p % (rows*cols)] << endl;
//
//								if (number_of_zeros.at(interacting_p) < 0){
//									cout << "Error!  N of zeros going less than 0 "<< endl;
//									exit(1);
//								}
//
//								tested_this_round.at(interacting_p) = true;
//							}
//
//							//pixel_index = 1000000000;
//
//							number_of_flips++;
//							//cout << "Flipped for " << p << endl;
//							//cout << "number of flips " << number_of_flips << endl;
//							//exit(1);
//						}
//
//
//						list_of_interacting_pixels.clear();
//
//					}
//
//				}
//
//			}
//		}
//	}
//
//}
//
//void ReconstructionStructure::VisualHullSetEmpties(bool verbose, ofstream& out){
//
//
//	// at this point, this is the visual hull
//	int number_of_flips = 0;
//	int_type_t p = 0;
//	int_type_t interacting_p_within_image;
//	int_type_t interacting_p;
//	double plane_value;
//	int64_t first_derivative = 0;
//	bool current_value;
//	bool current_value_interacting;
//	int camera_number;
//	int within_image_index;
//	vector<int_type_t> list_of_interacting_pixels;
//	vector<double> lambdaV(3, 0);
//	vector<double> V(3, 0);
//	vector<double> Xhm(4, 1);
//	vector<double> x(3, 1);
//	bool valid_point;
//	bool in_image;
//	bool flip;
//
//	vector<vector<double> > Pvv(3, vector<double>(4, 0));
//	vector<vector<double> > Pivv(3, vector<double>(3, 0));
//	vector<double> Cvv(3);
//	vector<double> Dvv(3);
//	vector<double> lambDvv(3);
//	vector<double> Xvv(3);
//	vector<double> xvv(3);
//	vector<double> Rvv(3);
//
//	int in_count = 0;
//
//
//	vector<bool> tested_this_round(number_pixels, false);
//
//	//__gnu_parallel::random_shuffle(random_pixels.begin(), random_pixels.end());
//
//
//	// randomly pick a plane ... scramble first ....
//	for (int i = 0; i < 3; i++){
//		//	__gnu_parallel::random_shuffle(all_planes[i].begin(), all_planes[i].end());
//	}
//
//
//	int number_rigid = 0;
//	bool some_positive;
//	// an east test -- flip all to occupied, the write.  See what happens.
//	for (int axis = 0; axis < 3; axis++){
//		// for beginning tests, choose the center plane...
//		//for (int plane_index = all_planes[axis].size()/2- 2 , np = all_planes[axis].size()/2 + 3; plane_index < np; plane_index++){
//		for (int plane_index = 0, np = all_planes[axis].size(); plane_index < np; plane_index++){
//			plane_value = all_planes[axis][plane_index];
//
//
//			// this can be parallel b/c not dependent ....
//			cout << "axis is " << axis << endl;
//			cout << "value is " << plane_value << endl;
//
//
//			// pick the pixel in the middle of the first image .... (r/2)*cols + cols/2
//			//int start_index = (rows/2)*cols + cols/2;
//			// there's still an indexing error somewhere ....
//			for (int pixel_index = 0; pixel_index < number_pixels; pixel_index++){
//				p = random_pixels[pixel_index];
//
//				//				if (coefficients[p] < 0){
//				//					cout << "LEss zero coefficient! " << int(coefficients[p]) <<endl; //char ch; cin >> ch;
//				//				}
//
//				if (pixel_index % 1000000 == 0 && verbose){
//					cout << "EMPTIES p: " << p << "  pixel index " << pixel_index << " out of " << number_pixels << " with number in " << in_count << endl;
//					cout << "Number rigid " << number_rigid << endl;
//					//cout << "Number of flips " << number_of_flips << endl;
//					//cout << "Number of zeros already " << number_of_zeros[p] << endl;
//					//cout << "V " << V[0] << " " << V[1] << " " << V[2] << endl;
//				}
//
//
//				//if (!tested_this_round[p])
//				{
//
//					// what's the originating camera?
//					camera_number = p / (rows*cols);
//
//					// within image index ...
//					within_image_index = p % (rows*cols);
//
//					current_value =
//							occupied_point_representation.at(camera_number*3 + axis).at(rows*cols*plane_index
//									+ within_image_index);
//
//					V[0] = pixel_vectors[camera_number][within_image_index*3];
//					V[1] = pixel_vectors[camera_number][within_image_index*3 + 1];
//					V[2] = pixel_vectors[camera_number][within_image_index*3 + 2];
//
//					valid_point = ComputeIntersectionOfRayWithPlane(camera_centers[camera_number],
//							V, axis, plane_value, BB, lambdaV, Xhm);
//					Xhm[3] = 1;
//
//					list_of_interacting_pixels.clear();
//					if (valid_point){
//						in_count++;
//
//						// make a list of all interacting terms with this one
//						list_of_interacting_pixels.push_back(p);
//
//						// project to all other cameras .... record pixel.
//						for (int cam_index = 0; cam_index < number_cameras; cam_index++){
//
//							if (cam_index != camera_number){
//
//								in_image = ProjectPointAndReturnIndex(projection_matrices[cam_index],
//										Xhm, x, rows, cols, interacting_p_within_image);
//
//								if (in_image){
//
//
//									//char fg; cin >> fg;
//
//									interacting_p = rows*cols*cam_index + interacting_p_within_image;
//
//									list_of_interacting_pixels.push_back(interacting_p);
//								}
//
//							}
//						}
//
//						//exit(1);
//
//
//						//.... sum up the first d values according to the number zeros vector
//						first_derivative = 0;
//
//						flip = false;
//
//						some_positive = false;
//
//						if (list_of_interacting_pixels.size() > 1){
//
//							for (int pv = 0, pvn = list_of_interacting_pixels.size(); pv < pvn; pv++){
//
//								interacting_p = list_of_interacting_pixels[pv];
//
//								// occluded by the visual hull here ....
//								if (number_of_zeros.at(interacting_p) == 0 && coefficients.at(interacting_p) > 0){
//									some_positive = true;
//								}
//							}
//
//							if (!some_positive && rigid_label_pixels[p] == false){
//
//								rigid_label_pixels[p] =  true;
//								number_rigid++;
//
//							}
//						}
//
//
//
//					}
//
//				}
//
//			}
//		}
//	}
//
//}

//void ReconstructionStructure::OneRoundOfLocalSearch(bool verbose, ofstream& out){
//
//
//	// at this point, this is the visual hull
//	int number_of_flips = 0;
//	int_type_t p = 0;
//	int_type_t interacting_p_within_image;
//	int_type_t interacting_p;
//	double plane_value;
//	int64_t first_derivative_label_true = 0;
//	int64_t first_derivative_label_false = 0;
//	bool current_value;
//	bool current_value_interacting;
//	int camera_number;
//	int within_image_index;
//	vector<int_type_t> list_of_interacting_pixels_label_true;
//	vector<int_type_t> list_of_interacting_pixels_label_false;
//	vector<double> lambdaV(3, 0);
//	vector<double> V(3, 0);
//	vector<double> Xhm(4, 1);
//	vector<double> x(3, 1);
//	bool valid_point;
//	bool in_image;
//	bool flip;
//
//	vector<vector<double> > Pvv(3, vector<double>(4, 0));
//	vector<vector<double> > Pivv(3, vector<double>(3, 0));
//	vector<double> Cvv(3);
//	vector<double> Dvv(3);
//	vector<double> lambDvv(3);
//	vector<double> Xvv(3);
//	vector<double> xvv(3);
//	vector<double> Rvv(3);
//
//	bool one_cannot_be_changed = false;
//
//	int in_count = 0;
//
//
//	vector<bool> tested_this_round(number_pixels, false);
//
//	__gnu_parallel::random_shuffle(random_pixels.begin(), random_pixels.end());
//
//
//	// randomly pick a plane ... scramble first ....
//	for (int i = 0; i < 3; i++){
//		__gnu_parallel::random_shuffle(all_planes[i].begin(), all_planes[i].end());
//	}
//
//	// so we're going to choose according to plane, and then walk through all of the pixels for that plane.
//
//
//	// TODo there's a erroneous access somewhere in here ....
//	// AND we need to set the min number of cameras parameter
//	// AND check on whether the 1st Ds are being computed correctly ....
//	// PLSU the function value is going wonky ....
//
//
//
//	// To try -- compare the occupied versus the number of zeros business to make sure that is correct.
//	//vector<bool> list_of_current_values;
//
//	// an east test -- flip all to occupied, the write.  See what happens.
//	for (int axis = 0; axis < 3; axis++){
//		// for beginning tests, choose the center plane...
//		//for (int plane_index = all_planes[axis].size()/2- 2 , np = all_planes[axis].size()/2 + 3; plane_index < np; plane_index++){
//		for (int plane_index = 0, np = all_planes[axis].size(); plane_index < np; plane_index++){
//			plane_value = all_planes[axis][plane_index];
//
//
//			cout << "axis is " << axis << endl;
//			cout << "value is " << plane_value << endl;
//
//			// pick the pixel in the middle of the first image .... (r/2)*cols + cols/2
//			//int start_index = (rows/2)*cols + cols/2;
//			// there's still an indexing error somewhere ....
//			for (int pixel_index = 0; pixel_index < number_pixels; pixel_index++){
//				p = random_pixels[pixel_index];
//
//				//				if (coefficients[p] < 0){
//				//					cout << "LEss zero coefficient! " << int(coefficients[p]) <<endl; //char ch; cin >> ch;
//				//				}
//
//				if (pixel_index % 1000000 == 0){
//					cout << "p: " << p << "  pixel index " << pixel_index << " out of " << number_pixels << " with number in " << in_count << endl;
//					cout << "Number of flips " << number_of_flips << endl;
//					cout << "error " << sfis_error/127 << endl;
//					//cout << "V " << V[0] << " " << V[1] << " " << V[2] << endl;
//				}
//
//
//				//if (!tested_this_round[p])
//				{
//
//					// what's the originating camera?
//					camera_number = p / (rows*cols);
//
//					// within image index ...
//					within_image_index = p % (rows*cols);
//
//					// TODO -- later, see if changing the access makes this faster ..
//
//
//					//current_value =
//					//		occupied_point_representation[camera_number*3 + axis][rows*cols*plane_value + within_image_index];
//
//					// Compute the point
//					//					V[0] = pixel_vectors.at(camera_number).at(within_image_index*3);
//					//										V[1] = pixel_vectors.at(camera_number).at(within_image_index*3 + 1);
//					//										V[2] = pixel_vectors.at(camera_number).at(within_image_index*3 + 2);
//					V[0] = pixel_vectors[camera_number][within_image_index*3];
//					V[1] = pixel_vectors[camera_number][within_image_index*3 + 1];
//					V[2] = pixel_vectors[camera_number][within_image_index*3 + 2];
//
//
//
//					valid_point = ComputeIntersectionOfRayWithPlane(camera_centers[camera_number],
//							V, axis, plane_value, BB, lambdaV, Xhm);
//					Xhm[3] = 1;
//
//
//
//					list_of_interacting_pixels_label_false.clear();
//					list_of_interacting_pixels_label_true.clear();
//					if (valid_point){
//						in_count++;
//
//						current_value =
//								occupied_point_representation.at(camera_number*3 + axis).at(rows*cols*plane_index
//										+ within_image_index);
//
//
//						// make a list of all interacting terms with this one
//						if (current_value){
//							list_of_interacting_pixels_label_true.push_back(p);
//						}	else {
//							list_of_interacting_pixels_label_false.push_back(p);
//						}
//
//						first_derivative_label_false = 0;
//						first_derivative_label_true = 0;
//						// project to all other cameras .... record pixel.
//						for (int cam_index = 0; cam_index < number_cameras; cam_index++){
//
//							if (cam_index != camera_number){
//
//								in_image = ProjectPointAndReturnIndex(projection_matrices[cam_index],
//										Xhm, x, rows, cols, interacting_p_within_image);
//
//								if (in_image){
//
//									interacting_p = rows*cols*cam_index + interacting_p_within_image;
//
//
//									// only push back if the value is the same as self .... b/c the first d, changing value, etc. depends on this.
//									current_value_interacting = occupied_point_representation.at((interacting_p / (rows*cols))*3 +
//											axis).at(rows*cols*plane_index
//													+ interacting_p % (rows*cols));
//
//									if (current_value_interacting){
//										list_of_interacting_pixels_label_true.push_back(interacting_p);
//									}	else {
//										list_of_interacting_pixels_label_false.push_back(interacting_p);
//									}
//								}
//
//							}
//						}
//
//						//exit(1);
//
//
//						//.... sum up the first d values according to the number zeros vector
//
//						flip = false;
//
//						bool some_negative = false;
//						bool some_positive = false;
//						one_cannot_be_changed = false;
//						int8_t current_coeff = 0;
//
//						if (list_of_interacting_pixels_label_true.size() > 0){
//
//							for (int pv = 0, pvn = list_of_interacting_pixels_label_true.size(); pv < pvn; pv++){
//
//								interacting_p = list_of_interacting_pixels_label_true[pv];
//
//								current_coeff =  int(coefficients.at(interacting_p));
//
//								if ( number_of_zeros.at(interacting_p) == 0 ){
//									first_derivative_label_true += current_coeff;
//									if (current_coeff < 0){
//										some_negative =  true;
//									}
//
//									if (current_coeff > 0){
//										some_positive = true;
//									}
//								}
//
//								if (rigid_label_pixels[interacting_p]){
//									one_cannot_be_changed = true;
//								}
//
//
//							}
//						}
//
//						if (list_of_interacting_pixels_label_false.size() > 0){
//
//							for (int pv = 0, pvn = list_of_interacting_pixels_label_false.size(); pv < pvn; pv++){
//
//								interacting_p = list_of_interacting_pixels_label_false[pv];
//
//								// these are all false, so
//
//								current_coeff =  int(coefficients.at(interacting_p));
//								if ( number_of_zeros.at(interacting_p) == 1 ){
//									first_derivative_label_false += current_coeff;
//									if (current_coeff < 0){
//										some_negative =  true;
//									}
//
//									if (current_coeff > 0){
//										some_positive = true;
//									}
//								}
//
//
//							}
//						}
//
//
//
//						//						int64_t change_if_flip_to_true;
//						//						int64_t change_if_flip_to_false;
//						bool final_label;
//
//						//flip = (current_value && first_derivative >= 0) || (!current_value && first_derivative < 0);
//						if (list_of_interacting_pixels_label_false.size() + list_of_interacting_pixels_label_true.size() > 1){
//
//							// do we flip, and to what value?   All of the labels need to be the same.
//							//if we flip, it will be to make the value of the cost function go down.
//							flip = (first_derivative_label_true > 0 && !one_cannot_be_changed) || (first_derivative_label_false < 0);
//
//							// if one of the pixels has a rigid label true, we cannot change it to false.
//
//
//							if (flip){
//								//if (flip && some_negative && some_positive){
//								//this should take care of the problem -- why isn't it?
//
//								//								cout << "first d true : " << first_derivative_label_true << endl;
//								//								cout << "first d false : " << first_derivative_label_false << endl;
//								//
//								//								cout << " Number currently true " << list_of_interacting_pixels_label_true.size() << endl;
//								//								cout << " Number currently false " << list_of_interacting_pixels_label_false.size() << endl;
//
//
//
//								// we have to choose
//								if (one_cannot_be_changed){
//									final_label = true;
//								}	else {
//									if ((first_derivative_label_true > 0) && (first_derivative_label_false < 0)){
//										if (fabs(first_derivative_label_true) > fabs(first_derivative_label_false)){
//											final_label = false;
//										}	else {
//											final_label = true;
//										}
//									}	else {
//										if ((first_derivative_label_false < 0)){
//											final_label = true;
//										}	else {
//											final_label = false;
//										}
//									}
//								}
//
//								//cout << "Final label " << final_label << endl;
//								//if (sfis_error/127 < 10000){
//								//cout << "sfis error " << sfis_error/127 << endl;
//								//}
//
//								if (final_label == true){
//									// subtract the f
//									if (fabs(first_derivative_label_false) > sfis_error){
//										cout << "Error!  sfis error is going to go to zero " << endl;
//										cout << "first d false " << first_derivative_label_false << endl;
//										cout << sfis_error << endl;
//										exit(1);
//									}
//
//									sfis_error -= fabs(first_derivative_label_false);
//								}	else {
//									// subtract the f
//									if (fabs(first_derivative_label_true) > sfis_error){
//										cout << "Error!  sfis error is going to go to zero " << endl;
//										cout << "first d false " << first_derivative_label_true << endl;
//										cout << sfis_error << endl;
//										exit(1);
//									}
//
//									sfis_error -= fabs(first_derivative_label_true);
//
//								}
//
//
//								// this covers self as well.
//								if (final_label == true){
//									// change the false-labled ones ....
//									for (int pv = 0, pvn = list_of_interacting_pixels_label_false.size(); pv < pvn; pv++){
//
//										interacting_p = list_of_interacting_pixels_label_false[pv];
//
//										occupied_point_representation.at((interacting_p / (rows*cols))*3 +
//												axis).at(rows*cols*plane_index
//														+ interacting_p % (rows*cols)) = final_label;
//
//
//										if (number_of_zeros.at(interacting_p) == 0){
//											cout << "Error!  N of zeros going less than 0 "<< endl;
//											exit(1);
//										}
//										// changing a false one to a true -- zsubtract an occupied labeled one ...
//										number_of_zeros.at(interacting_p)--;
//
//									}
//
//									number_of_flips += list_of_interacting_pixels_label_false.size();
//								}	else {
//									for (int pv = 0, pvn = list_of_interacting_pixels_label_true.size(); pv < pvn; pv++){
//
//										interacting_p = list_of_interacting_pixels_label_true[pv];
//
//										occupied_point_representation.at((interacting_p / (rows*cols))*3 +
//												axis).at(rows*cols*plane_index
//														+ interacting_p % (rows*cols)) = final_label;
//
//										// changing a true on to a false one
//										number_of_zeros.at(interacting_p)++;
//
//									}
//
//									number_of_flips += list_of_interacting_pixels_label_true.size();
//								}
//							}
//
//						}
//					}
//
//				}
//
//			}
//		}
//	}
//
//}
//


void ReconstructionStructure::ReturnPointGivenCameraPlanePixel(uint cam, int_type_t plane_number, int_type_t within_cam_pixel, vector<double>& X){
	vector<double> cam_plane(4, 0);
	vector<double> local_plane(4, 0);
	vector<double> V(4, 0);
	cam_plane = plane_equations[cam];

	// compute the intersection
	local_plane = cam_plane;
	//local_plane[3] -= -inter_plane_space*double(label_value_for_camera[pixel_index]) + start_plane_distances[cam];
	local_plane[3] -= (start_plane_distances[cam] + inter_plane_space*double(plane_number));


	V[0] = pixel_vectors[cam][3*within_cam_pixel];
	V[1] = pixel_vectors[cam][3*within_cam_pixel + 1];
	V[2] = pixel_vectors[cam][3*within_cam_pixel + 2];


	RayPlaneIntersection(camera_centers[cam], V , local_plane,  X );

}


int_type_t ReconstructionStructure::ReturnIndexFromXYZIndices(int_type_t x, int_type_t y, int_type_t z){

	if (x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z >= number_voxels_grid){
		cout << "ERROR ON size " <<x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z << ", n voxesl " << number_voxels_grid << endl;
		cout << "x , y, z " << x << ", " << y << ",  " << z << endl;
		exit(1);
	}

	return x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z;

}

void ReconstructionStructure::ReturnXYZIndicesFromIndex(int_type_t voxel_index, int_type_t& x, int_type_t& y, int_type_t& z){

	int_type_t temp_index = voxel_index;
	x = temp_index/(number_voxels_per_dim[1]*number_voxels_per_dim[2]);

	temp_index -= x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]);
	y = temp_index/(number_voxels_per_dim[2]);

	temp_index -= y*(number_voxels_per_dim[2]);

	z = temp_index;

	if (x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z != voxel_index){
		cout << "ERROR on vox computation! " << endl << endl;
		exit(1);
	}
}

void ReconstructionStructure::WriteResultImagesGrid(string directory, vector<ImageData*>& cameras, int iteration, string prefix){

	cv::Mat results_image;
	string filename;
	//results_image.copySize(cameras[0]->imcv);
	int pixel_index = 0;
	uchar red, g, b;
	bool special_SfSPM = true; //         ǵĳ    о   true    Ϊһ     ص   Ŀ    ߱     ȷ    


	results_image = cv::Mat(rows, cols, CV_8UC3, cv::Scalar(0,0,0));


	int image_error = 0;
	cout << "number of cameras in Local Min Write " << number_cameras << endl;
	for (int i = 0; i < number_cameras; i++){
		if (i % 10 == 0){
			cout << "Writing " << i << endl;
		}

		// got through and set ....
		uchar* p = &results_image.data[0];

		for (int r = 0; r < 3*rows; r++){
			for (int c = 0; c < cols; c++){
				*p = 0; ++p; // *p = 0; ++p; *p = 0; ++p;
			}
		}

		p = &results_image.data[0];
		//  ɫ  ʾ  ȷ    ɫ  ʾ ۲ ͼ   и     ΪĿ 굫δ   ؽ ģ   б  ֣   ɫ  ʾ ۲ ͼ   и     Ϊ         ؽ ģ   б  ֡ 
		for (int r = 0; r < rows; r++){
			for (int c = 0; c < cols; c++, pixel_index++){

				p = &results_image.data[3*(r*cols + c)];
				b = *p;

				if (coefficients.at(i*rows*cols + r*cols + c) > 0){
					*p = 255;

					if (special_SfSPM){
						*p = 2*(coefficients.at(i*rows*cols + r*cols + c) + 128);
					}
				}


				if ( current_pixel_labeling.at(pixel_index) == false){
					*p = 255;

					if (coefficients.at(i*rows*cols + r*cols + c) < 0){
						image_error++;
					}
				}// ۲ ͼ   и     Ϊ         ؽ ģ   б  ֡ 

				b = *p;

				//++p;
				p = &results_image.data[3*(r*cols + c) + 1];

				if (coefficients.at(i*rows*cols + r*cols + c) > 0){
					*p = 255;

					if (special_SfSPM){
						*p = 2*(coefficients.at(i*rows*cols + r*cols + c) + 128);
					}
				}

				if ( current_pixel_labeling.at(pixel_index) == true){
					*p = 0;
					//					if (cameras[i]->imcv.data[3*(r*cols + c)] > 150){
					//						image_error++;
					//					}
					//
					if (coefficients.at(i*rows*cols + r*cols + c) > 0){
						image_error++;
					}
				}//  ʾ      δ  ģ   б   
				g = *p;

				//++p;

				p = &results_image.data[3*(r*cols + c) + 2];
				//*p = 255;

				if (coefficients.at(i*rows*cols + r*cols + c) > 0){
					*p = 255;
					if (special_SfSPM){
						*p = 2*(coefficients.at(i*rows*cols + r*cols + c) + 128);
					}
				}


				if ( !is_in_undistorted_image.at(i*rows*cols + r*cols + c)){
					*p = 150;
				}

				red = *p;

				++p;

//				if (b == 0 && g == 0 && red == 255){
//					//						cout << "Found bad items ... " << endl;
//					//						cout << "row, col " << r << ", " << c << endl;
//					//						cout << " Current labeling " <<  current_pixel_labeling.at(pixel_index)  << endl << endl;
//					//						cout << "part of original silhouette " << (coefficients.at(i*rows*cols + r*cols + c) > 0 )<< endl;
//				}

			}
		}

		p = &results_image.data[0];


		filename = directory + "iter" + ToString<int>(iteration) + "im" + ToString<int>(i)+ ".png";
		filename = directory + prefix + ToString<int>(iteration) + "im" + ToString<int>(i)+ ".png";
		cv::imwrite(filename.c_str(), results_image);
	}

	cout << "Image error " << image_error << endl;
}

void ReconstructionStructure::SfIS_VH_grid5(vector<int_type_t>& voxels_to_process){

	int_type_t voxel_index = 0;


	if (terms.size() > 0){//terms  ʲô  ˼         أ   
		char ch;

		bool some_negative;
		bool some_positive;

		int_type_t last_term_position = 0;
		int_type_t number_positive = 0;
		int_type_t interact_index  = 0;

#pragma omp parallel for private(voxel_index, last_term_position, some_negative, some_positive, interact_index)
		for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

			if (x_index % 5 == 0){
				cout << "VH grid " << x_index << " out of " << number_voxels_per_dim[0] << endl;
			}
			for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
				for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++){

					voxel_index =  x_index*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;

					if (voxel_index == number_voxels_grid - 1){
						last_term_position = terms.size();//last_term_position  ʲô  ˼  
					}	else {
						last_term_position = start_index_for_terms[voxel_index + 1];//start_index_for_terms  ʲô  ˼    Ӧ     voxel     ص   ʼ   
					}

					some_negative = false;
					some_positive = false;//            ʲô  ˼  
					//				for (int_type_t interact_index = start_index_for_terms[voxel_index]; interact_index < last_term_position
					//				&& (some_negative == false || some_positive == false); interact_index++){
					number_positive = 0;
					for (interact_index = start_index_for_terms[voxel_index];
							interact_index < last_term_position && !(some_negative && some_positive); interact_index++){
						//some_negatibe or some_positive      һ    falseѭ   Ż    
						//				for (int_type_t interact_index = start_index_for_terms[voxel_index];
						//						interact_index < last_term_position
						//				&& (some_negative == false || some_positive == false); interact_index++){

						if (coefficients[terms[interact_index]] < 0){//Ŀ      >0,        <0
							some_negative = true;//ֻҪ  voxel  Ӧ  pixel   б      ؾ ΪTrue
						}	else {
							some_positive = true;//ֻҪ  voxel  Ӧ  pixel    Ŀ     ؾ ΪTrue
							//number_positive++;
						}
					}//Ӧ þ     ֪    voxel  Ӧ    Щ   صķֲ      ȫ    orȫĿ  or   У 

					if (!some_negative && last_term_position != start_index_for_terms[voxel_index])//  ֤          forѭ    ͬʱsome_negative    false.˵  interact_index  Ӧ     ص ȫ  Ŀ      .
					{
						//cout << "Some positive! " << endl;
						configuration_grid[voxel_index] = false;//   ر ǩ    ʼȫΪtrue  ʲô  ˼    false  ʾ  voxel  ռ ݵģ 
						set_by_first_d[voxel_index] = true;//  ʼΪfalse           ʾʲô  ˼  

						for (interact_index = start_index_for_terms[voxel_index];
								interact_index < last_term_position; interact_index++){
#pragma omp critical
							{
								visual_hull_pixels[terms[interact_index]] = true;//ʲô  ˼    
								current_pixel_labeling[terms[interact_index]] = false;//  ʾ  Щpixel  Ӧ     ض   ռ   ˣ   
								number_zeros_grid[terms[interact_index]]++;
							}
						}

					}


					if (some_negative && some_positive){//Ӧ   Ǳ ʾ                 ض   Ŀ       ڲ           silhouette ڲ                    ض   Ŀ       ⲿ          silhouette ⲿ               Ŀ       ڣ         silhouette ཻ            ˼    
						

					}	else {
						set_by_first_d[voxel_index] = true;//ʲô  ˼      ΪTrue    ʾ ܹ ȷ   Ƿ ռ ݣ ΪFalseʱ           жϣ 
					}
				}
			}
		}



	}
	else {//terms().size<=0
		// arr version


		if (terms_arr.size() > 0){
			char ch;

			bool some_negative;
			bool some_positive;

			int_type_t last_term_position = 0;
			int_type_t number_positive = 0;
			int_type_t interact_index  = 0;
			int_type_t arr_index_0, arr_index_1;
			int_type_t term_index, starting_term_position;
			int_type_t* term_index_ptr = 0;


#pragma omp parallel for private(voxel_index, last_term_position, some_negative, some_positive, interact_index, arr_index_0, arr_index_1, starting_term_position, term_index, term_index_ptr)
			for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

				if (x_index % 5 == 0){
					cout << "VH grid " << x_index << " out of " << number_voxels_per_dim[0] << endl;
				}
				for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
					for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++){

						voxel_index =  x_index*number_voxels_per_dim[1]*number_voxels_per_dim[2] + y_index*number_voxels_per_dim[2] + z_index;

						starting_term_position = start_index_for_terms_arr[voxel_index];
						last_term_position = start_index_for_terms_arr[voxel_index + 1];

						some_negative = false;
						some_positive = false;
						//				for (int_type_t interact_index = start_index_for_terms[voxel_index]; interact_index < last_term_position
						//				&& (some_negative == false || some_positive == false); interact_index++){
						number_positive = 0;
						arr_index_0 = starting_term_position/term_increment;
						arr_index_1 = starting_term_position % term_increment;
						//term_index_ptr = &terms_arr[arr_index_0][arr_index_1];
						for (interact_index = starting_term_position;
								interact_index < last_term_position && !(some_negative && some_positive); interact_index++){

							term_index = terms_arr[arr_index_0][arr_index_1];
							//term_index = *term_index_ptr;

							if (coefficients[term_index] < 0){
								some_negative = true;
							}	else {
								some_positive = true;
								//number_positive++;
							}

							arr_index_1++;
							//++term_index_ptr;
							if (arr_index_1 == term_increment){
								arr_index_0++;
								arr_index_1 = 0;
								//term_index_ptr = &terms_arr[arr_index_0][arr_index_1];
							}
						}

						if (!some_negative && last_term_position != starting_term_position)
						{
							configuration_grid[voxel_index] = false;
							set_by_first_d[voxel_index] = true;

							arr_index_0 = starting_term_position/term_increment;
							arr_index_1 = starting_term_position % term_increment;
							for (interact_index = starting_term_position;
									interact_index < last_term_position; interact_index++){
								term_index = terms_arr[arr_index_0][arr_index_1];
#pragma omp critical
								{
									visual_hull_pixels[term_index] = true;
									current_pixel_labeling[term_index] = false;
									number_zeros_grid[term_index]++;
								}
								arr_index_1++;

								if (arr_index_1 == term_increment){
									arr_index_0++;
									arr_index_1 = 0;
								}
							}

						}

						//number_of_positive_terms.at(voxel_index) = number_positive;
						//				if (!some_positive || !so){
						//					set_by_first_d[voxel_index] = true; // forever empty ....
						//				}

						if (some_negative ^ some_positive){//Ϊʲô              Ϊfalse  true
							set_by_first_d[voxel_index] = true;
						}


						//						if (some_negative && some_positive){
						//							//#pragma omp critical
						//							//					{
						//							//					voxels_to_process.push_back(voxel_index);
						//							//					}
						//						}	else {
						//							set_by_first_d[voxel_index] = true;
						//						}
					}
				}
			}
		}	else {
			// parallel arr version

			char ch;

			bool some_negative;
			bool some_positive;

			int_type_t last_term_position = 0;
			int_type_t number_positive = 0;
			int_type_t interact_index  = 0;
			int_type_t arr_index_0, arr_index_1;
			int_type_t term_index, starting_term_position;
			int_type_t* term_index_ptr = 0;
			int_type_t thread_ID;
			int_type_t counter;
			
#pragma omp parallel private(counter, voxel_index, thread_ID, last_term_position, some_negative, some_positive, interact_index, arr_index_0, arr_index_1, starting_term_position, term_index, term_index_ptr)
			{
				thread_ID = omp_get_thread_num();
				voxel_index = 0;
				counter = 0;

				for (int_type_t x_index = 0; x_index < number_voxels_per_dim[0]; x_index++){

					if (x_index % 20 == 0){
						cout << "VH grid arr" << x_index << " out of " << number_voxels_per_dim[0] << " thread ID " << thread_ID << endl;
					}
					for (int_type_t y_index = 0; y_index < number_voxels_per_dim[1]; y_index++){
						for (int_type_t z_index = 0; z_index < number_voxels_per_dim[2]; z_index++, voxel_index++){

							if (voxel_index % number_divisions == thread_ID){


								starting_term_position = start_indices_terms_per_div[thread_ID][counter];
								last_term_position = start_indices_terms_per_div[thread_ID][counter + 1];
								counter++;

								if (set_by_first_d[voxel_index] == false){
									some_negative = false;
									some_positive = false;
									//				for (int_type_t interact_index = start_index_for_terms[voxel_index]; interact_index < last_term_position
									//				&& (some_negative == false || some_positive == false); interact_index++){
									number_positive = 0;
									arr_index_0 = starting_term_position/term_increment;
									arr_index_1 = starting_term_position % term_increment;
									//term_index_ptr = &terms_arr[arr_index_0][arr_index_1];
									for (interact_index = starting_term_position;
											interact_index < last_term_position && !(some_negative && some_positive); interact_index++){

										term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
										//term_index = *term_index_ptr;

										if (coefficients[term_index] < 0){
											some_negative = true;
										}	else {
											some_positive = true;
											//number_positive++;
										}

										arr_index_1++;
										//++term_index_ptr;
										if (arr_index_1 == term_increment){
											arr_index_0++;
											arr_index_1 = 0;
											//term_index_ptr = &terms_arr[arr_index_0][arr_index_1];
										}
									}

									if (!some_negative && last_term_position != starting_term_position)
									{
										configuration_grid[voxel_index] = false;
										set_by_first_d[voxel_index] = true;

										arr_index_0 = starting_term_position/term_increment;
										arr_index_1 = starting_term_position % term_increment;
										for (interact_index = starting_term_position;
												interact_index < last_term_position; interact_index++){
											term_index = terms_arr_per_div[thread_ID][arr_index_0][arr_index_1];
#pragma omp critical
											{
												number_zeros_grid[term_index]++;
											}
											arr_index_1++;

											if (arr_index_1 == term_increment){
												arr_index_0++;
												arr_index_1 = 0;
											}
										}

									}

									//number_of_positive_terms.at(voxel_index) = number_positive;
									//				if (!some_positive || !so){
									//					set_by_first_d[voxel_index] = true; // forever empty ....
									//				}

									if (some_negative ^ some_positive){
										set_by_first_d[voxel_index] = true;
									}
								}
							}
						}
					}
				}
			}


			for (int_type_t p = 0; p < number_pixels; p++){
				if (number_zeros_grid[p] > 0){
					visual_hull_pixels[p] = true;
					current_pixel_labeling[p] = false; // technically, can do this later ...
				}
			}
		}
	}




	if (terms_arr_per_div.size() == 0){

		for (voxel_index = 0; voxel_index < number_voxels_grid; voxel_index++){
			if (!set_by_first_d[voxel_index]){
				voxels_to_process.push_back(voxel_index);
			}
		}
	}	else {
		cout << "Number divisions" << number_divisions << endl;
		for (voxel_index = 0; voxel_index < number_voxels_grid; voxel_index++){
			if (!set_by_first_d[voxel_index]){
				// will be sending back the "counter" value, do get the voxel value out of that, mult by number divisions and add thread ID
				multi_dim_voxels_to_process[voxel_index % number_divisions][voxel_index/number_divisions] = true;
			}
		}
	}


	// assess error:

	sfis_error = function_constant;
	for (int_type_t p = 0; p < number_pixels; p++){
		if (current_pixel_labeling[p] == true){
			sfis_error += coefficients[p];
		}
	}
}
