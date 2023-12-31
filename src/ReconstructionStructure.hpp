/*
 * ReconstructionStructure.hpp
 *
 *  Created on: July 22, 2020
 *      Author: Jinglong
 */

#ifndef RECONSTRUCTIONSTRUCTURE_HPP_
#define RECONSTRUCTIONSTRUCTURE_HPP_

#include "ImageData.hpp"
#include "Includes.hpp"

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

#include <iomanip>
#include <locale>

class ReconstructionStructure{
public:

	vector<bool> configuration_grid;
	vector<bool> set_by_first_d;
	vector<int_type_t> number_zeros_grid;

	vector<bool> visual_hull_pixels;
	vector<bool> current_pixel_labeling;

	vector<int_type_t> terms;
	vector<int_type_t> start_index_for_terms;
	vector<offset_int_type_t> offsets_for_terms;

	int_type_t term_increment;
	vector< int_type_t* > terms_arr;
	int_type_t* start_index_for_terms_arr;
	vector<int16_t*> tile_number;

	int_type_t number_divisions;
	vector< int_type_t* > start_indices_terms_per_div;
	vector< vector< int_type_t* > > terms_arr_per_div;
	vector< bool* > multi_dim_voxels_to_process;
	vector< int_type_t > voxel_ordering;

	vector< vector< int_type_t* > >  start_indices_terms_per_div_tiles;
	vector< vector< vector< int_type_t* > > > terms_arr_per_div_tiles;
	vector< int_type_t* > number_terms_per_div;
	vector< int_type_t** > terms_ptr_per_div;


	vector<vector<int_type_t> > terms_by_image;
	vector<vector<int_type_t> > start_index_for_terms_by_image;
	vector<vector<offset_int_type_t> > offsets_for_terms_by_image;

	vector<double> initial_offset;
	double inter_voxel_distance;
	vector<int_type_t> number_voxels_per_dim;
	int_type_t number_voxels_grid;


	vector<int32_t> camera_label_value;
	vector<int64_t> label_value_for_camera;
	vector< vector<double> > pixel_vectors;
	vector<bool> is_in_undistorted_image;
	vector<int_type_t> random_pixels;
	vector<int_type_t> random_planes;
	vector<int8_t> coefficients;
	vector< vector< vector<double> > > projection_matrices;
	vector< vector< vector<double> > > inverse_projection_matrices;
	vector< vector<double> > camera_centers;

	// new
	vector<double> start_plane_distances; // this is label '0'
	vector<double> end_plane_distances;
	vector<int> number_of_intervals_per_cam; // (end_plane - start_plane)/inter-plane_space

	vector< vector<double> > plane_equations;

	vector< vector<double> > BB;
	int_type_t number_pixels;
	int_type_t number_voxels;
	int_type_t number_cameras;
	double inter_plane_space;
	int rows;
	int cols;
	int64_t function_constant;
	uint64_t sfis_error;

	ReconstructionStructure();

	~ReconstructionStructure();

	void ConstructSfISCostFunction(vector<ImageData*>& cameras);

	void WriteResults(string directory, vector<ImageData*>& cameras, int iteration);


	///*****************************************
	void ReturnPointGivenCameraPlanePixel(unsigned int cam, int_type_t plane_number, int_type_t within_cam_pixel, vector<double>& X);

	void CreateGrid(vector< vector<double> >& boundingbox, vector<ImageData*>& cameras, double division, int downsample);

	void PrintFeaturesToFile(std::ofstream& out);

	void SfIS_VH_grid5(vector<int_type_t>& voxels_to_process);

	void WriteResultImagesGrid(string directory, vector<ImageData*>& cameras, int iteration, string prefix = "iter");

	int_type_t ReturnIndexFromXYZIndices(int_type_t x, int_type_t y, int_type_t z);

	void QuickFirstDs(vector<int_type_t>& voxels_to_process);

	void GenerateAndWriteSurfaceInPlyFormat(string outdir, int iteration, string prefix = "iter", int* cs=NULL);

	void GenerateAndWriteSurfaceInPlyFormatSmooth(string outdir, int iteration, int number_times);

	void GenerateSubdividedStructure(ReconstructionStructure& SRS);

	void FillInOneIteration();

	void WritePlyFile(string outfile,
			vector<int_type_t>& subdims,
			vector<double>& points,
			vector<int_type_t>& faces,
			vector<int_type_t>& color);

	void ReturnXYZIndicesFromIndex(int_type_t voxel_index, int_type_t& x, int_type_t& y, int_type_t& z);

	void ReturnSetOfFaces(int_type_t x_index, int_type_t y_index, int_type_t z_index, vector<double>& points, vector<int_type_t>& faces);

	void ComputeProjectionInformationSpheresArrRepresentationParallel();

	int_type_t PropogateLabelsAccordingToChild(ReconstructionStructure& RS_child, int division, int inner_distance, int distance);

	void ConstructSfISCostFunction(vector<ImageData*>& cameras, int downsample);

	void UpdateNumberZeros();

	void CheckProjectionForErrors();

	int_type_t MarkShellAsToProcess();

	void MarkToTestAsFalseFirstDsOnly();

	void ResetProjectionMatricesAndNumberZeros(vector<ImageData*>& cameras, int downsample);

	int_type_t PropogateLabelsAndMark(int distance);

	void DeallocBigArrays();

	int_type_t LocalMinSearchParallel6Arr(bool* terms_changed_this_round, bool first);

	int_type_t SameLevelExpansion(int distance);

	void ComputeProjectionInformationSpheresArrRepresentationParallelForOnlyThoseMarked();

	void ConstructSfISCostFunctionUsingAncestor(ReconstructionStructure* ancestor);

	void ComputeProjectionInformationSpheresTiles();

	void CreateGridUsingAncestor(ReconstructionStructure* RS, vector< vector<double> >& boundingbox, vector<ImageData*>& cameras, double division, int downsample);

	int_type_t LocalMinSearchParallel6ArrTiles(bool* terms_changed_this_round, bool first);

	void CreateGridMinimal(vector< vector<double> >& boundingbox, double division);

	int_type_t CompareWithModel(ReconstructionStructure* RS_exp, double& tn, double& tp, double& fn, double& fp );

	int_type_t QuicklyAlterZeros(bool* terms_changed_this_round, bool first);

	int_type_t NumberZeros();

	int_type_t CleanIsolatedVoxels();

	void WriteSurfaceInPlyFormat(vector<double>& points, vector<int_type_t>& faces, string outdir, int iteration, string prefix);

	void ReturnReconstructionFaces(vector<double>& points, vector<int_type_t>& faces, vector<int_type_t>& point_face_rep);

	void CreateGridUsingAncestorWithProjection(ReconstructionStructure* RS, vector< vector<double> >& boundingbox, double division, int downsample);

	void GenerateAndWriteSurfaceInPlyFormat(string outdir, int iteration, string prefix, vector<int>& cs);

	void AddProjectionMatrices(vector<ImageData*>& cameras);

};



template<class T>
std::string FormatWithCommas(T value)
{
	int_type_t uvalue = value;
	bool negative = false;
	if (value < 0){
		negative = true;
		uvalue = -value;
	}


	string s;
	  int cnt = 0;
	  do
	  {
	    s.insert(0, 1, char('0' + uvalue % 10));
	    uvalue /= 10;
	    if (++cnt == 3 && uvalue)
	    {
	      s.insert(0, 1, ',');
	      cnt = 0;
	    }
	  } while (uvalue);

	  if (negative){
		  s = "-" + s;
	  }
	  return s;
}

void PrintVector(vector<double>& p);

bool IsLegalPoint(ReconstructionStructure& RS, int source_cam, int within_cam_index, pair<int, int>cam_label_pair);

#endif /* RECONSTRUCTIONSTRUCTURE_HPP_ */
