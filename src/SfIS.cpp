/*
 * NewRep.cpp
 *
 *  Created on: Aug 17, 2020
 *      Author: Jinglong
 */

#include "SfIS.hpp"
#include "Includes.hpp"
#include "Helper0.hpp"
#include "Set.hpp"
#include "ReconstructionStructure.hpp"
#include "DirectoryFunctions.hpp"


#ifdef _MSC_VER
#include <time.h>
#include <windows.h>
int gettimeofday(struct timeval *tp, void *tzp);
#else
#include <sys/time.h>
#endif


#ifdef _MSC_VER
int gettimeofday(struct timeval *tp, void *tzp)
{
	time_t clock;
	struct tm tm;
	SYSTEMTIME wtm;
	GetLocalTime(&wtm);
	tm.tm_year = wtm.wYear - 1900;
	tm.tm_mon = wtm.wMonth - 1;
	tm.tm_mday = wtm.wDay;
	tm.tm_hour = wtm.wHour;
	tm.tm_min = wtm.wMinute;
	tm.tm_sec = wtm.wSecond;
	tm.tm_isdst = -1;
	clock = mktime(&tm);
	tp->tv_sec = clock;
	tp->tv_usec = wtm.wMilliseconds * 1000;
	return (0);
}
#endif

int LocalMinSearch6(int demo_number){

	char ch;
	bool verbose  = true;
	string source_file;

	cout << "This code is written assuming that the number of terms and the number of variables does not exceed 2^32 = 4,294,967,296 " << endl;

	float division = 100.0;

	string write_directory ="";
	double downsample = 1.0;

	string b_costs; //    
	string i_costs; //    
	string smooth;  string div;
	string short_dir_name;
	string base_dir;
	string div_dir;

	vector<double> pA(3, 0);
	vector<double> pB(3, 1000);
	// fill in based on BB
	//Point_3 dummy_point;

	ifstream in;
	int master_number = 0;
	vector<vector<double> > BB;
	string filename;

	ReconstructionStructure* RS_split = 0;
	ReconstructionStructure* RS_parent = 0;
	//GraphClass GC;
	double starttime = 0;
	double endtime = 0;
	double multiplier = 1;

	vector<int_type_t> voxels_to_process;
	vector<int_type_t> pixels_to_process;
	int_type_t number_pinned;
	int_type_t test_error;

	vector<ReconstructionStructure*> RSes;


	/////////////////////////***********************************************************
	////////// Preliminary stuff ///////////////////////////////////////////////////////



	int number_splits = SetItemsForNewRep(demo_number, source_file,
			base_dir, write_directory, downsample, division, pA, pB, multiplier);

	BB.push_back(pA);
	BB.push_back(pB);



	write_directory = write_directory + "6LS" + ToString<float>(division) + "_DS" + ToString<float>(downsample) + "_SP" + ToString<float>(number_splits) + "/";


	//string command= "rm -r " + write_directory;
	int r = 0;// system (command.c_str());
	//mkdir (write_directory.c_str(), S_IRWXU) ;

	string shuangyinhao = "\"";
	string object_dir = write_directory + "object_files";
	string command_obj_dir = "mkdir " + shuangyinhao + object_dir+ shuangyinhao;
	r = system(command_obj_dir.c_str());

	object_dir = write_directory + "vh_files";
	command_obj_dir= "mkdir " + shuangyinhao+object_dir+ shuangyinhao;
	r = system(command_obj_dir.c_str());

	object_dir = write_directory + "smoothed_files";
	command_obj_dir= "mkdir " + shuangyinhao + object_dir + shuangyinhao;
	r = system(command_obj_dir.c_str());

	object_dir = write_directory + "connected_component";
	command_obj_dir= "mkdir " + shuangyinhao + object_dir + shuangyinhao;
	r = system(command_obj_dir.c_str());

	string write_dir = write_directory;
	string current_write_file = write_directory + "trace.txt";
	//string short_dir_name = short_dir_name;

	//BSP* bsp = 0;

	ReconstructionStructure R_model;

	std::ofstream out;
	out.open(current_write_file.c_str());


	cout << "Source file  :" << source_file << endl;

	///////*************************************************///
	///******* LOAD THE CAMERAS ****************************///


	ImageData* id0 = 0;
	vector<ImageData*> cameras;

	in.open(source_file.c_str());
	in >> master_number;
	cout << "master Number " << master_number << endl;

	vector<int> downsampling;


	struct timeval tvA;
	struct timeval tvB;

	struct timeval tvC;
	struct timeval tvD;

	struct timeval tvE;
	struct timeval tvF;


	filename = write_directory + "bounding_box.ply";
	WriteBoundingBox(filename, BB);

	gettimeofday( & tvA, 0);

	LoadImageDataVector(cameras, source_file, multiplier, downsample, true);

	for (int i = 0; i < int(cameras.size()); i++){
		cameras[i]->index  = i;
	}

	cout << "After load " << endl;


	// don't need the artifical cameras

	out << "***********************" << endl;
	out << "Bounding box " << endl;
	for (int i  = 0; i < 3; i++){
		out << BB[0][i] << " " ;
	}

	out << endl;
	for (int i  = 0; i < 3; i++){
		out << BB[1][i] << " " ;
	}
	out << endl;


	cout << "After add " << endl;
	cout << "Number of cameras "<< cameras.size() << endl;

	//cout << "Number divisions for this computer " << RSnumber_divisions << endl;
	//cout << "Term increment " << RS.term_increment << endl;
	int_type_t number_to_test; 

	// the real value in this is allowing a downsampling run to go first ....
	int original_divider = 2; int original_distance_search = 8;
	int divider = original_divider; int distance_search = original_distance_search;
	int_type_t last_error;  int_type_t max_rounds_one_level = 1;
	int_type_t counter_at_one_level = 0;
	int_type_t count_splits = 0;

	gettimeofday( & tvA, 0);
	int_type_t s = 0;

	/////////////////////////////////////// Compute preliminary reconstruction ///////////////////////////////////////////

	count_splits = 0;

	// first round -- do refular reconstruction ....
	RS_split = new ReconstructionStructure();
	RSes.push_back(RS_split);

	// now work with the RS-child pointer ....
	RS_split->CreateGrid(BB, cameras, division, downsample);
	cout << "After creating grid ...." << endl;


	out << "s  = " << s << endl;
	//exit(1);

	gettimeofday( & tvC, 0);
	RS_split->ComputeProjectionInformationSpheresArrRepresentationParallel();
	gettimeofday( & tvD, 0);
	starttime = double(tvC.tv_usec)/1000000.0;
	endtime = double(tvD.tv_sec - tvC.tv_sec) + double(tvD.tv_usec)/1000000.0;
	out << "Time for back-projection ... " << endtime - starttime << endl;

	// can just clone this for later ones ....
	RS_split->ConstructSfISCostFunction(cameras, downsample);



	out << "Error without any label " << FormatWithCommas<int_type_t>(RS_split->sfis_error/127) << endl;
	bool* terms_changed_this_round = new bool[RS_split->coefficients.size()];// ߴ Ϊnumber_pixel  bool ͱ   

	gettimeofday( & tvC, 0);
	RS_split->PrintFeaturesToFile(out);

	cout << "Before VH " << endl;
	RS_split->SfIS_VH_grid5(voxels_to_process);

	//	RS_split->WriteResultImagesGrid(write_directory, cameras, -1);
	//	RS_split->GenerateAndWriteSurfaceInPlyFormat(write_directory, -1);

	out << "Error from visual hull " << FormatWithCommas<int_type_t>(RS_split->sfis_error/127) << endl;
	DoSearch(RS_split, out, terms_changed_this_round);

	gettimeofday( & tvD, 0);
	starttime = double(tvC.tv_usec)/1000000.0;
	endtime = double(tvD.tv_sec - tvC.tv_sec) + double(tvD.tv_usec)/1000000.0;
	out << "Time for minimization ... " << endtime - starttime << endl;


	last_error = RS_split->sfis_error;

	RSes[s]->DeallocBigArrays();


	s++;
	division /= 2.0;
	//////////////////////////////////// Splitting //////////////////////////////////////////////
	out << "Starting a split, division is " << division << endl;

	// so idea for the ICRA paper is to do a calibration refinement before the final division.

	// later -- fix up the backprojection of splits so that several new arrays are made each time .. then
	// for minimization, we just have to look up which channel has the data.

	int_type_t total_number_to_test;
	while (count_splits < number_splits){
		total_number_to_test = 0;

		gettimeofday( & tvE, 0);

		RS_split = new ReconstructionStructure();
		RSes.push_back(RS_split);//          ڸ ʲô    

		// now work with the RS-child pointer ....
		//RS_split->CreateGrid(BB, cameras, division, downsample);
		RS_split->CreateGridUsingAncestor(RSes[0], BB, cameras, division, downsample);
		cout << "After creating grid ...." << endl;
		cout << "TIle number " << RS_split->tile_number.size() << endl;//tile_numberָʲô  

		out <<"=========================================================" << endl;
		out << "Voxel size " << division << endl;
		out << "s  = " << s << endl;

		out << "distance searched is " << max(distance_search/s, int_type_t(1)) << endl;
		// then flip back after projection multi dim to process for those with firs d = true.
		number_to_test = RSes[s - 1]->PropogateLabelsAccordingToChild(*RS_split, divider, 2, max(distance_search/s, int_type_t(1)));
		cout << "TIle number2 " << RS_split->tile_number.size() << endl;

		out << "Number to test " << FormatWithCommas<int_type_t>(number_to_test) << " out of "
				<< FormatWithCommas<int_type_t>(RS_split->number_voxels_grid)
				<< "  proportion tested " << double(number_to_test)/double(RS_split->number_voxels_grid) << endl;

		total_number_to_test = number_to_test;
		if (number_to_test == 0){
			cout << "Number to test is zero ... " << endl;
		}
		gettimeofday( & tvC, 0);

		RS_split->ComputeProjectionInformationSpheresArrRepresentationParallel();
		gettimeofday( & tvD, 0);
		starttime = double(tvC.tv_usec)/1000000.0;
		endtime = double(tvD.tv_sec - tvC.tv_sec) + double(tvD.tv_usec)/1000000.0;
		out << "Time for back-projection ... " << endtime - starttime << endl;

		RS_split->ConstructSfISCostFunctionUsingAncestor(RSes[0]);

		RS_split->UpdateNumberZeros();

		RS_split->MarkToTestAsFalseFirstDsOnly();

		out << "Error without any label change " << FormatWithCommas<int_type_t>(RS_split->sfis_error/127) << endl;

		RS_split->PrintFeaturesToFile(out);

		DoSearch(RS_split, out, terms_changed_this_round);

		last_error = RS_split->sfis_error;

		int_type_t loop_counter = 1;
		bool some_change = true;
		//some_change = false;
		while (some_change){

			out << "Working on extending the current split .... " << division << " iterator " << loop_counter << endl;

			number_to_test = RSes[s]->SameLevelExpansion(distance_search*2);

			out << "Number to test " << FormatWithCommas<int_type_t>(number_to_test) << " out of "
					<< FormatWithCommas<int_type_t>(RS_split->number_voxels_grid)
					<< "  new proportion tested " << double(number_to_test)/double(RS_split->number_voxels_grid) << endl;

			total_number_to_test += number_to_test;
			out << "Total number to test " << FormatWithCommas<int_type_t>(total_number_to_test) << " out of "
					<< FormatWithCommas<int_type_t>(RS_split->number_voxels_grid)
					<< "  proportion tested " << double(total_number_to_test)/double(RS_split->number_voxels_grid) << endl;



			if (number_to_test > 0){
				//cout << "Before " << endl; cin >> ch;
				gettimeofday( & tvC, 0);

				RS_split->ComputeProjectionInformationSpheresArrRepresentationParallelForOnlyThoseMarked();
				gettimeofday( & tvD, 0);
				starttime = double(tvC.tv_usec)/1000000.0;
				endtime = double(tvD.tv_sec - tvC.tv_sec) + double(tvD.tv_usec)/1000000.0;
				out << "Time for back-projection ... " << endtime - starttime << endl;

				//cout << "After " << endl;
				//cin >> ch;

				DoSearch(RS_split, out, terms_changed_this_round);

				//if (loop_counter > 1){
				if (last_error == RS_split->sfis_error){
					some_change = false;
				}

				loop_counter++;

				last_error = RS_split->sfis_error;
			}	else {
				some_change = false;
			}
		}

		gettimeofday( & tvF, 0);
		starttime = double(tvE.tv_usec)/1000000.0;
		endtime = double(tvF.tv_sec - tvE.tv_sec) + double(tvF.tv_usec)/1000000.0;
		out << "Time for minimization of a split... " << endtime - starttime << endl;

		RSes[s]->DeallocBigArrays();

		count_splits++; s++;
		division /= 2.0;
	}

	delete [] terms_changed_this_round;
	gettimeofday( & tvB, 0);

	starttime = double(tvA.tv_usec)/1000000.0;
	endtime = double(tvB.tv_sec - tvA.tv_sec) + double(tvB.tv_usec)/1000000.0;
	out << endl << "Time for all of the operations ... " << endtime - starttime << endl;

	cout << "Writing ... " << endl;
	for (int_type_t i = 0; i < RSes.size(); i++){
		cout << "Writing  " << i << endl;
		RSes[i]->WriteResultImagesGrid(write_directory, cameras, i*100);
		RSes[i]->GenerateAndWriteSurfaceInPlyFormat(write_directory, i*100);

		string proposed_file = write_directory + "/smoothed_files/" + ToString<int>(i*100) + ".txt";

		ofstream out_model;
		out_model.open(proposed_file.c_str());

		for (int_type_t j = 0; j < RSes[i]->number_voxels_grid; j++){
			if (RSes[i]->configuration_grid[j] == 0){
				out_model << j << " ";
			}
		}

		out_model << endl << -1 << endl;

		out_model.close();
	}



	out.close();
	cout << "After delete " << endl;

	return 0;
}

void DoSearch(ReconstructionStructure* RS_split, ofstream& out, bool* terms_changed_this_round ){
	int_type_t number_to_test = 0;

	number_to_test = 0;
	for (int_type_t i = 0; i < RS_split->number_divisions; i++){
		for (int_type_t j = 0; j < RS_split->number_voxels_grid/RS_split->number_divisions  + 1; j++){
			number_to_test += RS_split->multi_dim_voxels_to_process[i][j];
		}
	}


	out << "Number voxels to process in next iter " << FormatWithCommas<int_type_t>(number_to_test) << endl;
	cout << "Number voxels to process in next iter " << FormatWithCommas<int_type_t>(number_to_test) << endl;


	int_type_t number_changed = 500;
	//bool* terms_changed_this_round = new bool[RS_split->coefficients.size()];
	//if (s == 0)
	{

		// should only have
		for (int iteration  = 1; iteration < 50 && number_changed > 2; iteration++){//Ϊʲô  ô      
			if (RS_split->start_indices_terms_per_div_tiles.size() == 0){
				number_changed = RS_split->LocalMinSearchParallel6Arr(terms_changed_this_round, false);
			}	else {
				//number_changed = RS_split->LocalMinSearchParallel6ArrTiles(terms_changed_this_round, false);
			}

			number_to_test = 0;
			for (int_type_t i = 0; i < RS_split->number_divisions; i++){
				for (int_type_t j = 0; j < RS_split->number_voxels_grid/RS_split->number_divisions  + 1; j++){
					number_to_test += RS_split->multi_dim_voxels_to_process[i][j];
				}
			}
			out << iteration << " Error after testing one round:   " << FormatWithCommas<int_type_t>(RS_split->sfis_error/127) << "   number changed "
					<< number_changed << " next round " << FormatWithCommas<int_type_t>(number_to_test) << endl;
			cout << iteration << " Error after testing one round:   " << FormatWithCommas<int_type_t>(number_to_test) << endl;
		}
		out << endl << endl;
	}


	//cout << "This version " << (RS_split->start_indices_terms_per_div_tiles.size() == 0) << endl;
	//char ch; cin >> ch;
}



void WriteBoundingBox(string outfile, vector< vector<double> >& bounding_box){

	// each vertex needs a color ....

	cout << "Writing to " << outfile << endl;
	std::ofstream out;
	out.open(outfile.c_str());

	out << "ply" << endl;
	out << "format ascii 1.0" << endl;
	out << "element vertex " << 8 << endl;
	out << "property float x" << endl;
	out << "property float y" << endl;
	out << "property float z" << endl;
	out << "property uchar red" << endl;
	out << "property uchar green" << endl;
	out << "property uchar blue" << endl;
	out << "property uchar alpha" << endl;
	out << "element face " << 6 << endl;
	out << "property list uchar int vertex_indices"<< endl;
	//	out << "element edge " << edge_indices.size() << endl;
	//	out << "property int vertex1" << endl;
	//	out << "property int vertex2" << endl;
	//	out << "property uchar red" << endl;
	//	out << "property uchar green" << endl;
	//	out << "property uchar blue" << endl;
	//	out << "property uchar alpha" << endl;
	out << "end_header" << endl;

	unsigned int zero = 0;

	// P0
	out << bounding_box[0][0] << " " << bounding_box[0][1] << " " << bounding_box[0][2] << " 255 150 150 100" << endl;

	// P1
	out << bounding_box[0][0] << " " << bounding_box[0][1] << " " << bounding_box[1][2] << " 255 150 150 100" << endl;

	//P2
	out << bounding_box[0][0] << " " << bounding_box[1][1] << " " << bounding_box[0][2] << " 255 150 150 100" << endl;

	//P3
	out << bounding_box[1][0] << " " << bounding_box[0][1] << " " << bounding_box[0][2] << " 255 150 150 100" << endl;

	//P4
	out << bounding_box[0][0] << " " << bounding_box[1][1] << " " << bounding_box[1][2] << " 255 255 150 100" << endl;

	// P5
	out << bounding_box[1][0] << " " << bounding_box[0][1] << " " << bounding_box[1][2] << " 255 255 150 100" << endl;

	// P6
	out << bounding_box[1][0] << " " << bounding_box[1][1] << " " << bounding_box[0][2] << " 150 150 255 100" << endl;

	//P7
	out << bounding_box[1][0] << " " << bounding_box[1][1] << " " << bounding_box[1][2] << " 150 150 255 100" << endl;


	// Faces
	out << "4 6 2 4 7 " << endl;

	out << "4 5 3 6 7 " << endl;

	out << "4 1 5 3 0 " << endl;

	out << "4 4 1 0 2 " << endl;

	out << "4 6 2 0 3 " << endl;

	out << "4 1 4 7 5 " << endl;


	out.close();
}
