/*
 * SfIS.hpp
 *
 *  Created on: Oct 20, 2020
 *      Author: Jinglong
 */

#ifndef SFIS_HPP_
#define SFIS_HPP_

#include "ImageData.hpp"
#include "Includes.hpp"
#include "ReconstructionStructure.hpp"

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

int LocalMinSearch6(int demo_number);

void DoSearch(ReconstructionStructure* RS_split, ofstream& out, bool* terms_changed_this_round );

//int LocalMinSearch(int demo_number, int type_method);
//
//int LocalMinSearchSparse(int demo_number, int type_method);
//
//int StereoOrSfIS(int demo_number, int type_method);
//
void WriteBoundingBox(string outfile, vector< vector<double> >& bounding_box);
//
//
//
//void WriteCaliFile(vector<ImageData*>& cameras, string dst_dir, double downsample);






#endif /* SFIS_HPP_ */
