/*
 * Set.hpp
 *
 *  Created on: July 19, 2020
 *      Author: Jinglong
 */

#ifndef SET_HPP_
#define SET_HPP_
#include "ImageData.hpp"
#include "Includes.hpp"

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

int SetItemsForNewRep(int demo_number, string& source_file,
		string& base_dir, string& write_directory, double& downsample, float& division,
		vector<double>& pA, vector<double>& pB, double& multiplier);


#endif /* SET_HPP_ */
