/*
 * Set.cpp
 *
 *  Created on: July 19, 2020
 *      Author: Jinglong
 */

#include "Set.hpp"

int SetItemsForNewRep(int demo_number, string& source_file,
		string& base_dir, string& write_directory, double& downsample, float& division,
		vector<double>& pA, vector<double>& pB, double& multiplier){

	multiplier = 1;//平移矩阵T的单位和体素分辨率单位需要统一，不一致时改变该参数

	int number_splits = 0;

	//double PFA_pix, PM_pix;
	switch (demo_number){
	case 0: {
		base_dir =  "./Sphere";

		source_file = base_dir + "/cali.txt";

		write_directory = base_dir + "/Experiment/";

		downsample = 1;

		division = 8;
		//division = 0.0200;

		pA = { -70, -70, -70 };
		pB = { 70, 70, 70 };

		number_splits = 4;

		/*if (1)
		{
			source_file = base_dir + "/cali_asteroid.txt";

			write_directory = base_dir + "/Experiment/";

			downsample = 1.0;
			division = 14;

			pA = { -70, -70, -70 };
			pB = { 70, 70, 70 };

			number_splits = 4;
		}
		else
		{
			source_file = base_dir + "/cali.txt";

			write_directory = base_dir + "/Experiment/";

			downsample = 1.0;
			division = 20.0;

			pA = { -140, -120, -140 };
			pB = { 140, 120, 140 };

			number_splits = 1;
		}*/
		
	}	break;

	}

	return number_splits;
}

