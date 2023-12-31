/*
 * main.cpp
 * This code is a companion to Jinglong's master dissertation
*/

#include "main.hpp"
#include <sys/stat.h>
#include <time.h>
#include "SfIS.hpp"

using std::list;
using std::cout;
using std::endl;
using std::string;

#ifndef _WIN64
//#ifdef _DEBUG
//#       pragma comment(lib, "D:/Dissertation/SFS/opencv/opencv_world347d.lib")
//#else
//#       pragma comment(lib, "D:/Dissertation/SFS/opencv/opencv_world347.lib")
//#endif
#else
#ifdef _DEBUG
#       pragma comment(lib, "D:/Dissertation/SFS/opencv/opencv_world347d.lib")
#else
#       pragma comment(lib, "D:/Dissertation/SFS/opencv/opencv_world347.lib")
#endif
#endif

int main(int argc, char **argv)
{
	char ch;
	int n = 0;

	if (argc > 1){
		n = atoi(argv[1]);
	}	else {
		cout << "No dataset number selected, assuming selection is 0" << endl;
	}


	LocalMinSearch6(n);

	getchar();

	return(0);
}

/*int main()
{
	

	LocalMinSearch6(0);

	getchar();

	return(0);
}*/
