/*
 * ImageData.cpp
 *
 *  Created on: Aug 6, 2019
 *      Author: Jinglong
 */
#include "ImageData.hpp"
//#include "StandardImageFunctions.hpp"
//#include "Undistort.hpp"





void DeleteImageData(ImageData* id){
	cout << id << endl;
	delete id;
}


ImageData::~ImageData(){


}

ImageData::ImageData(){}

void LoadImageDataVector(vector<ImageData*>& id_vector, string filename, double multiplier, float downsample_factor, bool image_load){

	cv::Mat view, rview, map1, map2;
	cv::Mat temp_im;
	ifstream filein;
	filein.open(filename.c_str());
	string name;
	int number_elements = 0;
	while (filein >> name){
		number_elements++;
	}
	filein.close();

	filein.open(filename.c_str());
	//continuous_im = 0;
	int master_num;//有几幅图像
	filein >> master_num;

	cv::Mat master_cameraMatrix = cv::Mat::eye(3, 3, CV_64F);//存储第0张图像的相机内参数
	cv::Mat master_distCoeffs = cv::Mat::zeros(8, 1, CV_64F);//存储第0张图像的相机即便参数

	cv::Mat cameraMatrix = cv::Mat::eye(3, 3, CV_64F);//存储相机内参数

	cv::Mat BigcameraMatrix = cv::Mat::eye(3, 3, CV_64F);
	cv::Mat distCoeffs = cv::Mat::zeros(8, 1, CV_64F);
	cv::Mat Rt = cv::Mat::zeros(3, 4, CV_64F); //存储[R|T]矩阵
	cv::Mat P;
	cv::Mat C = cv::Mat::zeros(3, 1, CV_64F);

	// can we compute how many elements are in this file ...?

	int number_items_per_image = (number_elements - 1)/master_num;//每幅图像有几个参数
	int number_distortion_elements = number_items_per_image - 1 - 9 - 9- 3;//相机畸变参数有几个
	cout << "number distortion elements " << number_distortion_elements << endl;
	//char ch; cin >> ch;
	cameraMatrix = cv::Mat::eye(3, 3, CV_64F);//相机内参数矩阵
	distCoeffs = cv::Mat::zeros(number_distortion_elements, 1, CV_64F);
	Rt = cv::Mat::zeros(3, 4, CV_64F);

	for (int i = 0; i < master_num; i++){



		filein >> name;
		for (int r = 0; r < 3; r++){
			filein >> cameraMatrix.at<double>(r, 0) >> cameraMatrix.at<double>(r, 1)  >> cameraMatrix.at<double>(r, 2);
		}//读取相机内参数

		for (int r = 0; r < 3; r++){
			filein >> Rt.at<double>(r, 0) >> Rt.at<double>(r, 1)  >> Rt.at<double>(r, 2);
		}//读取R

		for (int r = 0; r < 3; r++){
			filein >> Rt.at<double>(r, 3);
		}//读取T

		for (int r = 0; r < number_distortion_elements; r++){
			filein >> distCoeffs.at<double>(r, 0);
		}//读取相机畸变参数

//		cout << "cam matrix before" << endl << cameraMatrix << endl;
//		cameraMatrix = cv::getOptimalNewCameraMatrix(cameraMatrix, distCoeffs, temp_im.size(), 1, temp_im.size(), 0);
//		cout << "cam matrix after" << endl << cameraMatrix << endl;


		cout << "Name " << name << endl;
		// keep K the same ... alter in the projection stuff ....
		BigcameraMatrix = cameraMatrix.clone();
		//		cameraMatrix.at<double>(0, 0) /= downsample_factor;
		//		cameraMatrix.at<double>(0, 2) /= downsample_factor;
		//		cameraMatrix.at<double>(1, 1) /= downsample_factor;
		//		cameraMatrix.at<double>(1, 2) /= downsample_factor;

		if (i == 0){
			master_cameraMatrix = cameraMatrix;
			master_distCoeffs = distCoeffs;
		}

		Rt.at<double>(0, 3) *= multiplier;//统一矩阵T的单位和体素分辨率单位
		Rt.at<double>(1, 3) *= multiplier;
		Rt.at<double>(2, 3) *= multiplier;

		Matrix R(3, 3);
		ColumnVector t(3); ColumnVector Cc(3);
		for (int r = 0; r < 3; r++){
			for (int c = 0; c < 3; c++){
				R(r+ 1, c + 1) = Rt.at<double>(r, c);
			}//newmat库中的Matrix下标从1开始?

			t(r + 1) = Rt.at<double>(r, 3);
		}


		//cout << "cam matrix " << endl << cameraMatrix << endl;


		P = cameraMatrix*Rt;//相机矩阵

		//cout << "P " << endl << P <<endl;


		Cc = -R.i()*t;//i()是一个函数，R.i()表示取R的转置


		C.at<double>(0, 0) =  Cc(1);
		C.at<double>(1, 0) =  Cc(2);
		C.at<double>(2, 0) =  Cc(3);

		//cout << "Cc " << Cc << endl;
		cout << "C " << C << endl;

		ImageData* id = new ImageData();
		id->A = cameraMatrix.clone();
		id->C = C.clone();
		id->RT = Rt.clone();
		id->P = P.clone();
		id->k = distCoeffs.clone();
		id->name = name;

		id_vector.push_back(id);

		if (image_load){//image_load为true
			//id->im_original = cv::imread(name.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
			temp_im = cv::imread(name.c_str(), CV_LOAD_IMAGE_GRAYSCALE);


			// resize .....
			double d   = 0;
			double dist_sum = 0;

			if (i == 0){
//				cv::initUndistortRectifyMap(BigcameraMatrix, distCoeffs, cv::Mat(),
//						cv::getOptimalNewCameraMatrix(BigcameraMatrix, distCoeffs, temp_im.size(), 1, temp_im.size(), 0),
//						temp_im.size(), CV_16SC2, map1, map2);
				cv::initUndistortRectifyMap(BigcameraMatrix, distCoeffs, cv::Mat(),
										BigcameraMatrix,
										temp_im.size(), CV_16SC2, map1, map2);

			}	else {

				d = 0;
				for (int r = 0; r < number_distortion_elements; r++){
					d += fabs(distCoeffs.at<double>(r, 0) - master_distCoeffs.at<double>(r, 0));
					dist_sum += fabs(distCoeffs.at<double>(r, 0));
				}

				//if (dist_sum > 0.2){
					if (d > 0.2){
					cout << "New map ... " << endl;
//					cv::initUndistortRectifyMap(BigcameraMatrix, distCoeffs, cv::Mat(),
//							cv::getOptimalNewCameraMatrix(BigcameraMatrix, distCoeffs, temp_im.size(), 1, temp_im.size(), 0),
//							temp_im.size(), CV_16SC2, map1, map2);
					cv::initUndistortRectifyMap(BigcameraMatrix, distCoeffs, cv::Mat(),
												BigcameraMatrix,
												temp_im.size(), CV_16SC2, map1, map2);
				}
			}

			//cameraMatrix = cv::getOptimalNewCameraMatrix(BigcameraMatrix, distCoeffs, temp_im.size(), 1, temp_im.size(), 0);

			//cv::remap(temp_im, id->im_original, map1, map2, cv::INTER_LINEAR);
			cv::remap(temp_im, id->imcv, map1, map2, cv::INTER_LINEAR);

//			id->imcv = cv::Mat(id->im_original.rows/2, id->im_original.cols/2, CV_8UC1);
//
//			cv::resize(id->im_original, id->imcv, id->imcv.size());


			if (dist_sum > 0.2 || i == 0){


				temp_im= cv::Mat(id->imcv.rows, id->imcv.cols, CV_8UC1, cv::Scalar(255));

				if (dist_sum > 0.5){
					cv::remap(temp_im, id->undistorted_map_original, map1, map2, cv::INTER_LINEAR);
				}	else {
					id->undistorted_map_original = temp_im.clone();
				}

				//id->undistorted_map = cv::Mat(id->im_original.rows/2, id->im_original.cols/2, CV_8UC1);

				//cv::resize(id->undistorted_map_original, id->undistorted_map, id->imcv.size());
			}	else {
				id->undistorted_map_original = id_vector[i - 1]->undistorted_map_original.clone();
				//id->undistorted_map = id_vector[i - 1]->undistorted_map.clone();

			}


			//			cv::imwrite("temp.png", id->undistorted_map);
			//			cv::imwrite("im.png", id->imcv);
			//
			//			char ch; cin >> ch;
			//
			//			if (i == 0){
			//				cv::initUndistortRectifyMap(cameraMatrix, distCoeffs, cv::Mat(),
			//							cv::getOptimalNewCameraMatrix(cameraMatrix, distCoeffs, image_size, 1, image_size, 0),
			//							image_size, CV_16SC2, map1, map2);
			//			}

			// downsample ... then undistort ....



			//id->undistorted_map_original = cv::Mat(id->im_original.rows, id->im_original.cols, CV_8UC1, cv::Scalar(255));
		}
	}

	//cout << "Line 730" << endl;
	//char ch; cin >> ch;
	filein.close();

}


//cv::Mat view, rview, map1, map2;
//	cv::initUndistortRectifyMap(cameraMatrix, distCoeffs, cv::Mat(),
//			cv::getOptimalNewCameraMatrix(cameraMatrix, distCoeffs, image_size, 1, image_size, 0),
//			image_size, CV_16SC2, map1, map2);
//
//	for (int i = 0; i < int(external_images.size()) && i < 2; i++){
//		cout << "Writing external " << i << endl;
//		cv::remap(external_images[i], rview, map1, map2, cv::INTER_LINEAR);
//		filename  = write_directory + "/ext" + ToString<int>(i) + ".png";
//		cv::imwrite(filename.c_str(), rview);
//	}


//Point_2 ImageData::ImageOfCameraCenter(){
//	return Point_2(K.hm(0, 2), K.hm(1, 2));
//}

//void ProcessBinaryImageLowerFactor(IplImage* im, IplImage* result_image, int lower_factor){
//
//	int middle = 256/2;
//	// result image = im.
//	IplImage* im_copy = cvCreateImage( cvSize(im->width, im->height),
//			im->depth, im->nChannels );
//
//	cvCopy(im, im_copy);
//
//	vector<int> percentage(lower_factor, 0);
//
//	// vanishes to zero on lower_factor iter -- so do up to lower factor.
//	for (int i = 0; i < lower_factor; i++){
//		//percentage[i] = 192 - (sin((3.14/2.0) * double(i + 1)/double(lower_factor)) * 192.0);
//		// linear function
//		percentage[i] = double(middle) - double(middle)*(double(i + 1)/double(lower_factor));
//		//percentage[i] = 192 - (sin((3.14/2.0) * double(i + 1)/double(lower_factor)) * 192.0);
//	}
//
//	for (int i = 0; i < lower_factor; i++){
//
//		cvDilate(im_copy, im_copy);
//
//		for (int x = 0; x < im->width; x++){
//			for (int y = 0; y < im->height; y++){
//				if ( GetImageColor(im, y, x, 0) == 0 && GetImageColor(im_copy, y, x, 0) != 0 && GetImageColor(result_image, y, x, 0) == 0){
//
//					for (int ch = 0; ch < 3; ch++){
//						SetImageColor(result_image, y, x, ch, percentage[i]);
//					}
//				}
//			}
//		}
//	}
//}
//
//void ProcessBinaryImageUpperFactor(IplImage* im, IplImage* result_image, int lower_factor){
//
//	int middle = 256/2;
//	// result image = im.
//	IplImage* im_copy = cvCreateImage( cvSize(im->width, im->height),
//			im->depth, im->nChannels );
//
//	cvCopy(im, im_copy);
//
//	vector<int> percentage(lower_factor, 0);
//
//	// vanishes to zero on lower_factor iter -- so do up to lower factor.
//	for (int i = 0; i < lower_factor; i++){
//		//percentage[i] = 192 - (sin((3.14/2.0) * double(i + 1)/double(lower_factor)) * 192.0);
//		// linear function
//		percentage[i] = double(middle) + double(middle)*(double(i + 1)/double(lower_factor));
//		//percentage[i] = 192 - (sin((3.14/2.0) * double(i + 1)/double(lower_factor)) * 192.0);
//	}
//
//	for (int i = 0; i < lower_factor - 1; i++){
//
//		cvErode(im_copy, im_copy);
//
//		for (int x = 0; x < im->width; x++){
//			for (int y = 0; y < im->height; y++){
//				if ( GetImageColor(im, y, x, 0) == 255 && GetImageColor(im_copy, y, x, 0) == 0 && GetImageColor(result_image, y, x, 0) == 255){
//
//					for (int ch = 0; ch < 3; ch++){
//						SetImageColor(result_image, y, x, ch, percentage[i]);
//					}
//				}
//			}
//		}
//	}
//}


//ImageData::ImageData(string filename, int i, double distance, float downsample_factor, bool random_noise_example, double multiplier,
//		int dilation, int upper_factor, int lower_factor){
//	artificial = false;
//	farplane = 4000;
//	//farplane = distance;
//	//farplane = 500;
//	ifstream filein;
//	filein.open(filename.c_str());
//
//	//continuous_im = 0;
//	int master_num;
//	filein >> master_num;
//	//cout << "Master " << master_num << endl;
//
//	if (i <= master_num){
//		index = i;
//		for (int counter = 0; counter < i; counter++){
//			for (int j = 0; j < 18+3+4+1; j++){
//				filein >> name;
//			}
//		}
//
//		double k11, k12, k13, k21, k22, k23, k31, k32, k33;
//		double r11, r12, r13, r21, r22, r23, r31, r32, r33;
//		double t1, t2, t3;
//
//		double rd1, rd2, rd3, rd4;
//
//
//
//
//
//		filein >> name;
//		filein >> k11 >> k12 >> k13 >> k21 >> k22 >> k23 >> k31 >> k32 >> k33;
//		filein >> r11 >> r12 >> r13 >> r21 >> r22 >> r23 >> r31 >> r32 >> r33;
//		filein >> t1 >> t2 >> t3 >> rd1 >> rd2 >> rd3 >> rd4;
//
//		k11 = k11/downsample_factor;
//		k22 = k22/downsample_factor;
//
//		k13 = k13/downsample_factor;
//		k23 = k23/downsample_factor;
//
//		cout << "name " << name << endl;
//
//		Vector_3 t(t1, t2, t3);
//		// only for the temple sequence.... need to write with a new t.
//
//		t = multiplier*t;
//
//		cout << t << endl;
//
//		Aff_Matrix R(r11, r12, r13, r21, r22, r23, r31, r32, r33, 1);
//		Aff_Matrix negR(-r11, -r12, -r13, -r21, -r22, -r23, -r31, -r32, -r33, 1);
//
//		//cout << "Line 342"  << endl;
//
//		RT = Aff_Matrix(r11, r12, r13, t.x(), r21, r22, r23, t.y(), r31, r32, r33, t.z(), 1);
//		RTe = Aff_trans_exact(r11, r12, r13, t.x(), r21, r22, r23, t.y(), r31, r32, r33, t.z(), 1);
//		cout << "Line 346" << endl;
//
//		K = Aff_Matrix(k11, k12, k13, k21, k22, k23, k31, k32, k33, 1);
//		Ke = Aff_trans_exact(k11, k12, k13, k21, k22, k23, k31, k32, k33, 1);
//
//		//cout << "Line 351" << endl;
//		RTie = RTe.inverse();
//		Kie = Ke.inverse();
//		Pe = Ke*RTe;
//
//		P = K*RT;
//
//		//cout << "Line 353" << endl;
//
//		Aff_Matrix rinv;
//		rinv = negR.inverse();
//		//cout << "r inverse .. " << endl;
//		//	MatrixPrint(P);
//
//		//	Vector_3 zero_vec;
//
//		Aff_Matrix curr_RT = RT.inverse();
//		//cout << "rt inverse " << endl;
//		//MatrixPrint(curr_RT);
//		//Vector_3 zero_vec(curr_RT.hm ( 0, 3), curr_RT.hm ( 1, 3), curr_RT.hm ( 2, 3));
//		//zero_ve
//
//		//cout << "P " << endl;
//		//MatrixPrint(P);
//
//		// new stuff.
//		Aff_Matrix A(P.hm(0, 0), P.hm(0, 1), P.hm(0, 2), P.hm(1, 0), P.hm(1, 1), P.hm(1, 2), P.hm(2, 0), P.hm(2, 1), P.hm(2, 2), 1);
//		Aff_Matrix At(A.hm(0, 0), A.hm(1, 0), A.hm(2, 0), A.hm(0, 1), A.hm(1, 1), A.hm(2, 1), A.hm(0, 2), A.hm(1, 2), A.hm(2, 2), 1);
//
//		Aff_Matrix AtA = At*A;
//		Aff_Matrix AtAiAt = AtA.inverse()*At;
//
//		Vector_3 b(-P.hm(0, 3), -P.hm(1, 3), -P.hm(2, 3));
//
//		Vector_3 csolve = AtAiAt.transform(b);
//		//cout << "Line 381" << endl;
//
//		//C = Point_3(curr_RT.hm ( 0, 3), curr_RT.hm ( 1, 3), curr_RT.hm ( 2, 3));
//
//		P1 = Plane_3(P.hm(0, 0), P.hm(0, 1), P.hm(0, 2), P.hm(0, 3));
//		P2 = Plane_3(P.hm(1, 0), P.hm(1, 1), P.hm(1, 2), P.hm(1, 3));
//		PrincipalPlane = Plane_3(P.hm(2, 0), P.hm(2, 1), P.hm(2, 2), P.hm(2, 3));
//
//		P1e = Plane_3_exact(Pe.hm(0, 0), Pe.hm(0, 1), Pe.hm(0, 2), Pe.hm(0, 3));
//		P2e = Plane_3_exact(Pe.hm(1, 0), Pe.hm(1, 1), Pe.hm(1, 2), Pe.hm(1, 3));
//		PrincipalPlaneExact = Plane_3_exact(Pe.hm(2, 0), Pe.hm(2, 1), Pe.hm(2, 2), Pe.hm(2, 3));
//
//		//		CGAL::Object result = CGAL::intersection(P1, P2, PrincipalPlane);
//		//		Point_3 ipoint;
//
//		C = ComputePointFromThreePlanes(P1, P2, PrincipalPlane);
//		ComputePointFromThreePlanes(P1e, P2e, PrincipalPlaneExact, Ce);
//		cout << "C " << C << endl;
//
//		double f = 1.0/K.hm(1, 1);
//		f = 300;
//
//
//		//NearPlane = Plane_3(PrincipalPlane.a(), PrincipalPlane.b(), PrincipalPlane.c(), PrincipalPlane.d() - f);
//
//		//FarPlane = Plane_3(PrincipalPlane.a(), PrincipalPlane.b(), PrincipalPlane.c(), PrincipalPlane.d() - farplane);
//
//		if (name != "temp"){
//			im = cvLoadImage(name.c_str());
//			// ... also, undistort.
//
//			//if (rd1 != 0 || rd2 != 0)
//			{ // this is a raw image ....
//
//				//				cvErode(im, im);
//				//				cvDilate(im, im);
//
//				cout << "Image size is " << im->width << ", " << im->height << endl;
//
//				//				if (dilation > 0){
//				//					cvErode(im, im);
//				//					cvDilate(im, im);
//				//				}
//				//				continuous_im = cvCreateImage( cvSize(im->width, im->height),
//				//						im->depth, im->nChannels );
//				//				cvCopy(im, continuous_im);
//
////				IplImage* temp_image = cvCreateImage( cvSize(im->width, im->height),
////						im->depth, im->nChannels );
//
//				//cvCopy(im, temp_image);
//				//				if (upper_factor > 0 || lower_factor > 0){
//				//					if (lower_factor > 0){
//				//						ProcessBinaryImageLowerFactor(temp_image, continuous_im, lower_factor);
//				//					}
//				//
//				//					if (upper_factor > 0){
//				//						ProcessBinaryImageUpperFactor(temp_image, continuous_im, upper_factor);
//				//					}
//				//				}
//
//				//cvReleaseImage( & temp_image );
//
//
//				//				for (int x = 0; x < dilation; x++){
//				//					cvDilate(im, im);
//				//				}
//
//
//				IplImage* imsmall = cvCreateImage( cvSize(im->width/downsample_factor, im->height/downsample_factor),
//						im->depth, im->nChannels );
//
//				cvResize(im, imsmall);
//				cvReleaseImage( & im );
//				im = imsmall;
//
//				cout << "Reassign" << endl;
//
//				for (int r = 0; r < im->height; r++){
//					for (int c = 0; c < im->width; c++){
//
//						if (GetImageColor(im, r, c, 0) > 0){
//							SetImageColor(im, r, c, 0, 255);
//							SetImageColor(im, r, c, 1, 255);
//							SetImageColor(im, r, c, 2, 255);
//						}
//					}
//				}
//
//
//				IplImage* improc = cvCreateImage( cvSize(im->width, im->height),
//						im->depth, im->nChannels );
//				cout << "Im " << i << " loaded." << endl;
//
//
//				// color all black pixels another color ....
//				for (int x= 0; x < im->width; x++){
//					for (int y = 0; y < im->width; y++){
//						if (GetImageColor(im, y, x, 0) == 0){
//							SetImageColor(im, y, x, 0, 150);
//						}
//					}
//				}
//
//				//				IplImage* imgrey = cvCreateImage( cvSize(im->width, im->height),
//				//														im->depth, im->nChannels );
//				//
//				//				IplImage* imgrey_process = cvCreateImage( cvSize(im->width, im->height),
//				//																		im->depth, im->nChannels );
//
//
//				//cout << distort << endl;
//
//				// undistort
//				cv::Mat M = (cv::Mat_<double>(3,3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
//
//				cv::Mat intrinsic_matrix =  (cv::Mat_<double>(3, 3) << k11, k12,k13,k21, k22, k23, k31, k32, k33);
//				cv::Mat distortion_coeffs = (cv::Mat_<double>(4, 1) << rd1, rd2, rd3, rd4);
//
//
//				cv::Mat image_mat(im);
//				cv::Mat image_proc_mat(improc);
//				//				image_mat = cvLoadImageM(name.c_str(), true);
//				//				//for (int r = 0; )
//
//				//cvRectangle(CvArr* img, CvPoint pt1, CvPoint pt2, CvScalar color, int thickness=1, int lineType=8, int shift=0)
//
//
//
//
//				cv::undistort( image_mat, image_proc_mat, intrinsic_matrix, distortion_coeffs);
//
//
//
//
//
//
//				string temp = "/home/atabb/DemoData/Undistort/" + itoa(i) + "undistorted.ppm";
//				cvSaveImage(temp.c_str(), improc);
//
//				temp = "/home/atabb/DemoData/Undistort/" + itoa(i) + "distorted.ppm";
//				cvSaveImage(temp.c_str(), im);
//
//				cout << "Before release" << endl;
//				//cvReleaseImage( & im);
//				cout << "After release" << endl;
//
//
//				cvCopyImage(improc, im);
//				//cvReleaseImage( &improc );
//				//im = improc;
//
//				for (int x= 0; x < im->width; x++){
//					for (int y = 0; y < im->width; y++){
//						if (GetImageColor(im, y, x, 0) == 0){
//							SetImageColor(im, y, x, 0, 150);
//							SetImageColor(im, y, x, 1, 150);
//							SetImageColor(im, y, x, 2, 150);
//						}	else {
//							if (GetImageColor(im, y, x, 0) == 150){
//								SetImageColor(im, y, x, 0, 0);
//							}	else {
//								SetImageColor(im, y, x, 0, 255);
//							}
//						}
//					}
//				}
//
//				cout << "After transfer" << endl;
//				////////////////////////////////////////////////////////////////////////////////////////
//				//				imsmall = cvCreateImage( cvSize(continuous_im->width/downsample_factor, continuous_im->height/downsample_factor),
//				//						im->depth, im->nChannels );
//				//				cvResize(continuous_im, imsmall);
//				//				cvReleaseImage( & continuous_im );
//				//				continuous_im = imsmall;
//				//
//				//				IplImage* improcc = cvCreateImage( cvSize(imsmall->width, imsmall->height),
//				//						imsmall->depth, imsmall->nChannels );
//				//				cout << "Im " << i << " loaded." << endl;
//				//
//				//				//cout << distort << endl;
//				//
//				//				// undistort
//				//				cv::Mat image_matc(continuous_im);
//				//				cv::Mat image_proc_matc(improcc);
//				//				//				image_mat = cvLoadImageM(name.c_str(), true);
//				//				//				//for (int r = 0; )
//				//
//				//
//				//				cv::undistort( image_matc, image_proc_matc, intrinsic_matrix, distortion_coeffs);
//				//
//				//				temp = "/home/atabb/DemoData/Undistort/" + itoa(i) + "undistorted_continuous.ppm";
//				//				cvSaveImage(temp.c_str(), improcc);
//				//
//				//				cvReleaseImage( & continuous_im);
//				//
//				//
//				//				continuous_im = improcc;
//				//				cvSaveImage(temp.c_str(),continuous_im);
//
//
//
//			}
//
//
//
//
//			// make a 0-1 image
//			//			for (int r = 0; r < im->height; r++){
//			//				for (int c = 0; c < im->width; c++){
//			//
//			//					if (GetImageColor(im, r, c, 0) > 200){
//			//						SetImageColor(im, r, c, 0, 255);
//			//						SetImageColor(im, r, c, 1, 255);
//			//						SetImageColor(im, r, c, 2, 255);
//			//					}	else {
//			//						if (GetImageColor(im, r, c, 0) < 55){
//			//							SetImageColor(im, r, c, 0, 0);
//			//							SetImageColor(im, r, c, 1, 0);
//			//							SetImageColor(im, r, c, 2, 0);
//			//						}	else {
//			//
//			//							if (random_noise_example){
//			//								if (GetImageColor(im, r, c, 0) > 150){
//			//									SetImageColor(im, r, c, 0, 0);
//			//									SetImageColor(im, r, c, 1, 0);
//			//									SetImageColor(im, r, c, 2, 0);
//			//								}	else {
//			//									SetImageColor(im, r, c, 0, 255);
//			//									SetImageColor(im, r, c, 1, 255);
//			//									SetImageColor(im, r, c, 2, 255);
//			//
//			//								}
//			//							}	else {
//			//								if (GetImageColor(im, r, c, 0) > 150){
//			//									SetImageColor(im, r, c, 0, 255);
//			//									SetImageColor(im, r, c, 1, 255);
//			//									SetImageColor(im, r, c, 2, 255);
//			//								}	else {
//			//									SetImageColor(im, r, c, 0, 0);
//			//									SetImageColor(im, r, c, 1, 0);
//			//									SetImageColor(im, r, c, 2, 0);
//			//
//			//								}
//			//
//			//							}
//			//						}
//			//
//			//					}
//			//
//			//
//			//				}
//			//			}
//
//
//
//
//
//
//
//
//			//			if (upper_factor == 0 && lower_factor == 0){
//			//				cvCopy(im, continuous_im);
//			//			}
//
//
//			//			IplImage* destination = cvCreateImage( cvSize(int(im->width/downsample_factor),
//			//									int(im->height/downsample_factor)), im->depth, im->nChannels );
//			//
//			//							cvResize(im, destination);
//			//
//			//							cvReleaseImage( & im );
//			//							im = destination;
//
//
//
//
//			string temp = name + "downsample.png";
//			cvSaveImage(temp.c_str(), im);
//
//			//			temp = name + "continuous.png";
//			//			cvSaveImage(temp.c_str(), continuous_im);
//		}	else {
//			im = 0;
//		}
//
//		FillInContainingPlanes();
//
//
//
//
//	} else {
//		cout << "Error!  Called ImageData(string filename, int i) with an index greater than the number of avaliable images." << endl;
//		char ch;  cin >> ch;
//	}
//
//	cout << "Line 730" << endl;
//	filein.close();
//}


//ImageData::ImageData(string filename, int i, double multiplier, float downsample_factor){
//
//	//cout << "COming through this image data .... " << endl;
//	//char ch;
//	//cin >> ch;
//
//	artificial = false;
//	farplane = 4000;
//	//farplane = distance;
//	//farplane = 500;
//	ifstream filein;
//	filein.open(filename.c_str());
//
//	//continuous_im = 0;
//	int master_num;
//	filein >> master_num;
//	//cout << "Master " << master_num << endl;
//
//	if (i <= master_num){
//		index = i;
//		for (int counter = 0; counter < i; counter++){
//			for (int j = 0; j < 18+3+4+1; j++){
//				filein >> name;
//			}
//		}
//
//		double k11, k12, k13, k21, k22, k23, k31, k32, k33;
//		double r11, r12, r13, r21, r22, r23, r31, r32, r33;
//		double t1, t2, t3;
//
//		//double rd1, rd2, rd3, rd4;
//
//
//
//
//
//		filein >> name;
//		filein >> k11 >> k12 >> k13 >> k21 >> k22 >> k23 >> k31 >> k32 >> k33;
//		filein >> r11 >> r12 >> r13 >> r21 >> r22 >> r23 >> r31 >> r32 >> r33;
//		filein >> t1 >> t2 >> t3 >> rd1 >> rd2 >> rd3 >> rd4;
//
//		cout << setw(10) << "rd1 from Image " << rd1 << endl;
//
//		// keep K the same ... alter in the projection stuff ....
////		k11 = k11/downsample_factor;
////		k22 = k22/downsample_factor;
////
////		k13 = k13/downsample_factor;
////		k23 = k23/downsample_factor;
//
//		cout << "name " << name << endl;
//
//		Vector_3 t(t1, t2, t3);
//		// only for the temple sequence.... need to write with a new t.
//
//		t = multiplier*t;
//
//		cout << t << endl;
//
//		Aff_Matrix R(r11, r12, r13, r21, r22, r23, r31, r32, r33, 1);
//		Aff_Matrix negR(-r11, -r12, -r13, -r21, -r22, -r23, -r31, -r32, -r33, 1);
//
//		//cout << "Line 342"  << endl;
//
//		RT = Aff_Matrix(r11, r12, r13, t.x(), r21, r22, r23, t.y(), r31, r32, r33, t.z(), 1);
//		RTe = Aff_trans_exact(r11, r12, r13, t.x(), r21, r22, r23, t.y(), r31, r32, r33, t.z(), 1);
//		cout << "Line 346" << endl;
//
//		K = Aff_Matrix(k11, k12, k13, k21, k22, k23, k31, k32, k33, 1);
//		Ke = Aff_trans_exact(k11, k12, k13, k21, k22, k23, k31, k32, k33, 1);
//
//		//cout << "Line 351" << endl;
//		RTie = RTe.inverse();
//		Kie = Ke.inverse();
//		Pe = Ke*RTe;
//
//		P = K*RT;
//
//		//cout << "Line 353" << endl;
//
//		Aff_Matrix rinv;
//		rinv = negR.inverse();
//		//cout << "r inverse .. " << endl;
//		//	MatrixPrint(P);
//
//		//	Vector_3 zero_vec;
//
//		Aff_Matrix curr_RT = RT.inverse();
//		//cout << "rt inverse " << endl;
//		//MatrixPrint(curr_RT);
//		//Vector_3 zero_vec(curr_RT.hm ( 0, 3), curr_RT.hm ( 1, 3), curr_RT.hm ( 2, 3));
//		//zero_ve
//
//		//cout << "P " << endl;
//		//MatrixPrint(P);
//
//		// new stuff.
//		Aff_Matrix A(P.hm(0, 0), P.hm(0, 1), P.hm(0, 2), P.hm(1, 0), P.hm(1, 1), P.hm(1, 2), P.hm(2, 0), P.hm(2, 1), P.hm(2, 2), 1);
//		Aff_Matrix At(A.hm(0, 0), A.hm(1, 0), A.hm(2, 0), A.hm(0, 1), A.hm(1, 1), A.hm(2, 1), A.hm(0, 2), A.hm(1, 2), A.hm(2, 2), 1);
//
//		Aff_Matrix AtA = At*A;
//		Aff_Matrix AtAiAt = AtA.inverse()*At;
//
//		Vector_3 b(-P.hm(0, 3), -P.hm(1, 3), -P.hm(2, 3));
//
//		Vector_3 csolve = AtAiAt.transform(b);
//		//cout << "Line 381" << endl;
//
//		//C = Point_3(curr_RT.hm ( 0, 3), curr_RT.hm ( 1, 3), curr_RT.hm ( 2, 3));
//
//		P1 = Plane_3(P.hm(0, 0), P.hm(0, 1), P.hm(0, 2), P.hm(0, 3));
//		P2 = Plane_3(P.hm(1, 0), P.hm(1, 1), P.hm(1, 2), P.hm(1, 3));
//		PrincipalPlane = Plane_3(P.hm(2, 0), P.hm(2, 1), P.hm(2, 2), P.hm(2, 3));
//
//		P1e = Plane_3_exact(Pe.hm(0, 0), Pe.hm(0, 1), Pe.hm(0, 2), Pe.hm(0, 3));
//		P2e = Plane_3_exact(Pe.hm(1, 0), Pe.hm(1, 1), Pe.hm(1, 2), Pe.hm(1, 3));
//		PrincipalPlaneExact = Plane_3_exact(Pe.hm(2, 0), Pe.hm(2, 1), Pe.hm(2, 2), Pe.hm(2, 3));
//
//		//		CGAL::Object result = CGAL::intersection(P1, P2, PrincipalPlane);
//		//		Point_3 ipoint;
//
//		//		C = ComputePointFromThreePlanes(P1, P2, PrincipalPlane);
//		//		ComputePointFromThreePlanes(P1e, P2e, PrincipalPlaneExact, Ce);
//		C = Point_3(csolve.x(), csolve.y(), csolve.z());
//		cout << "C " << C << endl;
//		cout << "C " << C << endl;
//
//		double f = 1.0/K.hm(1, 1);
//		f = 300;
//
//
//		//NearPlane = Plane_3(PrincipalPlane.a(), PrincipalPlane.b(), PrincipalPlane.c(), PrincipalPlane.d() - f);
//
//		//FarPlane = Plane_3(PrincipalPlane.a(), PrincipalPlane.b(), PrincipalPlane.c(), PrincipalPlane.d() - farplane);
//
//		if (name != "temp"){
//			cout << "name " << name << endl;
//			IplImage* raw_image = cvLoadImage(name.c_str());
//			// ... also, undistort.
//			cout << "Loaded ..." << endl;
//
//			cout << "Image size is " << raw_image->width << ", " << raw_image->height << endl;
//
//			//			if (downsample_factor != 1){
//			//				cout << "Downsample not currently set in line 909 of Image Data" << endl;
//			//				exit(1);
//			//			}
//
//			//
//			//			IplImage* imsmall = cvCreateImage( cvSize(raw_image->width/downsample_factor, raw_image->height/downsample_factor),
//			//					raw_image->depth, raw_image->nChannels );
//			//
//			//			cvResize(raw_image, imsmall);
//			//
//			//			cout << "Line 909" << endl;
//			//
//			//			for (int r = 0; r < imsmall->height; r++){
//			//				for (int c = 0; c < imsmall->width; c++){
//			//
//			//					if (GetImageColor(imsmall, r, c, 0) > 0){
//			//						SetImageColor(imsmall, r, c, 0, 255);
//			//						SetImageColor(imsmall, r, c, 1, 255);
//			//						SetImageColor(imsmall, r, c, 2, 255);
//			//					}
//			//				}
//			//			}
//
//			//im = cvCreateImage( cvSize(raw_image->width/downsample_factor, raw_image->height/downsample_factor),
//			//		raw_image->depth, raw_image->nChannels );
//
//			cout << "Undistort map ... " << endl;
//			im = cvCreateImage( cvSize(raw_image->width, raw_image->height),
//					raw_image->depth, raw_image->nChannels );
//			cvZero(im);
//
//			cout << "Line 925" << endl;
//
//			IplImage* imtemp = cvCreateImage( cvSize(raw_image->width, raw_image->height),
//					raw_image->depth, raw_image->nChannels );
//
//			IplImage* imtemp_undistorted = cvCreateImage( cvSize(raw_image->width, raw_image->height),
//					raw_image->depth, raw_image->nChannels );
//
//			cvRectangle( imtemp, cvPoint(0, 0), cvPoint(im->width - 1, im->height - 1), CV_RGB(255, 255, 255), -1 );
//
//			//			// color all black pixels another color ....
//			//						for (int x= 0; x < im->width; x++){
//			//							for (int y = 0; y < im->height; y++){
//			//								if (GetImageColor(imsmall, y, x, 0) == 0){
//			//									SetImageColor(imtemp, y, x, 0, 150);
//			//								}
//			//							}
//			//						}
//
//
//
//
//			//cout << distort << endl;
//
//			// undistort
//			cv::Mat M = (cv::Mat_<double>(3,3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
//
//			// DIFFERENT -- use the original k and image, then resize
////			k11 = k11*downsample_factor;
////			k22 = k22*downsample_factor;
////
////			k13 = k13*downsample_factor;
////			k23 = k23*downsample_factor;
//
//			//			k11 = k11/downsample_factor;
//			//			k22 = k22/downsample_factor;
//			//
//			//			k13 = k13/downsample_factor;
//			//			k23 = k23/downsample_factor;
//			// END DIFFERENT
//			cv::Mat intrinsic_matrix =  (cv::Mat_<double>(3, 3) << k11, k12,k13,k21, k22, k23, k31, k32, k33);
//			cv::Mat distortion_coeffs = (cv::Mat_<double>(4, 1) << rd1, rd2, rd3, rd4);
//
//
//			//cv::Mat image_mat(imsmall);
//			cv::Mat image_mat(raw_image);
//			cv::Mat image_proc_mat(im);
//			cv::Mat temp_mat(imtemp);
//			cv::Mat image_proc_mat2(imtemp_undistorted);
//
//			cv::undistort( image_mat, image_proc_mat, intrinsic_matrix, distortion_coeffs);
//			cv::undistort( temp_mat, image_proc_mat2, intrinsic_matrix, distortion_coeffs);
//
//			cout << "After undistort .... " << endl;
//			//imcv = cv::Mat(im);
//
//			// need to resize ....
//			im_original.copySize(image_proc_mat);
//			image_proc_mat.copyTo(im_original);
//
//			undistorted_map_original.copySize(image_proc_mat2);
//			image_proc_mat2.copyTo(undistorted_map_original);
//
//			if (downsample_factor > 1.0){
//
//
//				cout << "Size of im_original " << im_original.rows << endl;
//
//				// try to save on loadig time ...
//				// dilate
////
////				int dilation_type = cv::MORPH_CROSS;
////				int dilation_size = downsample_factor;
////
////				cv::Mat element = cv::getStructuringElement( dilation_type,
////						cv::Size( 2*dilation_size + 1, 2*dilation_size+1 ),
////						cv::Point( dilation_size, dilation_size ) );
////				/// Apply the dilation operation
////				cv::dilate( image_proc_mat, image_proc_mat, element );
//
//
//
//
//				cv::Mat small_im(im->height/downsample_factor, im->width/downsample_factor, CV_8UC3, cv::Scalar(0,0,0));
//
//				cv::resize(image_proc_mat, image_proc_mat, small_im.size());
//
//				cv::resize(image_proc_mat2, image_proc_mat2, small_im.size());
//			}
//
//
//
//			imcv.copySize(image_proc_mat);
//			image_proc_mat.copyTo(imcv);
//
//			cout << "Size of imcv " << imcv.rows << endl;
//
//
//			undistorted_map.copySize(image_proc_mat);
//			image_proc_mat2.copyTo(undistorted_map);
//
//
//			//undistorted_map.clone(image_proc_mat2);
//			//undistorted_map = image_proc_mat2;
//
//
//			cout << "After transfer" << endl;
//
//			string temp;
////			temp = name + "undistorted.png";
////			cvSaveImage(temp.c_str(), im);
////
////			cv::imwrite(temp.c_str(), imcv);
////
////			temp = name + "undistorted_map.png";
////			cv::imwrite(temp.c_str(), undistorted_map);
//
//
//			//cvReleaseImage( & imsmall );
//			cvReleaseImage( & raw_image );
//			cvReleaseImage( & imtemp );
//			cvReleaseImage( & imtemp_undistorted );
//
//			//			temp = name + "continuous.png";
//			//			cvSaveImage(temp.c_str(), continuous_im);
//		}	else {
//			im = 0;
//		}
//
//		//FillInContainingPlanes();
//
//
//
//
//	} else {
//		cout << "Error!  Called ImageData(string filename, int i) with an index greater than the number of avaliable images." << endl;
//		char ch;  cin >> ch;
//	}
//
//	cout << "Line 730" << endl;
//	filein.close();
//}

//ImageData::ImageData(string filename, int i, double multiplier, bool resize_image, float downsample_factor){
//
//	//cout << "COming through this image data .... " << endl;
//	//char ch;
//	//cin >> ch;
//
//	artificial = false;
//	farplane = 4000;
//	//farplane = distance;
//	//farplane = 500;
//	ifstream filein;
//	filein.open(filename.c_str());
//
//	//continuous_im = 0;
//	int master_num;
//	filein >> master_num;
//	//cout << "Master " << master_num << endl;
//
//	if (i <= master_num){
//		index = i;
//		for (int counter = 0; counter < i; counter++){
//			for (int j = 0; j < 18+3+4+1; j++){
//				filein >> name;
//			}
//		}
//
//		double k11, k12, k13, k21, k22, k23, k31, k32, k33;
//		double r11, r12, r13, r21, r22, r23, r31, r32, r33;
//		double t1, t2, t3;
//
//		//double rd1, rd2, rd3, rd4;
//
//
//
//
//
//		filein >> name;
//		filein >> k11 >> k12 >> k13 >> k21 >> k22 >> k23 >> k31 >> k32 >> k33;
//		filein >> r11 >> r12 >> r13 >> r21 >> r22 >> r23 >> r31 >> r32 >> r33;
//		filein >> t1 >> t2 >> t3 >> rd1 >> rd2 >> rd3 >> rd4;
//
//		cout << setw(10) << "rd1 from Image " << rd1 << endl;
//
//		// keep K the same ... alter in the projection stuff ....
//		k11 = k11/downsample_factor;
//		k22 = k22/downsample_factor;
//
//		k13 = k13/downsample_factor;
//		k23 = k23/downsample_factor;
//
//		cout << "name " << name << endl;
//
//		Vector_3 t(t1, t2, t3);
//		// only for the temple sequence.... need to write with a new t.
//
//		t = multiplier*t;
//
//		cout << t << endl;
//
//		Aff_Matrix R(r11, r12, r13, r21, r22, r23, r31, r32, r33, 1);
//		Aff_Matrix negR(-r11, -r12, -r13, -r21, -r22, -r23, -r31, -r32, -r33, 1);
//
//		//cout << "Line 342"  << endl;
//
//		RT = Aff_Matrix(r11, r12, r13, t.x(), r21, r22, r23, t.y(), r31, r32, r33, t.z(), 1);
//		RTe = Aff_trans_exact(r11, r12, r13, t.x(), r21, r22, r23, t.y(), r31, r32, r33, t.z(), 1);
//		cout << "Line 346" << endl;
//
//		K = Aff_Matrix(k11, k12, k13, k21, k22, k23, k31, k32, k33, 1);
//		Ke = Aff_trans_exact(k11, k12, k13, k21, k22, k23, k31, k32, k33, 1);
//
//		//cout << "Line 351" << endl;
//		RTie = RTe.inverse();
//		Kie = Ke.inverse();
//		Pe = Ke*RTe;
//
//		P = K*RT;
//
//		//cout << "Line 353" << endl;
//
//		Aff_Matrix rinv;
//		rinv = negR.inverse();
//		//cout << "r inverse .. " << endl;
//		//	MatrixPrint(P);
//
//		//	Vector_3 zero_vec;
//
//		Aff_Matrix curr_RT = RT.inverse();
//		//cout << "rt inverse " << endl;
//		//MatrixPrint(curr_RT);
//		//Vector_3 zero_vec(curr_RT.hm ( 0, 3), curr_RT.hm ( 1, 3), curr_RT.hm ( 2, 3));
//		//zero_ve
//
//		//cout << "P " << endl;
//		//MatrixPrint(P);
//
//		// new stuff.
//		Aff_Matrix A(P.hm(0, 0), P.hm(0, 1), P.hm(0, 2), P.hm(1, 0), P.hm(1, 1), P.hm(1, 2), P.hm(2, 0), P.hm(2, 1), P.hm(2, 2), 1);
//		Aff_Matrix At(A.hm(0, 0), A.hm(1, 0), A.hm(2, 0), A.hm(0, 1), A.hm(1, 1), A.hm(2, 1), A.hm(0, 2), A.hm(1, 2), A.hm(2, 2), 1);
//
//		Aff_Matrix AtA = At*A;
//		Aff_Matrix AtAiAt = AtA.inverse()*At;
//
//		Vector_3 b(-P.hm(0, 3), -P.hm(1, 3), -P.hm(2, 3));
//
//		Vector_3 csolve = AtAiAt.transform(b);
//		//cout << "Line 381" << endl;
//
//		//C = Point_3(curr_RT.hm ( 0, 3), curr_RT.hm ( 1, 3), curr_RT.hm ( 2, 3));
//
//		P1 = Plane_3(P.hm(0, 0), P.hm(0, 1), P.hm(0, 2), P.hm(0, 3));
//		P2 = Plane_3(P.hm(1, 0), P.hm(1, 1), P.hm(1, 2), P.hm(1, 3));
//		PrincipalPlane = Plane_3(P.hm(2, 0), P.hm(2, 1), P.hm(2, 2), P.hm(2, 3));
//
//		P1e = Plane_3_exact(Pe.hm(0, 0), Pe.hm(0, 1), Pe.hm(0, 2), Pe.hm(0, 3));
//		P2e = Plane_3_exact(Pe.hm(1, 0), Pe.hm(1, 1), Pe.hm(1, 2), Pe.hm(1, 3));
//		PrincipalPlaneExact = Plane_3_exact(Pe.hm(2, 0), Pe.hm(2, 1), Pe.hm(2, 2), Pe.hm(2, 3));
//
//		//		CGAL::Object result = CGAL::intersection(P1, P2, PrincipalPlane);
//		//		Point_3 ipoint;
//
//		//		C = ComputePointFromThreePlanes(P1, P2, PrincipalPlane);
//		//		ComputePointFromThreePlanes(P1e, P2e, PrincipalPlaneExact, Ce);
//		C = Point_3(csolve.x(), csolve.y(), csolve.z());
//		cout << "C " << C << endl;
//		cout << "C " << C << endl;
//
//		double f = 1.0/K.hm(1, 1);
//		f = 300;
//
//
//		//NearPlane = Plane_3(PrincipalPlane.a(), PrincipalPlane.b(), PrincipalPlane.c(), PrincipalPlane.d() - f);
//
//		//FarPlane = Plane_3(PrincipalPlane.a(), PrincipalPlane.b(), PrincipalPlane.c(), PrincipalPlane.d() - farplane);
//
//		if (name != "temp"){
//			cout << "name " << name << endl;
//			IplImage* raw_image = cvLoadImage(name.c_str());
//			// ... also, undistort.
//			cout << "Loaded ..." << endl;
//
//			cout << "Image size is " << raw_image->width << ", " << raw_image->height << endl;
//
//			//			if (downsample_factor != 1){
//			//				cout << "Downsample not currently set in line 909 of Image Data" << endl;
//			//				exit(1);
//			//			}
//
//			//
//			//			IplImage* imsmall = cvCreateImage( cvSize(raw_image->width/downsample_factor, raw_image->height/downsample_factor),
//			//					raw_image->depth, raw_image->nChannels );
//			//
//			//			cvResize(raw_image, imsmall);
//			//
//			//			cout << "Line 909" << endl;
//			//
//			//			for (int r = 0; r < imsmall->height; r++){
//			//				for (int c = 0; c < imsmall->width; c++){
//			//
//			//					if (GetImageColor(imsmall, r, c, 0) > 0){
//			//						SetImageColor(imsmall, r, c, 0, 255);
//			//						SetImageColor(imsmall, r, c, 1, 255);
//			//						SetImageColor(imsmall, r, c, 2, 255);
//			//					}
//			//				}
//			//			}
//
//			//im = cvCreateImage( cvSize(raw_image->width/downsample_factor, raw_image->height/downsample_factor),
//			//		raw_image->depth, raw_image->nChannels );
//			im = cvCreateImage( cvSize(raw_image->width, raw_image->height),
//					raw_image->depth, raw_image->nChannels );
//			cvZero(im);
//
//			cout << "Line 925" << endl;
//
//			IplImage* imtemp = cvCreateImage( cvSize(raw_image->width, raw_image->height),
//					raw_image->depth, raw_image->nChannels );
//
//			IplImage* imtemp_undistorted = cvCreateImage( cvSize(raw_image->width, raw_image->height),
//					raw_image->depth, raw_image->nChannels );
//
//			cvRectangle( imtemp, cvPoint(0, 0), cvPoint(im->width - 1, im->height - 1), CV_RGB(255, 255, 255), -1 );
//
//			//			// color all black pixels another color ....
//			//						for (int x= 0; x < im->width; x++){
//			//							for (int y = 0; y < im->height; y++){
//			//								if (GetImageColor(imsmall, y, x, 0) == 0){
//			//									SetImageColor(imtemp, y, x, 0, 150);
//			//								}
//			//							}
//			//						}
//
//
//
//
//			//cout << distort << endl;
//
//			// undistort
//			cv::Mat M = (cv::Mat_<double>(3,3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
//
//			// DIFFERENT
//			k11 = k11*downsample_factor;
//			k22 = k22*downsample_factor;
//
//			k13 = k13*downsample_factor;
//			k23 = k23*downsample_factor;
//
//			//			k11 = k11/downsample_factor;
//			//			k22 = k22/downsample_factor;
//			//
//			//			k13 = k13/downsample_factor;
//			//			k23 = k23/downsample_factor;
//			// END DIFFERENT
//			cv::Mat intrinsic_matrix =  (cv::Mat_<double>(3, 3) << k11, k12,k13,k21, k22, k23, k31, k32, k33);
//			cv::Mat distortion_coeffs = (cv::Mat_<double>(4, 1) << rd1, rd2, rd3, rd4);
//
//
//			//cv::Mat image_mat(imsmall);
//			cv::Mat image_mat(raw_image);
//			cv::Mat image_proc_mat(im);
//			cv::Mat temp_mat(imtemp);
//			cv::Mat image_proc_mat2(imtemp_undistorted);
//
//			cv::undistort( image_mat, image_proc_mat, intrinsic_matrix, distortion_coeffs);
//			cv::undistort( temp_mat, image_proc_mat2, intrinsic_matrix, distortion_coeffs);
//
//			//imcv = cv::Mat(im);
//
//			// need to resize ....
//
//			if (downsample_factor > 1.0){
//					image_proc_mat.clone(im_original);
//				// dilate
//
//				  int dilation_type = cv::MORPH_CROSS;
//				  int dilation_size = downsample_factor;
//
//				  cv::Mat element = cv::getStructuringElement( dilation_type,
//				                                       cv::Size( 2*dilation_size + 1, 2*dilation_size+1 ),
//				                                       cv::Point( dilation_size, dilation_size ) );
//				  /// Apply the dilation operation
//				  cv::dilate( image_proc_mat, image_proc_mat, element );
//
//
//
//
//				cv::Mat small_im(im->height/downsample_factor, im->width/downsample_factor, CV_8UC3, cv::Scalar(0,0,0));
//
//				cv::resize(image_proc_mat, image_proc_mat, small_im.size());
//
//				cv::resize(image_proc_mat2, image_proc_mat2, small_im.size());
//			}
//
//
//
//			imcv.copySize(image_proc_mat);
//			image_proc_mat.copyTo(imcv);
//
//
//			undistorted_map.copySize(image_proc_mat);
//			image_proc_mat2.copyTo(undistorted_map);
//
//
//			//undistorted_map.clone(image_proc_mat2);
//			//undistorted_map = image_proc_mat2;
//
//
//			cout << "After transfer" << endl;
//
//			string temp;
//			temp = name + "undistorted.png";
//			cvSaveImage(temp.c_str(), im);
//
//			cv::imwrite(temp.c_str(), imcv);
//
//			temp = name + "undistorted_map.png";
//			cv::imwrite(temp.c_str(), undistorted_map);
//
//
//			//cvReleaseImage( & imsmall );
//			cvReleaseImage( & raw_image );
//			cvReleaseImage( & imtemp );
//			cvReleaseImage( & imtemp_undistorted );
//
//			//			temp = name + "continuous.png";
//			//			cvSaveImage(temp.c_str(), continuous_im);
//		}	else {
//			im = 0;
//		}
//
//		//FillInContainingPlanes();
//
//
//
//
//	} else {
//		cout << "Error!  Called ImageData(string filename, int i) with an index greater than the number of avaliable images." << endl;
//		char ch;  cin >> ch;
//	}
//
//	cout << "Line 730" << endl;
//	filein.close();
//}

//ImageData::ImageData(string filename, int i, double multiplier, float downsample_factor, bool no_image){
//
//	//cout << "COming through this image data .... " << endl;
//	//char ch;
//	//cin >> ch;
//
//	artificial = false;
//	farplane = 4000;
//	//farplane = distance;
//	//farplane = 500;
//	ifstream filein;
//	filein.open(filename.c_str());
//
//	//continuous_im = 0;
//	int master_num;
//	filein >> master_num;
//	//cout << "Master " << master_num << endl;
//
//	if (i <= master_num){
//		index = i;
//		for (int counter = 0; counter < i; counter++){
//			for (int j = 0; j < 18+3+4+1; j++){
//				filein >> name;
//			}
//		}
//
//		double k11, k12, k13, k21, k22, k23, k31, k32, k33;
//		double r11, r12, r13, r21, r22, r23, r31, r32, r33;
//		double t1, t2, t3;
//
//		//double rd1, rd2, rd3, rd4;
//
//
//
//
//
//		filein >> name;
//		filein >> k11 >> k12 >> k13 >> k21 >> k22 >> k23 >> k31 >> k32 >> k33;
//		filein >> r11 >> r12 >> r13 >> r21 >> r22 >> r23 >> r31 >> r32 >> r33;
//		filein >> t1 >> t2 >> t3 >> rd1 >> rd2 >> rd3 >> rd4;
//
//		cout << setw(10) << "rd1 from Image " << rd1 << endl;
//
//		// keep K the same ... alter in the projection stuff ....
//		k11 = k11/downsample_factor;
//		k22 = k22/downsample_factor;
//
//		k13 = k13/downsample_factor;
//		k23 = k23/downsample_factor;
//
//		cout << "name " << name << endl;
//
//		Vector_3 t(t1, t2, t3);
//		// only for the temple sequence.... need to write with a new t.
//
//		t = multiplier*t;
//
//		cout << t << endl;
//
//		Aff_Matrix R(r11, r12, r13, r21, r22, r23, r31, r32, r33, 1);
//		Aff_Matrix negR(-r11, -r12, -r13, -r21, -r22, -r23, -r31, -r32, -r33, 1);
//
//		//cout << "Line 342"  << endl;
//
//		RT = Aff_Matrix(r11, r12, r13, t.x(), r21, r22, r23, t.y(), r31, r32, r33, t.z(), 1);
//		RTe = Aff_trans_exact(r11, r12, r13, t.x(), r21, r22, r23, t.y(), r31, r32, r33, t.z(), 1);
//		cout << "Line 346" << endl;
//
//		K = Aff_Matrix(k11, k12, k13, k21, k22, k23, k31, k32, k33, 1);
//		Ke = Aff_trans_exact(k11, k12, k13, k21, k22, k23, k31, k32, k33, 1);
//
//		//cout << "Line 351" << endl;
//		RTie = RTe.inverse();
//		Kie = Ke.inverse();
//		Pe = Ke*RTe;
//
//		P = K*RT;
//
//		//cout << "Line 353" << endl;
//
//		Aff_Matrix rinv;
//		rinv = negR.inverse();
//		//cout << "r inverse .. " << endl;
//		//	MatrixPrint(P);
//
//		//	Vector_3 zero_vec;
//
//		Aff_Matrix curr_RT = RT.inverse();
//		//cout << "rt inverse " << endl;
//		//MatrixPrint(curr_RT);
//		//Vector_3 zero_vec(curr_RT.hm ( 0, 3), curr_RT.hm ( 1, 3), curr_RT.hm ( 2, 3));
//		//zero_ve
//
//		//cout << "P " << endl;
//		//MatrixPrint(P);
//
//		// new stuff.
//		Aff_Matrix A(P.hm(0, 0), P.hm(0, 1), P.hm(0, 2), P.hm(1, 0), P.hm(1, 1), P.hm(1, 2), P.hm(2, 0), P.hm(2, 1), P.hm(2, 2), 1);
//		Aff_Matrix At(A.hm(0, 0), A.hm(1, 0), A.hm(2, 0), A.hm(0, 1), A.hm(1, 1), A.hm(2, 1), A.hm(0, 2), A.hm(1, 2), A.hm(2, 2), 1);
//
//		Aff_Matrix AtA = At*A;
//		Aff_Matrix AtAiAt = AtA.inverse()*At;
//
//		Vector_3 b(-P.hm(0, 3), -P.hm(1, 3), -P.hm(2, 3));
//
//		Vector_3 csolve = AtAiAt.transform(b);
//		//cout << "Line 381" << endl;
//
//		//C = Point_3(curr_RT.hm ( 0, 3), curr_RT.hm ( 1, 3), curr_RT.hm ( 2, 3));
//
//		P1 = Plane_3(P.hm(0, 0), P.hm(0, 1), P.hm(0, 2), P.hm(0, 3));
//		P2 = Plane_3(P.hm(1, 0), P.hm(1, 1), P.hm(1, 2), P.hm(1, 3));
//		PrincipalPlane = Plane_3(P.hm(2, 0), P.hm(2, 1), P.hm(2, 2), P.hm(2, 3));
//
//		P1e = Plane_3_exact(Pe.hm(0, 0), Pe.hm(0, 1), Pe.hm(0, 2), Pe.hm(0, 3));
//		P2e = Plane_3_exact(Pe.hm(1, 0), Pe.hm(1, 1), Pe.hm(1, 2), Pe.hm(1, 3));
//		PrincipalPlaneExact = Plane_3_exact(Pe.hm(2, 0), Pe.hm(2, 1), Pe.hm(2, 2), Pe.hm(2, 3));
//
//		//		CGAL::Object result = CGAL::intersection(P1, P2, PrincipalPlane);
//		//		Point_3 ipoint;
//
//		//		C = ComputePointFromThreePlanes(P1, P2, PrincipalPlane);
//		//		ComputePointFromThreePlanes(P1e, P2e, PrincipalPlaneExact, Ce);
//		C = Point_3(csolve.x(), csolve.y(), csolve.z());
//		cout << "C " << C << endl;
//		cout << "C " << C << endl;
//
//		double f = 1.0/K.hm(1, 1);
//		f = 300;
//
//
//		//NearPlane = Plane_3(PrincipalPlane.a(), PrincipalPlane.b(), PrincipalPlane.c(), PrincipalPlane.d() - f);
//
//		//FarPlane = Plane_3(PrincipalPlane.a(), PrincipalPlane.b(), PrincipalPlane.c(), PrincipalPlane.d() - farplane);
//
//		im = 0;
//
//
//
//
//	} else {
//		cout << "Error!  Called ImageData(string filename, int i) with an index greater than the number of avaliable images." << endl;
//		char ch;  cin >> ch;
//	}
//
//	cout << "Line 730" << endl;
//	filein.close();
//}


//void ImageData::LoadAndUnDistortThisImage(string filename, IplImage* mask, int downsample){
//
//	name = filename;
//
//	cout << "name " << name << endl;
//	IplImage* raw_image = cvLoadImage(name.c_str());
//	// ... also, undistort.
//	cout << "Loaded ..." << endl;
//
//	cout << "Image size is " << raw_image->width << ", " << raw_image->height << endl;
//
//	if (im == 0){
//		im = cvCreateImage( cvSize(raw_image->width, raw_image->height),
//				raw_image->depth, raw_image->nChannels );
//		cvZero(im);
//	}	else {
//		//		if (downsample != 1){
//		//			cvReleaseImage( & im );
//		//			im = cvCreateImage( cvSize(raw_image->width, raw_image->height),
//		//					raw_image->depth, raw_image->nChannels );
//		//			cvZero(im);
//		//		}
//	}
//
//	cout << "Image size is " << im->width << ", " << im->height << endl;
//
//
//	IplImage* imtemp = cvCreateImage( cvSize(raw_image->width, raw_image->height),
//			raw_image->depth, raw_image->nChannels );
//
//
//	cvRectangle( imtemp, cvPoint(0, 0), cvPoint(im->width - 1, im->height - 1), CV_RGB(255, 255, 255), -1 );
//
//
//	// undistort
//	cv::Mat M = (cv::Mat_<double>(3,3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
//
//	// END DIFFERENT
//	cv::Mat intrinsic_matrix =  (cv::Mat_<double>(3, 3) << K.m(0, 0)*downsample, K.m(0, 1), K.m(0, 2)*downsample, K.m(1, 0),
//			K.m(1, 1)*downsample, K.m(1, 2)*downsample, K.m(2, 0), K.m(2, 1), K.m(2, 2));
//	cv::Mat distortion_coeffs = (cv::Mat_<double>(4, 1) << rd1, rd2, rd3, rd4);
//
//
//	//cv::Mat image_mat(imsmall);
//	cv::Mat image_mat(raw_image);
//	cv::Mat image_proc_mat(im);
//	cv::Mat temp_mat(imtemp);
//	cv::Mat image_proc_mat2(mask);
//
//	cv::undistort( image_mat, image_proc_mat, intrinsic_matrix, distortion_coeffs);
//	cv::undistort( temp_mat, image_proc_mat2, intrinsic_matrix, distortion_coeffs);
//
//	string temp;
//
//	// if temp is < 150, then mark as 150. ... not really necessary.
//	//	for (int r = 0; r < raw_image->height; r++){
//	//		for (int c = 0; c < raw_image->width; c++){
//	//
//	//			if (GetImageColor(mask, r, c, 0) < 100){
//	//				SetImageColor(im, r, c, 0, 150);
//	//				SetImageColor(im, r, c, 1, 150);
//	//				SetImageColor(im, r, c, 2, 150);
//	//			}	else {
//	//				//				if (GetImageColor(im, r, c, 0) > 0){
//	//				//					SetImageColor(im, r, c, 0, 255);
//	//				//					SetImageColor(im, r, c, 1, 255);
//	//				//					SetImageColor(im, r, c, 2, 255);
//	//				//				}
//	//			}
//	//		}
//	//	}
//
//
//	cout << "After transfer" << endl;
//
//
//	temp = name + "undistorted.tif";
//	//cvSaveImage(temp.c_str(), im);
//
//	temp = name + "mask.png";
//	//cvSaveImage(temp.c_str(), mask);
//
//	//cvReleaseImage( & imsmall );
//	cvReleaseImage( & raw_image );
//	cvReleaseImage( & imtemp );
//	//cvReleaseImage( & imtemp_undistorted );
//
//
//
//}
//
//void ImageData::FillInContainingPlanes(){
//	containing_planes.push_back(PrincipalPlane);
//
//	// now, make all of the planes ....
//	int offset = 1;
//	int end_rows = im->height - 1 + offset;
//	int end_cols = im->width - 1 + offset;
//	int start_rows = 0 - offset;
//	int start_cols = 0 - offset;
//
//	Point_2 luc(start_cols, start_rows);
//	Point_2 llc(start_cols, end_rows);
//	Point_2 ruc(end_cols, start_rows);
//	Point_2 rlc(end_cols, end_rows);
//
//	containing_planes.push_back(CalculatePlaneFromImageSegment(Segment_2(llc, luc)));
//
//	containing_planes.push_back(CalculatePlaneFromImageSegment(Segment_2(luc, ruc)));
//
//	containing_planes.push_back(CalculatePlaneFromImageSegment(Segment_2(ruc, rlc)));
//
//	containing_planes.push_back(CalculatePlaneFromImageSegment(Segment_2(rlc, llc)));
//
//}

//void ImageData::Print(){
//
//	cout << "Printing ImageData" << endl;
//	cout << "Name : " << name << endl;
//	cout << "Index : " << index << endl;
//	cout << "K :" << endl;
//	MatrixPrint(K);
//
//	cout << "RT :" << endl;
//	MatrixPrint(RT);
//
//	cout << "P :" << endl;
//	MatrixPrint(P);
//
//	cout << "Principal Plane " << endl;
//	PlanePrint(PrincipalPlane);
//
//	cout<< "P1 " << endl;
//	PlanePrint(P1);
//
//	cout << "P2 " << endl;
//	PlanePrint(P2);
//
//	//	cout << "Ki :" << endl;
//	//	MatrixPrint(Ki);
//	//
//	//	cout << "RTi :" << endl;
//	//	MatrixPrint(RTi);
//	//
//	//	cout << "Pi :" << endl;
//	//	MatrixPrint(Pi);
//
//	cout << "C :" << endl;
//	cout << C << endl;
//
//}

