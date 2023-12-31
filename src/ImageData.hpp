/*
 * ImageData.hpp
 *
 *  Created on: Jul 10, 2019
 *      Author: Jinglong
 */

#ifndef IMAGEDATA_HPP_
#define IMAGEDATA_HPP_
//  -DCGAL_USE_F2C -DCGAL_USE_F2C
//#include <cv.h>
//#include <cvaux.h>
//#include <highgui.h>

#include <iostream>
#include <fstream>
#include "newmatap.h"     // newmat advanced functions
#include "newmatio.h"     // newmat headers including output functions
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/highgui/highgui.hpp"

using std::string;

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
using std::vector;

//void Compute3DPlaneFromImagePoints(Point_3& ip1, Point_3& ip2, Point_3& returnpoint, Aff_Matrix& P){
//
//
//
//}
using std::ifstream;

class ImageData{
public:
	string name;
	int index;
	cv::Mat imcv;
	//cv::Mat im_original;
	//cv::Mat undistorted_map;
	cv::Mat undistorted_map_original;
	//cv::Mat mat_im;
	//IplImage* continuous_im;

	cv::Mat A;
	cv::Mat P;
	cv::Mat RT;
	cv::Mat C;
	cv::Mat k;

	~ImageData();

	ImageData();
	//ImageData(string filename, int i, double multiplier, float downsample_factor);

	//ImageData(string filename, int i, double multiplier, float downsample_factor, bool no_image);

	//ImageData(string filename, int i, double multiplier, bool resize_image, float downsample_factor);

	//void Print();

	//void LoadAndUnDistortThisImage(string filename, IplImage* mask, int downsample);

//	bool CalculateLocationOfVertex(const Point& id_vertex);
//
//	Point CalculateTransformedVertex(const Point& id_vertex);
//
//	int CalculateTransformedPlaneOrientation(const Plane_3& id_plane);
//
//	Plane_3 CalculatePlaneFromImageSegment(const Segment_2& image_segment);
//
//	Plane_3_exact CalculatePlaneFromImageSegment(const Segment_2_exact& image_segment);

//	inline Point_2 Project3DPointTo2D(const Point_3& id_vertex){
//		Point_3 x_proj_unnorm = P.transform(id_vertex);
//
//		Point_2 x_proj(x_proj_unnorm.x()/x_proj_unnorm.z(), x_proj_unnorm.y()/ x_proj_unnorm.z());
//
//		return x_proj;
//	}
//
//	inline Point_2_exact Project3DPointTo2DExact(const Point_3& id_vertex){
//		Point_3_exact x_proj_unnorm = Pe.transform(Point_3_exact(id_vertex.x(), id_vertex.y(), id_vertex.z()));
//
//		Point_2_exact x_proj(x_proj_unnorm.x()/x_proj_unnorm.z(), x_proj_unnorm.y()/ x_proj_unnorm.z());
//
//		return x_proj;
//	}
//
//	inline Point_2_exact Project3DPointTo2DExact(const Point_3_exact& id_vertex){
//		Point_3_exact x_proj_unnorm = Pe.transform(id_vertex);
//
//		Point_2_exact x_proj(x_proj_unnorm.x()/x_proj_unnorm.z(), x_proj_unnorm.y()/ x_proj_unnorm.z());
//
//		return x_proj;
//	}

//	Line_2 Project3DPlaneTo2DLine(const Plane_3& p);
//
//	//void FillInContainingPlanes();
//
//	Point_2 ImageOfCameraCenter();
//
//	Point_2_exact CalculateTransformedVertexExact(const Point_3_exact& id_vertex);
//
//	Segment_2_exact BackprojectImageSegment(Point_3_exact X1, Point_3_exact X2);
};


template<class T>
string ToString(T arg)
{
	std::ostringstream s;

	s << arg;

	return s.str();

}


template<class T>
T FromString(const std::string& s)
{
	std::istringstream stream (s);
	T t;
	stream >> t;
	return t;


}


void LoadImageDataVector(vector<ImageData*>& id_vector, string filename,double multiplier, float downsample_factor, bool image_load);


void DeleteImageData(ImageData* id);

string ReturnModelFileName(string infile);

//inline Point_2 Project3DPointTo2D(ImageData& id, const Point_3& id_vertex){
//	Point_3 x_proj_unnorm = id.P.transform(id_vertex);
//
//	Point_2 x_proj(x_proj_unnorm.x()/x_proj_unnorm.z(), x_proj_unnorm.y()/ x_proj_unnorm.z());
//
//
//
//	return x_proj;
//
//}

//
//inline Point_2 Project3DPointTo2DThreadsafe(ImageData& id, const Point_3& id_vertex){
//
//	double x, y, z;
//
//	x = id.P.hm(0,0)*id_vertex.x() + id.P.hm(0,1)*id_vertex.y() + id.P.hm(0,2)*id_vertex.z() + id.P.hm(0, 3);
//	y = id.P.hm(1,0)*id_vertex.x() + id.P.hm(1,1)*id_vertex.y() + id.P.hm(1,2)*id_vertex.z() + id.P.hm(1, 3);
//	z = id.P.hm(2,0)*id_vertex.x() + id.P.hm(2,1)*id_vertex.y() + id.P.hm(2,2)*id_vertex.z() + id.P.hm(2, 3);
//
//	Point_2 x_proj(x/z, y/z);
//
//	return x_proj;
//
//}
//
//inline Point_2 Project3DPointTo2DThreadsafe(Aff_Matrix& P, const Point_3& id_vertex){
//
//	double x, y, z;
//
//	x = P.hm(0,0)*id_vertex.x() + P.hm(0,1)*id_vertex.y() + P.hm(0,2)*id_vertex.z() + P.hm(0, 3);
//	y = P.hm(1,0)*id_vertex.x() + P.hm(1,1)*id_vertex.y() + P.hm(1,2)*id_vertex.z() + P.hm(1, 3);
//	z = P.hm(2,0)*id_vertex.x() + P.hm(2,1)*id_vertex.y() + P.hm(2,2)*id_vertex.z() + P.hm(2, 3);
//
//	Point_2 x_proj(x/z, y/z);
//
//	return x_proj;
//
//}

//Point_2 ImageData::Project3DPointTo2D(const Point_3& id_vertex){
//	Point_3 x_proj_unnorm = P.transform(id_vertex);
//
//	Point_2 x_proj(x_proj_unnorm.x()/x_proj_unnorm.z(), x_proj_unnorm.y()/ x_proj_unnorm.z());
//
//	//trans = x_proj;
//
//	return x_proj;
//
//}

//Point_2 Project3DPointTo2D(ImageData& id, const Point_3& id_vertex){
//	Point_3 x_proj_unnorm = id.P.transform(id_vertex);
//
//	Point_2 x_proj(x_proj_unnorm.x()/x_proj_unnorm.z(), x_proj_unnorm.y()/ x_proj_unnorm.z());
//
//	//trans = x_proj;
//
//	return x_proj;
//}

#endif /* IMAGEDATA_HPP_ */
