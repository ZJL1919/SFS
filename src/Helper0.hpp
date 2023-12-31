/*
 * Helper0.hpp
 *
 *  Created on: Nov 19, 2019
 *      Author: Jinglong
 */

#ifndef HELPER0_HPP_
#define HELPER0_HPP_

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



inline void MultiplyMatrixVector(const vector< vector<double> >& M, int rows, int cols, const vector<double>& v, vector<double>& X){

	int c;
	// don't test to save time
	for (int r = 0; r < rows; r++){
		X[r] = 0;

		for (c = 0; c < cols; c++){
			X[r] += M[r][c]*v[c];
		}
	}
}

inline void MultiplySquareMatrixMatrix(const vector< vector<double> >& M0, const vector< vector<double> >& M1,
		int rows, vector< vector<double> >& R){

	int c;
	double r;
	// don't test to save time
	for (int r = 0; r < rows; r++){


		for (c = 0; c < rows; c++){
			// row r in M0 * col c in M1.
			R[r][c] = 0;

			for (int inc = 0; inc < rows; inc++){
				R[r][c] += M0[r][inc]*M1[inc][c];
			}


		}
	}
}


inline void SubtractVectorFromVector(const vector<double>& A, const vector<double>& B, vector<double>& C, int rows){
	for (int r = 0; r < rows; r++){
		C[r] = A[r] - B[r];
	}
}

inline void AddVectorToVector(const vector<double>& A, const vector<double>& B, vector<double>& C, int rows){
	for (int r = 0; r < rows; r++){
		C[r] = A[r] + B[r];
	}
}

inline void MultiplyVectorByScalar(const vector<double>& A, const double scalar, vector<double>& B, int rows){
	for (int r = 0; r < rows; r++){
		B[r] = scalar*A[r];
	}

}

inline double SquaredDistance(const vector<double>& A, const vector<double>& B, int rows){
	double d = 0;
	for (int r = 0; r < rows; r++){
		d += (A[r] - B[r])*(A[r] - B[r]);
	}
	return d;
}



inline bool ProjectPointAndReturnIndex(const vector< vector<double> >& P,
		vector<double>& X, vector<double>& x, int rows, int cols, int_type_t& pixel_index){

	bool in = false;
	int r, c;
	MultiplyMatrixVector(P, 3, 4, X, x);


	//cout << "big X " << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << endl;
	//cout << "image x " << x[0] << " " << x[1] << " " << x[2] <<  endl;
	x[0] /= x[2];  /// c
	x[1] /= x[2];  /// r
	x[2] = 1;

	//cout << "image x " << x[0] << " " << x[1] << endl;
	c = round(x[0]);//   ﷽      ת       ص    ꣬        
	r = round(x[1]);

	if (c >= 0 && c < cols && r >= 0 && r < rows){
		in = true;
		pixel_index = round(x[1])*cols + round(x[0]);//   ַ ʽ    pixel_index
	}	else {
		pixel_index = 0;//   ͶӰ  ͼ   ⣬  pixel_index  Ϊ0
	}

	return in;
}

inline void NormalizePlane(vector<double>& p){//ʲô  ˼    

	double mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

	for (int i = 0; i < 4; i++){
		p[i] /= mag;
	}
}

inline void NormalizeVector(vector<double>& p){

	double mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

	if (fabs(mag) > 0.00000001){
	for (int i = 0; i < 3; i++){
		p[i] /= mag;
	}
	}
}

inline double MagnitudeVector(vector<double>& p){
	double mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
	return mag;
}

inline double DotProduct(const vector<double>& a, const vector<double>& b, int size){

	double r = 0;

	for (int i = 0; i < size; i++){
		r += a[i]*b[i];
	}
	return r;
}

inline void CrossProduct(const vector<double>& a, const vector<double>& b, vector<double>& c){//   

	c[0] = a[1]*b[2] - b[1]*a[2];
	c[1] = -a[0]*b[2] + b[0]*a[2];
	c[2] = a[0]*b[1] - b[0]*a[1];

}

inline void RayPlaneIntersection(const vector<double>& C, const vector<double>& V, const vector<double>& P, vector<double>& X ){//ʲô  ˼    

	double dp0 = DotProduct(C, P, 4);
	double dp1 = DotProduct(V, P, 4);

	if (dp1 != 0){
		double lambda = -dp0/dp1;

		X[0] = C[0] + lambda*V[0];
		X[1] = C[1] + lambda*V[1];
		X[2] = C[2] + lambda*V[2];


	}	else {

		X[0] = C[0];
		X[1] = C[1];
		X[2] = C[2];
		cout << "Slight error -- cannot find appropriate lambda b/c dot product is zero: " << endl;
		cout << "top " << dp0 << "   bottom " << dp1 << endl;
	}
}

inline bool RayPlaneIntersectionBool(const vector<double>& C, const vector<double>& V, const vector<double>& P, vector<double>& X ){

	double dp0 = DotProduct(C, P, 4);
	double dp1 = DotProduct(V, P, 4);

	if (dp1 != 0){
		double lambda = -dp0/dp1;

		X[0] = C[0] + lambda*V[0];
		X[1] = C[1] + lambda*V[1];
		X[2] = C[2] + lambda*V[2];

		if (lambda < 0){
			return false;
		}	else {
			return true;
		}

	}	else {

		X[0] = C[0];
		X[1] = C[1];
		X[2] = C[2];
		cout << "Slight error -- cannot find appropriate lambda b/c dot product is zero: " << endl;
		cout << "top " << dp0 << "   bottom " << dp1 << endl;

		return false;
	}
}

#endif /* HELPER0_HPP_ */
