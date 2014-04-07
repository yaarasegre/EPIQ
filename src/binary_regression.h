/*
 * regression.h
 *
 *  Created on: Aug 13, 2013
 *      Author: yaara
 */

#ifndef REGRESSION_H_
#define REGRESSION_H_

#include <cmath>
#include <iostream>
#include <stddef.h>
#include "typedefs.h"
using namespace std;

class Mat {

	long rows, cols;

	public:

		double** data;

    	inline Mat(){ data = NULL; rows = 0; cols = 0;}
    	Mat(long,long);
    	Mat(long,long, double** memory);
    	~Mat();

    	inline double& operator()(long a, long b) {
    		return data[a][b];
    	}

    	bool mulTransform(Mat const mat2, Mat& result) ;
    	void print() const;
    	inline const int  getrows() const { return rows;}
    	inline const int  getcols() const { return cols;}
    	bool inverse(Mat& result) ;

	private:

//    	bool gluInvertMatrix4(const long double m[16], long double invOut[16]);
    	bool inverse1(Mat& result);
    	bool inverse2(Mat& result);
    	bool inverse3(Mat& result);
    	bool inverse4(Mat& result);

};

/**
 * Linear regression
 */
bool findBeta(Mat& beta, SNP2** X_transpose, Mat  &Y,  long rows, long cols,
		Mat &xtx, Mat &inv, Mat &XtX_X) ;
void modelY(SNP2** X_tr,long rows, long cols, Mat &beta, Mat &result);
double RSS(SNP2** X_tr, Mat  &Y, long rows, long cols,
		Mat &xtx, Mat &inv, Mat &XtX_X, Mat &beta, Mat &predicted_y);
double multipleRegression (SNP2** Xnull_tr,SNP2** Xalt_tr, Mat  &Y,
		Mat &mat_3_3, Mat &mat_3_3b, Mat &mat_2_n, Mat &beta2,
		Mat &mat_4_4, Mat &mat_4_4b, Mat &mat_3_n, Mat &beta3,Mat &predicted_y);
double simpleRegression (SNP2** X_tr, Mat  &Y) ;
double simpleRegression (SNP2** X_tr, Mat  &Y,  //in
		Mat &mat_1_1, Mat &mat_1_1b, Mat &mat_1_n, Mat &beta, Mat &predicted_y) ; // buffers
/*
 * Matrix mult
 */
void XtX( Mat& result, SNP2** X_transpose, long rows, long cols);
bool mul(Mat &res, Mat const &mat1, Mat const &mat2) ;
bool mul( Mat &res, Mat const &mat1, SNP2** const mat2, long rows, long cols) ;

#endif /* REGRESSION_H_ */
