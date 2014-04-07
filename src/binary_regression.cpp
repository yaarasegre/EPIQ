/*
 * regression.cpp
 *
 *  Created on: Aug 13, 2013
 *      Author: yaara
 */
#include "binary_regression.h"


bool findBeta(Mat& beta, SNP2** X_transpose, Mat  &Y,  long rows, long cols,
		Mat &xtx, Mat &inv, Mat &XtX_X) {

//	cout << "x transpose: " << endl;
//	for(int i = 0; i<rows;i++) {
//		for(int j = 0; j<cols;j++) {
//			int a = X_transpose[j][i];
//			cout << a << " " ;
//		}
//		cout << endl;
//	}
//	cout  << endl;

	XtX(xtx, X_transpose, rows, cols);
//	cout << "x transpose x: " << endl;
//	xtx.print();
	if (!xtx.inverse(inv))
		return false; // singular
	if (!mul(XtX_X, inv, X_transpose, cols, rows))
		return false;
	if (!mul(beta, XtX_X, Y))
		return false;

//	cout<< "x";
//	X.print();
//	cout << "xtx";
//	xtx.print();
//	cout << "inv xtx";
//	inv.print();
//	cout << "inv xtx * xt";
//	XtX_X.print();
//	cout << "beta";
//	beta.print();

	return true;
}

/**
 * Predict Y:
 * Y^ = X*beta^
 * Result is written to Mat &result.
 */
void modelY(SNP2** X_tr,long rows, long cols, Mat &beta, Mat &result) {

	long double sum ;

	for (long i = 0 ; i < rows ; i++ ) {
		sum = 0;
		for (long k = 0 ; k < cols ; k++ )
		{
				sum = sum + X_tr[k][i] * beta.data[k][0];
		}
		result.data[i][0] = sum;
	}
}

/**
 * Linear regression RSS
 */
double RSS(SNP2** X_tr, Mat  &Y, long rows, long cols,
		Mat &xtx, Mat &inv, Mat &XtX_X, Mat &beta, Mat &predicted_y) {

	if (!findBeta(beta, X_tr, Y, rows, cols, xtx, inv, XtX_X))
		return -1;

	modelY(X_tr, rows, cols, beta, predicted_y);

	long double sum = 0;
//	cout << "predicted y" << endl;
	for (long i = 0; i < Y.getrows(); ++i) {
		sum += pow(Y(i,0) - predicted_y(i,0), 2);
//		cout << predicted_y(i,0);
	}
	return sum;
}

/**
 * Multiple regression.
 * return value: -2ln GLR statistic. Assumes first col of X contains ones!
 * Parameters after Y are buffers for the functions.
 */
double multipleRegression (SNP2** Xnull_tr,SNP2** Xalt_tr, Mat  &Y, Mat &mat_3_3, Mat &mat_3_3b, Mat &mat_2_n, Mat &beta2,
		Mat &mat_4_4, Mat &mat_4_4b, Mat &mat_3_n, Mat &beta3,Mat &predicted_y) {
//	double rssalt= log(RSS(Xalt, Y));
//	double rssnull = log(RSS(Xnull, Y));
//	 int yr =  Y.getrows();
	double rss_null = RSS(Xnull_tr, Y, Y.getrows(), 3, mat_3_3, mat_3_3b, mat_2_n, beta2, predicted_y);
	if (rss_null == -1)
		return 0;


	double rss_alt = RSS(Xalt_tr, Y, Y.getrows(), 4, mat_4_4, mat_4_4b, mat_3_n, beta3, predicted_y);
	if (rss_alt == -1)
			return 0;
	double ret =  -(double)Y.getrows() * (log(rss_alt) - log(rss_null));
	return ret;
}

/**
 * Simple regression. Assumes first col of X contains ones!
 * Null model: y = 0 + err.
 * return value: LRT statistic
 * Parameters after Y are buffers for the functions.
 */
double simpleRegression (SNP2** X_tr, Mat  &Y,
		Mat &mat_2_2, Mat &mat_2_2b, Mat &mat_2_n, Mat &beta, Mat &predicted_y) {

	double rss_alt = RSS(X_tr, Y, Y.getrows(), 2, mat_2_2, mat_2_2b, mat_2_n, beta, predicted_y);
	if (rss_alt == -1)
		return 0;

	// calc RSS null
	double rss_null = 0;
	for (long i=0; i<Y.getrows(); i++ ) {
		rss_null += pow(Y(i,0),2); // assuming mean(Y) is 0!!!
	}
	double ret =  -(double)Y.getrows() * (log(rss_alt) - log(rss_null));
	return ret;
}


/**
 * Simple regression wrapper
 */
double simpleRegression (SNP2** X_tr, Mat  &Y) {

	// buffers
	Mat mat_2_2(2,2);
	Mat mat_2_2b(2,2);
	Mat mat_2_n(2,Y.getrows());
	Mat beta(2,2);
	Mat predicted_y(Y.getrows(),1);

	return  simpleRegression (X_tr, Y, mat_2_2, mat_2_2b, mat_2_n, beta, predicted_y) ;

}
/**
 * Constructor
 */
Mat::Mat (long n, long m) {
	// allocate matrix
	rows = n; cols = m;
	data = NULL;
	data = new double*[n];
	if (data)
	{
		data[0] = new double[n * m]; // allocate all at once
	    for (long i = 0; i < n; ++i)
	    	data[i] = data[0] + i * m; // first element in row points at the row
	} else {
		cerr << "Matrix allocation failed " << endl;
	}
}

/**
 * Constructor 2
 */
Mat::Mat (long n, long m, double** memory) {
	// allocate matrix
	rows = n; cols = m;
	data = memory;
}

/*
 * Destructor
 */
Mat::~Mat () {
	if (data != NULL && rows !=0) {
		if (data[0] != NULL)
			delete [] data[0];
		delete [] data;
	}
}


void Mat::print() const {

	cout << endl;
	for (long i = 0; i < rows; i++) {
		for (long j = 0; j < cols; j++) {
			cout << data[i][j] << "\t" ;
//			cerr << i;
		}
		cout << endl;
	}
	cout << endl;

}
/**
 * Matrix multiply with transform of mat2
 */
bool Mat::mulTransform(Mat const mat2,Mat& result) {

	if (this->cols != mat2.cols) {
		return false;
	}
	long double sum = 0;

	for (long i = 0 ; i < this->rows ; i++ )
	    {
	      for (long k = 0 ; k < mat2.rows ; k++ )
	      {
	        for (long j = 0 ; j < this->cols  ; j++ )
	        {
	          sum = sum + this->data[i][j] * mat2.data[k][j];
	        }

	        result.data[i][k] = sum;
	        sum = 0;
	      }
	    }


	return true;
}




/**
 * Invert  matrices up to 4x4
 */
bool Mat::inverse(Mat& result)  {

	bool success = false;
	if (cols == 1 && rows == 1) {
			success = inverse1(result);
	}
	if (cols == 2 && rows == 2) {
		success = inverse2(result);
	}
	else if(cols == 3 && rows == 3)  {
		success = inverse3(result);
	}
	else if(cols == 4 && rows == 4)  {
		success = inverse4(result);
	}
	return success;
}

bool Mat::inverse1(Mat& result) {

	if (data[0][0] == 0)
		return false;

	result.data[0][0] = 1/data[0][0];
	return true;
}


bool Mat::inverse2(Mat& result) {

	long double determinant =    data[0][0]*data[1][1]-data[1][0]*data[0][1];
	if (determinant == 0)
		return false;

	long double invdet = 1/determinant;
	result.data[0][0] =  data[1][1] * invdet;
	result.data[0][1] = -data[0][1] * invdet;
	result.data[1][0] = -data[1][0] * invdet;
	result.data[1][1] =  data[0][0] * invdet;

	return true;
}


bool Mat::inverse3(Mat& result) {

	long double determinant =    +data[0][0]*(data[1][1]*data[2][2]-data[2][1]*data[1][2])
	                        -data[0][1]*(data[1][0]*data[2][2]-data[1][2]*data[2][0])
	                        +data[0][2]*(data[1][0]*data[2][1]-data[1][1]*data[2][0]);
	if (determinant == 0)
		return false;

	long double invdet = 1/determinant;
	result.data[0][0] =  (data[1][1]*data[2][2]-data[2][1]*data[1][2]) * invdet;
	result.data[0][1] = -(data[0][1]*data[2][2]-data[0][2]*data[2][1]) * invdet;
	result.data[0][2] =  (data[0][1]*data[1][2]-data[0][2]*data[1][1]) * invdet;
	result.data[1][0] = -(data[1][0]*data[2][2]-data[1][2]*data[2][0]) * invdet;
	result.data[1][1] =  (data[0][0]*data[2][2]-data[0][2]*data[2][0]) * invdet;
	result.data[1][2] = -(data[0][0]*data[1][2]-data[1][0]*data[0][2]) * invdet;
	result.data[2][0] =  (data[1][0]*data[2][1]-data[2][0]*data[1][1]) * invdet;
	result.data[2][1] = -(data[0][0]*data[2][1]-data[2][0]*data[0][1]) * invdet;
	result.data[2][2] =  (data[0][0]*data[1][1]-data[1][0]*data[0][1]) * invdet;

	return true;
}


/**
 * inverse matrix (Ugly)
 *
 */
bool Mat::inverse4(Mat& result)
// const long double data[16], long double invOut[16])
{
	long double inv[16], det;
	int i,j;

	long double tmp[16];
	for (i = 0; i < 4; i++)
		for ( j = 0; j < 4; j++)
			tmp[i*4+j] = data[i][j];


	inv[0] = tmp[5]  * tmp[10] * tmp[15] -
			tmp[5]  * tmp[11] * tmp[14] -
			tmp[9]  * tmp[6]  * tmp[15] +
			tmp[9]  * tmp[7]  * tmp[14] +
			tmp[13] * tmp[6]  * tmp[11] -
			tmp[13] * tmp[7]  * tmp[10];

	inv[4] = -tmp[4]  * tmp[10] * tmp[15] +
			tmp[4]  * tmp[11] * tmp[14] +
			tmp[8]  * tmp[6]  * tmp[15] -
			tmp[8]  * tmp[7]  * tmp[14] -
			tmp[12] * tmp[6]  * tmp[11] +
			tmp[12] * tmp[7]  * tmp[10];

	inv[8] = tmp[4]  * tmp[9] * tmp[15] -
			tmp[4]  * tmp[11] * tmp[13] -
			tmp[8]  * tmp[5] * tmp[15] +
			tmp[8]  * tmp[7] * tmp[13] +
			tmp[12] * tmp[5] * tmp[11] -
			tmp[12] * tmp[7] * tmp[9];

	inv[12] = -tmp[4]  * tmp[9] * tmp[14] +
			tmp[4]  * tmp[10] * tmp[13] +
			tmp[8]  * tmp[5] * tmp[14] -
			tmp[8]  * tmp[6] * tmp[13] -
			tmp[12] * tmp[5] * tmp[10] +
			tmp[12] * tmp[6] * tmp[9];

	inv[1] = -tmp[1]  * tmp[10] * tmp[15] +
			tmp[1]  * tmp[11] * tmp[14] +
			tmp[9]  * tmp[2] * tmp[15] -
			tmp[9]  * tmp[3] * tmp[14] -
			tmp[13] * tmp[2] * tmp[11] +
			tmp[13] * tmp[3] * tmp[10];

	inv[5] = tmp[0]  * tmp[10] * tmp[15] -
			tmp[0]  * tmp[11] * tmp[14] -
			tmp[8]  * tmp[2] * tmp[15] +
			tmp[8]  * tmp[3] * tmp[14] +
			tmp[12] * tmp[2] * tmp[11] -
			tmp[12] * tmp[3] * tmp[10];

	inv[9] = -tmp[0]  * tmp[9] * tmp[15] +
			tmp[0]  * tmp[11] * tmp[13] +
			tmp[8]  * tmp[1] * tmp[15] -
			tmp[8]  * tmp[3] * tmp[13] -
			tmp[12] * tmp[1] * tmp[11] +
			tmp[12] * tmp[3] * tmp[9];

	inv[13] = tmp[0]  * tmp[9] * tmp[14] -
			tmp[0]  * tmp[10] * tmp[13] -
			tmp[8]  * tmp[1] * tmp[14] +
			tmp[8]  * tmp[2] * tmp[13] +
			tmp[12] * tmp[1] * tmp[10] -
			tmp[12] * tmp[2] * tmp[9];

	inv[2] = tmp[1]  * tmp[6] * tmp[15] -
			tmp[1]  * tmp[7] * tmp[14] -
			tmp[5]  * tmp[2] * tmp[15] +
			tmp[5]  * tmp[3] * tmp[14] +
			tmp[13] * tmp[2] * tmp[7] -
			tmp[13] * tmp[3] * tmp[6];

	inv[6] = -tmp[0]  * tmp[6] * tmp[15] +
			tmp[0]  * tmp[7] * tmp[14] +
			tmp[4]  * tmp[2] * tmp[15] -
			tmp[4]  * tmp[3] * tmp[14] -
			tmp[12] * tmp[2] * tmp[7] +
			tmp[12] * tmp[3] * tmp[6];

	inv[10] = tmp[0]  * tmp[5] * tmp[15] -
			tmp[0]  * tmp[7] * tmp[13] -
			tmp[4]  * tmp[1] * tmp[15] +
			tmp[4]  * tmp[3] * tmp[13] +
			tmp[12] * tmp[1] * tmp[7] -
			tmp[12] * tmp[3] * tmp[5];

	inv[14] = -tmp[0]  * tmp[5] * tmp[14] +
			tmp[0]  * tmp[6] * tmp[13] +
			tmp[4]  * tmp[1] * tmp[14] -
			tmp[4]  * tmp[2] * tmp[13] -
			tmp[12] * tmp[1] * tmp[6] +
			tmp[12] * tmp[2] * tmp[5];

	inv[3] = -tmp[1] * tmp[6] * tmp[11] +
			tmp[1] * tmp[7] * tmp[10] +
			tmp[5] * tmp[2] * tmp[11] -
			tmp[5] * tmp[3] * tmp[10] -
			tmp[9] * tmp[2] * tmp[7] +
			tmp[9] * tmp[3] * tmp[6];

	inv[7] = tmp[0] * tmp[6] * tmp[11] -
			tmp[0] * tmp[7] * tmp[10] -
			tmp[4] * tmp[2] * tmp[11] +
			tmp[4] * tmp[3] * tmp[10] +
			tmp[8] * tmp[2] * tmp[7] -
			tmp[8] * tmp[3] * tmp[6];

	inv[11] = -tmp[0] * tmp[5] * tmp[11] +
			tmp[0] * tmp[7] * tmp[9] +
			tmp[4] * tmp[1] * tmp[11] -
			tmp[4] * tmp[3] * tmp[9] -
			tmp[8] * tmp[1] * tmp[7] +
			tmp[8] * tmp[3] * tmp[5];

	inv[15] = tmp[0] * tmp[5] * tmp[10] -
			tmp[0] * tmp[6] * tmp[9] -
			tmp[4] * tmp[1] * tmp[10] +
			tmp[4] * tmp[2] * tmp[9] +
			tmp[8] * tmp[1] * tmp[6] -
			tmp[8] * tmp[2] * tmp[5];

	det = tmp[0] * inv[0] + tmp[1] * inv[4] + tmp[2] * inv[8] + tmp[3] * inv[12];

	if (det == 0)
		return false;

	det = 1.0 / det;

	for (i = 0; i < 4; i++)
		for ( j = 0; j < 4; j++)
			result.data[i][j] = inv[i*4+j] * det;

	return true;
}

/**
 * Matrix X transpose * X
 * Takes advantage of the fact that X is binary. (does && instead of *)
 */
void XtX( Mat& result, SNP2** X_transpose, long rows, long cols) {

	long double sum = 0;
	double mul;

	for (long i = 0 ; i < cols ; i++ ) {
		for (long k = 0 ; k < cols ; k++ )
		{
			for (long j = 0 ; j < rows  ; j++ )
			{
				mul = (X_transpose[i][j] && X_transpose[k][j] ? 1.0 : 0.0);
				sum = sum + mul;
//				sum = sum + X_transpose[i][j] * X_transpose[k][j];
			}

			result.data[i][k] = sum;
			sum = 0;
		}
	}
}



///**
// * Matrix mult
// */
//bool mul(Mat &res, Mat const &mat1, Mat const &mat2) {
//
//	if (mat1.getcols() != mat2.getrows()) {
//		return false;
//	}
//
//	for (long i = 0 ; i < mat1.getrows() ; i++ )
//		for (long k = 0 ; k < mat2.getcols() ; k++ )
//			res.data[i][k] = 0.0;
////	cout << "i = 1:" << mat1.getrows();
////	cout << " j = 1:" << mat1.getcols(); 1000
////	cout << " k = 1:" << mat2.getcols()<<endl;
//
//	for (long i = 0 ; i < mat1.getrows() ; i++ )
//		for (long k = 0 ; k < mat2.getcols() ; k++ )
//			for (long j = 0 ; j < mat1.getcols()  ; j++ )
//				res.data[i][k] += mat1.data[i][j] * mat2.data[j][k];
//	return true;
//}

/**
 * Matrix mult
 */
bool mul(Mat &res, Mat const &mat1, Mat const &mat2) {

	if (mat1.getcols() != mat2.getrows()) {
		return false;
	}
	double sum = 0;

	for (long i = 0 ; i < mat1.getrows() ; i++ )
	    {
	      for (long k = 0 ; k < mat2.getcols() ; k++ )
	      {
	        for (long j = 0 ; j < mat1.getcols()  ; j++ )
	        {
	          sum = sum + mat1.data[i][j] * mat2.data[j][k];
	        }

	        res.data[i][k] = sum;
	        sum = 0;
	      }
	    }
	return true;

}


///**
// * Matrix multiplication  override
// */
//bool mul( Mat &res, Mat const &mat1, SNP2** const mat2, long rows, long cols) {
//
//	if (mat1.getcols() != rows)
//		return false;
//
//	// init
//	for (long i = 0 ; i < mat1.getrows() ; i++ )
//		for (long k = 0 ; k < cols ; k++ )
//			res.data[i][k] = 0.0;
////	cout << "~i = 1:" << mat1.getrows(); // ~i = 1:4 j = 1:4 k = 1:1000
////	cout << " j = 1:" << mat1.getcols();
////	cout << " k = 1:" <<cols<<endl;
//	// mul
//	double sum = 0;
//	for (long i = 0 ; i < mat1.getrows() ; i++ )
//	{
//		for (long k = 0 ; k < cols ; k++ )
//		{
//			for (long j = 0 ; j < mat1.getcols()  ; j++ )
//			{
//				sum += (mat2[j][k]==0? 0: mat1.data[i][j]); // (same as mat1.data[i][j] * mat2.data[j][k])
//			}
//
//			res.data[i][k] = sum;
//			sum = 0;
//		}
//	}
//	return true;
//}

/**
 * Matrix multiplication  override
 */
bool mul( Mat &res, Mat const &mat1, SNP2** const mat2, long rows, long cols) {

	if (mat1.getcols() != rows)
		return false;

	// init
	for (long i = 0 ; i < mat1.getrows() ; i++ )
		for (long k = 0 ; k < cols ; k++ )
			res.data[i][k] = 0.0;
//	cout << "~i = 1:" << mat1.getrows(); // ~i = 1:4 j = 1:4 k = 1:1000
//	cout << " j = 1:" << mat1.getcols();
//	cout << " k = 1:" <<cols<<endl;
	// mul

	for (long k = 0 ; k < cols ; k++ )
		for (long i = 0 ; i < mat1.getrows() ; i++ )
			for (long j = 0 ; j < mat1.getcols()  ; j++ )
				res.data[i][k] += (mat2[j][k]==0? 0: mat1.data[i][j]); // (same as mat1.data[i][j] * mat2.data[j][k])
	return true;

}



