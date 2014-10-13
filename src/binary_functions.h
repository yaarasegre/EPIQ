/*
 * functions.h
 *
 *  Created on: Aug 7, 2013
 *      Author: yaara
 */

#ifndef BINARY_FUNCTIONS_H_
#define BINARY_FUNCTIONS_H_

#include <sys/time.h>
#include <iostream>
#include <string>
#include <iterator>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm> // sort
#include <random>
#include <limits>
#include <chrono>
#include <time.h>
#include <getopt.h>
#include "typedefs.h"
#include "binary_regression.h"
using namespace std;

#define SECS_PER_MIN 60
#define SECS_PER_HOUR 3600

struct _data {
	SNP2** snps;  // TODO: This matrix consumes a lot of memory. replace it with snpIndices
	vector<int> snpIndices;  // sparse matrix of SNPs. Cells are indices of locations where snps[i][j]==1.
							 // SNPs are separated by cells containing n.
	_snpid* snpIndicesMap;   // Maps where a new snp starts at snpIndices.
	SNP2** intSnps1;
	SNP2** intSnps2;
	double** intPheno;
	double* pheno;
	double* sqrt_abs_pheno;
//	double* sign_sqrt_pheno;
	double* sqrt_p;
	double* bin_assign;

	char*	chromosome;
	long long int* basePairPos;
	struct _buffers {
		// create buffers for lrt functions
		Mat* mat_3_3;
		Mat* mat_3_3b;
		Mat* mat_3_n;
		Mat* beta3;
		Mat* mat_4_4;
		Mat* mat_4_4b;
		Mat* mat_4_n;
		Mat* beta4;
		Mat* predicted_y;
		Mat* Y;
		SNP2* thrdCol;
		SNP2* ones;
	} buffers;
	_pairs_set candidatePairs;
	_validated_pairs validatedPairs;
};

struct _parms {
	string snpsFile;
	string outfile;
	string phenoFile;
	string mapFile;
	string thresholdsFile;
	string binEdgesFile;
	long n;
	double sqrt_n;
	long m;
	long L;
	int numBins;
	double *binEdges;
	double **bins_t;
	int	minIntersect;
	double t;
	double t_test_thresh;
	double chi2cutoff; //  significance cutoff for GLRT
	long long min_BP_distance;
	bool doPermutationTests;

	struct _test {
		bool validateImmediately;
		double beta;		// linear model coefficient
		double base_chi2cutoff; //  significance cutoff for base model (top marginal snps)
		bool debug; // for debug
		_snpid realSNP1; // for debug
		_snpid realSNP2; // for debug
		double p1; // for debug
		double p2; // for debug
		double realSNP1maf; // for debug
		double realSNP2maf; // for debug
		string realSNP1file; // for debug
		string realSNP2file; // for debug
		int currentIntPair; // for debug
		string FPRfile;
		bool quickRun;
		bool skipValidations;
		int realPairIntersect;
		int tests;
	} test;
};


// structs for sorting
struct _score {
    double score;
    _snpid snpid;
};

struct by_score_ascend {
    inline bool operator()(_score const &left, _score const &right) {
        return left.score < right.score;
    }
};

struct by_score_descend {
    inline bool operator()(_score const &left, _score const &right) {
        return left.score > right.score;
    }
};


// prototypes
void parseArgs(int argc, char* argv[], _parms &parms);
void init(_parms &parms, _data &data, double &lrt_score);
void quickInit(_parms &parms, _data &data, double &lrt_score);
void findSnpIndices(vector<int> &snpIndices, _snpid snpIndicesMap[], SNP2 **snps, _parms const &parms) ;
void restoreMat(vector<int> &indices, _parms const &parms);

SNP2** readSnps(string& fname, long n, long &numSnps);
void readPheno(_data &data, _parms &parms);
double** readCsv(string fname, long rows, long cols);
void readBinEdges(_data &data, _parms &parms) ;
bool readMap(_parms &parms, _data &data) ;

double** snpNormalize(SNP2** const snps, long n, long m);

void myscale(double* vec, long n);
double calcStd(double* vec, double mean, long n);

void printMat(SNP2 **mat, long n, long m);
void printMat(double **mat, long n, long m);

void oneIterationBins(long snpsTocheck, _parms &parms, _data &data,
		default_random_engine& generator, normal_distribution<double> normdist,
		double* r, double* ry1, double* ry2) ;

// functions to find pairs that pass score threshold

inline _hashed_pair hashPair(_snpid snp1, _snpid snp2) {
	uint64_t both = min(snp1, snp2);
	return (both<<32 | max(snp1, snp2));
}

inline _pair unhashPair(_hashed_pair both) {
	_pair p;
	p.first = both>>32;
	p.second = both & 0xFFFFFFFF;
	return p;
}

void findPairsBin(_parms &parms, _data &data, vector<_score> &sortedScores_a, vector<_score> &sortedScores_b,
		vector<_score> const &origScores_a, vector<_score> const &origScores_b, long m, double t_square);

void calcDotBins(vector<vector<_score>> &sortedScores_a, vector<vector<_score>> &sortedScores_b,
		vector<_score> &origScore1, vector<_score> &origScore2,
		long snpsTocheck, _parms &parms, _data &data,
		default_random_engine& generator, normal_distribution<double> normdist,
		double* r, double* ry1, double* ry2);

void randArr(long n, double arr[], default_random_engine& generator, normal_distribution<double> normdist) ;

/*
 * Multiple regression
 */
void validateAllPairs(_data &data, _parms &parms);
_validated_pair validatePair(_pair pair, _parms &parms, _data &data);

void genPheno(double pheno[], SNP2**snps, long n, long a, long b, double beta,
		default_random_engine& generator) ;
_pair chooseSnps(double pheno[], SNP2**snps, _parms &parms,  _data &data,
		default_random_engine& generator, double &score)  ;
_pair samplePairAndPheno(double pheno[], SNP2**snps,  double &score, _parms &parms,  _data &data,
		default_random_engine& generator);

/**
 * Test base model
 */
bool baseModel (SNP2** const snps, double* pheno,
		long n, long m,  _parms &parms, _data &data,_pair realPair, double baseLRT_thresh);
bool baseTestRealPair(set<_score, by_score_ascend> &best , SNP2** const snps, double* pheno,
		long n, _parms &parms, _data &data, _pair realPair, double baseLRT_thresh);
void simpleAll(double score[], SNP2** const snps, double* pheno, long n, long m);
double doSimpleReg(_snpid snpid, SNP2** const snps, double* pheno, long n);
void findBest(set<_score, by_score_ascend> &best, double scores[],  long m);

void cleanup(_parms &parms, _data &data);


void printSet(_pairs_set pairs);

/**
 * Inline functions
 */

/*
 * find mean of pheno where snp1*snp2 != 0,
 * Do it by going over the indices where snp1 !=0.
 * More efficient when maf(snp1) < maf(snp2).
 */
inline double calcMean(_data &data, _snpid snp1, _snpid snp2, _parms &parms, double &count) {
	int i1 = data.snpIndicesMap[snp1];
	double sum = 0;
	count = 0;
	while (data.snpIndices[i1] < parms.n) { // not end of SNP indices

		if (data.snps[snp2][data.snpIndices[i1]] != 0) {
			sum += data.pheno[data.snpIndices[i1]];
			count++;
		}
		i1++;
	}
	return sum/count;

}

/*
 * Calculate sample standard deviation.
 * Do it by going over the indices where snp1 !=0.
 * More efficient when maf(snp1) < maf(snp2).
 */
inline double sampleStdev(_data &data, _snpid snp1, _snpid snp2, _parms &parms, double mean, double count) {
	int i1 = data.snpIndicesMap[snp1];
	double sum = 0;
	while (data.snpIndices[i1] < parms.n) { // not end of SNP indices
		if (data.snps[snp2][data.snpIndices[i1]] != 0) {
			sum += pow(data.pheno[data.snpIndices[i1]] -mean, 2);
		}
		i1++;
	}
	return sqrt(sum/(count - 1));
}

/*
 * calculate one-sample t-test score for the test:
 * y_i*x_1i*x_2i ~ N(miu,sigma^2), H0: miu = 0, H1: miu != 0
 */
inline double t_test(_data &data, _parms &parms, _snpid snp1, _snpid snp2) {

	_snpid a;
	// find snp with lower maf (for faster t-test)
	if(data.sqrt_p[snp1] > data.sqrt_p[snp2]) {
		a = snp1;
		snp1 = snp2;
		snp2 = a;
	}
	double count;
	double mean = calcMean(data, snp1, snp2, parms, count);
	double stdev = sampleStdev(data, snp1, snp2, parms, mean, count);

	if (count == 0) {
		return 0;
	}
	if (count == 1) {
		return mean;
	}
	return (mean / (stdev / sqrt(count)));
}

/*
 * Insert candidate pairs that passed t-test
 */
inline void doInsert(_parms &parms, _data &data, _snpid snp1, _snpid snp2){

	// do t-test before insertion
	double t_score = t_test(data, parms, snp1, snp2);
	// TODO: different threshold for diff degrees of freedom
	if(t_score  >= parms.t_test_thresh || t_score  <= - parms.t_test_thresh) {

		if (parms.test.validateImmediately) {
			// check if found already
			int c = data.validatedPairs.count(hashPair(snp1,snp2));

			if (c == 0) {
				// Not found already!
				_validated_pair vpair = validatePair({snp1, snp2},  parms, data);

				if (vpair.LRTstat >= parms.chi2cutoff) {
					data.validatedPairs.insert({hashPair(snp1, snp2),vpair});
				}
			}
		}
		else {
			data.candidatePairs.insert(hashPair(snp1, snp2));
		}
	}
}

/*
 * Find zero-based bin number according to p.
 */
inline int  binNum(double p, _parms &parms) {

	if(parms.numBins == 1) {
		return 0;
	}

	//outliers
	if(p>parms.binEdges[parms.numBins] )
		return parms.numBins-1;
	if(p<parms.binEdges[0])
		return 0;

	// binary search
	int low = 0; int high = parms.numBins;
	int currBin;
	bool found = false;
	while (!found) {
		currBin = floor((low+high)/2);
		if (parms.binEdges[currBin] > p) {
			high = currBin;
		} else if  (parms.binEdges[currBin+1] < p) {
			low = currBin;
		} else {
			found = true;
		}
	}
	return currBin;
	//	return floor(p*parms.numBins);
}

/*
 * Return true if snps are neighbors
 */
inline bool areNeighbors(_data &data,  _parms parms, _snpid snp1,_snpid snp2) {
	return (data.chromosome[snp1] == data.chromosome[snp2] &&
			abs(data.basePairPos[snp1] - data.basePairPos[snp2]) < parms.min_BP_distance);
}

void printDuration(time_t start, time_t end) ;

#endif /* BINARY_FUNCTIONS_H_ */

//void oneIteration(_pairs_set& pairs, long snpsTocheck, _parms &parms, _data &data,
//		default_random_engine& generator, normal_distribution<double> normdist,
//		double* r, double* ry1, double* ry2);

//void findPairsA(_pairs_set &passedPairs, _parms &parms, _data &data, vector<_score>&, vector<_score>&,
//		vector<_score> const &origScores1, vector<_score> const &origScores2, long m);
//void findPairsB(_pairs_set &passedPairs, _parms &parms, _data &data, vector<_score>&, vector<_score>&,
//		vector<_score> const &origScores1, vector<_score> const &origScores2, long m);
//void findPairsC(_pairs_set &passedPairs, _parms &parms, _data &data, vector<_score>&, vector<_score>&,
//		vector<_score> const &origScores1, vector<_score> const &origScores2, long m);
//void findPairsD(_pairs_set &passedPairs, _parms &parms, _data &data, vector<_score>&, vector<_score>&,
//		vector<_score> const &origScores1, vector<_score> const &origScores2, long m);

