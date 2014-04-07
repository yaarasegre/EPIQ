//============================================================================
// Name        : projections.cpp
// Author      : Yaara Arkin
// Version     :
// Copyright   : 
// Description : EPIQ
//============================================================================
#include <iostream>
#include "binary_functions.h"
using namespace std;


int main(int argc, char* argv[])  {

	struct tm *current;
	time_t now, start, end;
	time(&start);
	current = localtime(&start);
	printf("%2.i:%2.i:%2.i: Started\n", current->tm_hour, current->tm_min, current->tm_sec);

	// paraeters: <snps file
	_parms parms;
	_data data;
	parseArgs(argc, argv, parms) ;

	// read snps file
	cout << endl<< "reading SNPs file " << parms.snpsFile << endl;
	data.snps = readSnps(parms.snpsFile, parms.n ,parms.m); // sets m!
	cout << "Found " << parms.m << " snps, " << parms.n << " samples." << endl;
	_snpid m = parms.m;
	long n = parms.n;

	double lrt_score;
	init(parms, data, lrt_score);

	srand (time(NULL)); // initialize random seed:
	default_random_engine generator;
	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	generator.seed(seed1);
	normal_distribution<double> normdist(0.0,1.0);


	/*
	 * do L iterations
	 */

	// create buffers
	double r[n];
	double ry1[n];
	double ry2[n];
	vector<_score> score1(m); vector<_score> score2(m);


//	if (parms.expectedFPRforPairsSet > 0)  // TODO maybe reserve
//		// reserve for expected # of pairs
//		pairs.reserve(round(parms.expectedFPRforPairsSet * (parms.m*(parms.m-1)/2)));

	cout << "Starting projections" << endl;

	for (long l = 0; l < parms.L; l++) {

		if (l%10==0) { // progress message
			time(&now); current = localtime(&now);
			long long int pairsFound = (parms.test.validateImmediately ?
					data.validatedPairs.size() : data.candidatePairs.size());
			printf("- %02d:%02d:%02d: Iteration %ld of %ld, %lld pairs found. ",
								current->tm_hour, current->tm_min, current->tm_sec, l, parms.L, pairsFound);
			printDuration(start, now);

//			printf("- %02d:%02d:%02d: Iteration %ld of %ld, at most %ld pairs found\n",
//					current->tm_hour, current->tm_min, current->tm_sec, l, parms.L, pairs.element_count());
			cout.flush();
		}

		// report all pairs that pass threshold t in one iteration
		oneIterationBins(m, parms, data, generator,  normdist, r, ry1, ry2);

	}


	/**
	 * ------------------------
	 * validate all pairs found
	 * ------------------------
	 */
	time(&now);
	current = localtime(&now);
	printf("%02d:%02d:%02d: ", current->tm_hour, current->tm_min, current->tm_sec);

	if (!parms.test.validateImmediately) {
		cout << "Found " << data.candidatePairs.size() << " pairs for validation using linear regression. " ;
		cout << endl <<  "Validating pairs..." << endl;
		validateAllPairs(data, parms);
	}

	// print pairs
	cout << "\nFound " << data.validatedPairs.size() << " pairs after validation." << endl;
	for (_validated_pairs::iterator it = data.validatedPairs.begin() ; it != data.validatedPairs.end(); ++it) {
		cout << '(' << it->second.snpid1 << "," << it->second.snpid2 << ") LRT score "
				<< it->second.LRTstat << endl ;
	}

	cout << "Cleaning up..." << endl;
	cleanup(parms, data);

	time(&end);
	current = localtime(&end);
	printf("%2.d:%2.d:%2.d: Done!\n", current->tm_hour, current->tm_min, current->tm_sec);

	printDuration(start, end);

	return 0;
}
