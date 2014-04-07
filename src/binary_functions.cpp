/*
 * functions.cpp
 *
 *  Created on: Aug 7, 2013
 *      Author: yaara
 */

#include "binary_functions.h"


void parseArgs(int argc, char* argv[], _parms &parms) {

	/**
	 * Default values
	 */
	parms.t = 0;
	parms.test.debug = false; //default
	parms.minIntersect = 0;
//	parms.seedOffset = 0;
	parms.test.quickRun = false;
	parms.test.skipValidations = false;
	parms.t_test_thresh = 3.8; // TODO get as input
	parms.numBins = 1;
	parms.phenoFile = "";
	parms.mapFile = "";
	parms.test.validateImmediately = true;
	parms.test.currentIntPair = 0;
	parms.min_BP_distance = 10000;
	parms.doPermutationTests = false;

	int nflag = 0;
//	int mflag = 0;
	int c;

	string str(256, ' ');

	while ((c = getopt (argc, argv, "n:m:f:p:o:M:e:D:L:t:N:F:i:C:b:1:2:T:B:W:X:Y:Z:dQVP")) != -1)
		switch (c)
		{
		case 'n':
			nflag = 1;
			parms.n = atoi(optarg);
			parms.sqrt_n = sqrt(parms.n);
			printf ("n = %s\n",optarg);

			break;
		case 'm':
			parms.m = atoi(optarg);
			printf ("m = %s\n",optarg);
			break;
		case 'L':
			parms.L = atoi(optarg);
			printf ("# iterations = %s\n",optarg);
			break;
		case 't':
			parms.t = atof(optarg);
			printf ("t = %s\n",optarg);
			break;
		case 'e':
			str.assign(optarg);
			parms.binEdgesFile.assign(str.c_str());
			printf ("Bins file = %s\n",optarg);
			break;
		case 'i':
			parms.minIntersect = atoi(optarg);
			printf ("Min intersect = %s\n",optarg);
			break;
		case 'C':
			parms.chi2cutoff = atof(optarg);
			printf ("LRT score cutoff = %s\n",optarg);
			break;
		case 'D':
			parms.min_BP_distance = atof(optarg);
			printf ("Minimal base-pair distance between SNPs = %s\n",optarg);
			break;
		case 'B':
			parms.test.base_chi2cutoff = atof(optarg);
			printf ("Base model LRT score cutoff = %s\n",optarg);
			break;
		case 'f':
			parms.snpsFile.assign(optarg);
			cout<<  "SNPs file = " <<parms.snpsFile << endl;
			break;
		case 'p':
			str.assign(optarg);
			parms.phenoFile.assign(str.c_str());
			cout<<  "Phenotype file = " <<parms.phenoFile << endl;
			break;
		case 'M':
			str.assign(optarg);
			parms.mapFile.assign(str.c_str());
			cout<<  "Map file = " <<parms.mapFile << endl;
			break;
		case 'F':
			str.assign(optarg);
			parms.thresholdsFile.assign(str.c_str());
			cout<<  "Bin's thresholds file = " <<parms.thresholdsFile << endl;
			break;
		case 'o':
			cout<<  "Output file = " << optarg << endl;
			parms.outfile.assign(optarg);
			break;
		case 'P':
			cout<<  "Performing permutation" << endl;
			parms.doPermutationTests = true;
			break;
		case 'b':
			parms.test.beta = atof(optarg);
			cout<<  "Beta = " <<parms.test.beta << endl;
			break;
		case '1':
			parms.test.p1 = atof(optarg);
			cout<<  "First interacting snp maf = " <<parms.test.p1 << endl;
			break;
		case '2':
			parms.test.p2 = atof(optarg);
			cout<<  "Second interacting snp maf = " <<parms.test.p2 << endl;
			break;
		case 'W':
			str.assign(optarg);
			parms.test.realSNP1file.assign(str.c_str());
			cout<<  "First interacting snps file = " <<parms.test.realSNP1file << endl;
			break;
		case 'X':
			parms.test.realSNP2file.assign(optarg);
			cout<<  "Second interacting snps file = " <<parms.test.realSNP2file << endl;
			break;
		case 'd':
			parms.test.debug = true;
			cout<<  "Running in debug mode!" << endl;
			break;
		case 'Q':
			parms.test.quickRun = true;
			cout<<  "Quick run! (only check interacting pair)" << endl;
			break;
		case 'V':
			parms.test.skipValidations = true;
			cout<<  "Skipping validations stage" << endl;
			break;
		case 'T':
			parms.test.tests = atoi(optarg);
			cout<<  "Number of tests = " << optarg << endl;
			break;
		case '?':
			if (isprint (optopt))
				cerr<< "Unknown option '-" << optopt << "'\n";
			else
				cerr << "invalid usage" << endl;
			exit(0);
		default:
			cerr << "invalid usage"<< endl;
			exit(0);
		}

	if (!nflag) {
		cerr << "Missing parameter -n"<< endl;
		exit(0);

	}// TODO finish - check that input is valid

	if (parms.test.validateImmediately) {
		cout << "Validating pairs immediately" << endl;
	}
	if (parms.t != 0 && parms.binEdgesFile.size()!=0) {
		cerr << "invalid combination of arguments, -t and -e colide." << endl;
		exit(0);
	}
}



/**
 * Do init stuff
 */
double** initBinThresholds(_parms &parms) {

	if (parms.t == 0) {
		// read t's from file
		return readCsv(parms.thresholdsFile, parms.numBins, parms.numBins);
	}
	else {
		// just 1 bin
		double ** mat = new double*[1];
		mat[0] = new double[1];
		mat[0][0] = parms.t;
		parms.numBins = 1;
		return mat;
	}
}

void init(_parms &parms, _data &data, double &lrt_score) {

	// Allocate stuff
	data.pheno = new double[parms.n];
	data.sqrt_abs_pheno = new double[parms.n];
	data.bin_assign = new double[parms.m];
	data.snpIndicesMap = new _snpid[parms.m];
	data.sqrt_p = new double[parms.m];

	// read map file
	readMap(parms, data);

	/*
	 * Buffers
	 */
	// create buffers for lrt functions
	data.buffers.mat_3_3 = new Mat(3,3); data.buffers.mat_3_3b = new Mat(3,3);
	data.buffers.mat_3_n = new Mat(3,parms.n); data.buffers.beta3 = new Mat(3,1);

	data.buffers.mat_4_4 = new Mat(4,4); data.buffers.mat_4_4b = new Mat(4,4);
	data.buffers.mat_4_n = new Mat(4,parms.n); data.buffers.beta4 = new Mat(4,1);
	data.buffers.predicted_y = new Mat(parms.n,1);
	// pheno
	data.buffers.Y = new Mat(parms.n,1);
	for (unsigned int i = 0; i < parms.n; i++) {
		data.buffers.Y->data[i][0] = data.pheno[i];
	}
	data.buffers.thrdCol = new SNP2[parms.n];
	data.buffers.ones = new SNP2[parms.n];
	fill_n(data.buffers.ones, parms.n, 1);


	/*
	 * Phenotype
	 */

	if (!parms.test.debug) {
		 // read phenotype
		cout << endl<< "reading phenotypes file " << parms.phenoFile << endl;
		readPheno(data, parms);
		for (int sample = 0; sample<parms.n; sample++){
			data.buffers.Y->data[sample][0] = data.pheno[sample];
		}
//		cout << "Pheno:" << endl;
//		for (int i = 0; i < parms.n; i++) {
//			cout << data.pheno[i] << endl;
//		}
	}
	else {
		if (parms.phenoFile == "") {
			//  choose pair and sample phenotype
			unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
			default_random_engine generator(seed1);
			samplePairAndPheno(data.pheno, data.snps, lrt_score, parms, data, generator);

			cout << "Interacting pair is "<< parms.test.realSNP1 << ", " << parms.test.realSNP2 <<
					". Intersect size is " << parms.test.realPairIntersect << endl;

		}
		else {
			parms.test.realSNP1 = 0;
			parms.test.realSNP2 = 1;
			// copy snps
			for (int sample = 0; sample<parms.n; sample++){
				data.snps[0][sample] = data.intSnps1[parms.test.currentIntPair][sample];
				data.snps[1][sample] = data.intSnps2[parms.test.currentIntPair][sample];
				data.pheno[sample] = data.intPheno[sample][parms.test.currentIntPair];
				data.buffers.Y->data[sample][0] = data.pheno[sample];
			}
			_validated_pair vp = validatePair({0,1},  parms, data);
			cout << "Pair LRT statistic: " << vp.LRTstat << ". cutoff is " << parms.chi2cutoff << endl;
			lrt_score = vp.LRTstat;
		}
	}

	/*
	 * BINS
	 */
	if (parms.t == 0) {
		// thresholds should be taken from file
		readBinEdges(data, parms);
	}
	parms.bins_t = initBinThresholds(parms);

	/*
	 * Calculate "MAFs" ( Pr[snp==1] ) for each SNP.
	 */
	for (long snp = 0; snp < parms.m; snp++) {
		double currp = 0;
		for (long i = 0; i < parms.n; i++) {
			currp += data.snps[snp][i];
		}
		currp/= double(parms.n);
		data.sqrt_p[snp] = sqrt(currp);

		data.bin_assign[snp] = binNum(currp, parms);
	}

	if(parms.test.debug) {
		parms.test.realSNP1maf = pow(data.sqrt_p[0],2);
		parms.test.realSNP2maf = pow(data.sqrt_p[1],2);

	}

	/*
	 * Calculate snp indices
	 */
	findSnpIndices(data.snpIndices, data.snpIndicesMap, data.snps,  parms) ;

	/*
	 * Calculate sqrt of pheno
	 */
	for (long sample = 0; sample < parms.n; sample++) {

		data.sqrt_abs_pheno[sample] = sqrt(abs(data.pheno[sample]));
//		data.sign_sqrt_pheno[sample] = copysign(1,data.pheno[sample]) * data.sqrt_abs_pheno[sample];
	}

}

/*
 * For multiple iterations in test_alg. assumes parms.test.debug == true
 */
void quickInit(_parms &parms, _data &data, double &lrt_score) {

	data.candidatePairs.clear();
	data.validatedPairs.clear();

	/*
	 * Phenotype
	 */

	if (parms.phenoFile == "") {
		//  choose pair and sample phenotype
		unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator(seed1);
		samplePairAndPheno(data.pheno, data.snps, lrt_score, parms, data, generator);

		cout << "Interacting pair is "<< parms.test.realSNP1 << ", " << parms.test.realSNP2 <<
				". Intersect size is " << parms.test.realPairIntersect << endl;

	}
	else {
		// read from a debug phenotype file
		parms.test.realSNP1 = 0;
		parms.test.realSNP2 = 1;
		// copy snps
		for (int sample = 0; sample<parms.n; sample++){
			data.snps[0][sample] = data.intSnps1[parms.test.currentIntPair][sample];
			data.snps[1][sample] = data.intSnps2[parms.test.currentIntPair][sample];
			data.pheno[sample] = data.intPheno[sample][parms.test.currentIntPair];
			data.buffers.Y->data[sample][0] = data.pheno[sample];
		}
		_validated_pair vp = validatePair({0,1},  parms, data);
		cout << "Pair LRT statistic: " << vp.LRTstat << ". cutoff is " << parms.chi2cutoff << endl;
		lrt_score = vp.LRTstat;
	}

	/*
	 * Calculate "MAFs" for first 2 SNP.
	 */
	for (long snp = 0; snp < 2; snp++) {
		double currp = 0;
		for (long i = 0; i < parms.n; i++) {
			currp += data.snps[snp][i];
		}
		currp/= double(parms.n);
		data.sqrt_p[snp] = sqrt(currp);

		data.bin_assign[snp] = binNum(currp, parms);
	}

	parms.test.realSNP1maf = pow(data.sqrt_p[0],2);
	parms.test.realSNP2maf = pow(data.sqrt_p[1],2);

	/*
	 * re-calculate snp indices
	 */
	findSnpIndices(data.snpIndices, data.snpIndicesMap, data.snps,  parms) ;

	/*
	 * Calculate sqrt of pheno
	 */
	for (long sample = 0; sample < parms.n; sample++) {

		data.sqrt_abs_pheno[sample] = sqrt(abs(data.pheno[sample]));
	}

}

	/*
	 * Init bloom filter parameters
	 */
//	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//	minstd_rand0 generator (seed);
//
//
//	parameters.projected_element_count    = round(parms.expectedFPRforPairsSet * (parms.m*(parms.m-1)/2));
//	parameters.false_positive_probability = 0.0005; // TODO: parameter
//	parameters.random_seed                = generator();
//
//	if (!parameters)
//	{
//		cerr << "Error - Invalid set of bloom filter parameters!" << std::endl;
//		exit (EXIT_FAILURE);;
//	}
//	parameters.compute_optimal_parameters();
//	cout << "Bloom filter number of hashes: " << parameters.optimal_parameters.number_of_hashes
//		  << " table size: " << parameters.optimal_parameters.table_size << endl;


/*
 * flatten snp matrix to a vector with indices of non-zero values.
 * snpIndices contains indices of non-zero values
 * snpIndicesMap is a nx1 vector that contains, for each SNP, it's starting location in snpIndices.
 */
void findSnpIndices(vector<int> &snpIndices, _snpid snpIndicesMap[], SNP2 **snps, _parms const &parms) {

	snpIndices.clear();
	for (_snpid snp = 0; snp < parms.m; snp++) {
		snpIndicesMap[snp] = snpIndices.size();
		for (int sample = 0; sample < parms.n; sample++) {
			if (snps[snp][sample] !=0) {
				snpIndices.push_back(sample);
			}
		}
		snpIndices.push_back(parms.n); // indicate end of snp
	}
}

/*
 * For debug
 */
void restoreMat(vector<int> &indices, _parms const &parms) {

	int sample = 0;
	for (vector<int>::iterator it = indices.begin() ; it != indices.end(); ++it) {
		if (*it == parms.n) {
			while(sample < parms.n ) {
				cout << "0 ";
				sample++;
			}
			cout << endl;
			sample = 0;
		}
		else {
			while(sample < *it ) {
				cout << "0 ";
				sample++;
			}
			cout << "1 ";
			sample++;
		}
	}
}




/*
 * Read SNPs file.
 * n is an input
 * numSnps is discovered from the file.
 * Snps stored in an **m x n** matrix, so each SNP is continuous in memory
 */
SNP2** readSnps(string& fname, long n, long &numSnps) {

	long filerows = 0;

	ifstream snpsfile(fname.data());
	string line;
	SNP2** snps = NULL;

	if (snpsfile.is_open()) {

		// get numSnps:
		snpsfile.seekg (0, snpsfile.end);
		numSnps = snpsfile.tellg()/n - 1;
		snpsfile.seekg(0, snpsfile.beg);

		// allocate matrix
		snps = new SNP2*[numSnps];
		if (snps)
		{
			snps[0] = new SNP2[numSnps*n]; // allocate all at once
			for (long i = 0; i < numSnps; ++i)
				snps[i] = snps[0] + i * n; // first element in row points at the row
		} else {
			cerr << "SNPs matrix allocation failed " << endl;
			exit(0);
		}

		while ( snpsfile.good() && filerows < n)
		{

			getline (snpsfile,line);

			if ((long)line.length()!= numSnps && filerows < n ) {
				cerr << "Invalid line length " << line.length() << " at line " <<
						filerows +1 << ". expected " << numSnps<<  endl;
				exit(0);
			} else {

				if(filerows % 500 == 499) { // DEBUG
					cout<< "reading line " << filerows +1 << endl;
				}
				// read each char and convert to SNP3
				for (long filecol = 0; filecol<numSnps; filecol++){

					// !!
					snps[filecol][filerows] = atoi(line.substr(filecol,1).c_str());

					if (snps[filecol][filerows] != 0 && snps[filecol][filerows]!= 1) {
						cerr << "Invalid SNP value " << line.substr(filecol,1) <<
								" at line " << filerows +1 << " filecol " << filecol << endl;
						exit(0);
					}

				}
			}
			filerows++;
		}
		snpsfile.close();
	}
	else {
		cerr << "Unable to open file " << fname << endl;
		exit(0);
	}

	return snps;
}


///**
// * read multiple phenotypes
// */
//double** readMultiPheno(_parms &parms, long cols) {
//
//	// allocate matrix
//	double ** multipheno = new double*[parms.n];
//	if (multipheno) {
//		multipheno[0] = new double[cols*parms.n]; // allocate all at once
//		for (long i = 0; i < parms.n; ++i)
//			multipheno[i] = multipheno[0] + i * cols; // first element in row points at the row
//	} else {
//		cerr << "Phenotypes matrix allocation failed " << endl;
//		exit(0);
//	}
//
//	long row = 0;
//
//	ifstream pfile(parms.phenoFile.data());
//	string line;
//
//	// read pheno
//	if (pfile.is_open()) {
//		while ( pfile.good() && row < parms.n)
//		{
//			getline (pfile,line);
//			string word;
//			stringstream stream(line,  ios_base::in);
//			int col = 0;
//			while( getline(stream, word, ',') ) {
//				multipheno[row][col] = atof(word.c_str());
//				col++;
//			}
//			if (col != cols) {
//					cerr << "ERROR: Expecting "<<cols<<" cols in phenotype file. Found " << col << " cols in row "
//							<< row+1 << " at " <<
//							parms.phenoFile << endl;
//					exit(0);
//				}
//			row++;
//		}
//		pfile.close();
//	}
//	else {
//		cerr << "Unable to open file " << parms.phenoFile << endl;
//		exit(0);
//	}
//
//	if (row != parms.n) {
//		cerr << "ERROR: Expecting "<<parms.n<<" rows in phenotype file. Found " << row << " rows in " <<
//				parms.phenoFile << endl;
//		exit(0);
//	}
//
//	return multipheno;
////	myscale(pheno, parms.n); assume scaled
//
//}

/**
 * read double csv
 */
double** readCsv(string fname, long rows, long cols) {

	// allocate matrix
	double ** mat = new double*[rows];
	if (mat) {
		mat[0] = new double[cols*rows]; // allocate all at once
		for (long i = 0; i < rows; ++i)
			mat[i] = mat[0] + i * cols; // first element in row points at the row
	} else {
		cerr << "Matrix allocation failed " << endl;
		exit(0);
	}
	long row = 0;

	ifstream file(fname.data());
	string line;

	// read file
	if (file.is_open()) {
		while ( file.good() && row < rows)
		{
			getline (file,line);
			string word;
			stringstream stream(line,  ios_base::in);
			int col = 0;
			while( getline(stream, word, ',') ) {
				mat[row][col] = atof(word.c_str());
				col++;
			}
			if (col != cols) {
					cerr << "ERROR: Expecting "<<cols<<" columns in file. Found " << col <<
							" columns in row " << row+1 << " at " <<fname << endl;
					exit(0);
				}
			row++;
		}
		file.close();
	}
	else {
		cerr << "Unable to open file " << fname << endl;
		exit(0);
	}
	if (row != rows) {
		cerr << "ERROR: Expecting "<<rows<<" rows. Found " << row << " rows in " <<
				fname << endl;
		exit(0);
	}
	return mat;
}


/**
 * read a map file with 3-4 columns:
 * chr snpid {morgans} Base-pair-position
 */
bool readMap(_parms &parms, _data &data) {

	data.chromosome = new char[parms.m];
	data.basePairPos = new long long int[parms.m];

	ifstream file(parms.mapFile.data());
	string line;
	_snpid row=0;

	// read file
	if (file.is_open()) {
		while (file.good() && row < parms.m)
		{
			getline (file,line);
			// split line
			istringstream iss(line);
			vector<string> tokens{istream_iterator<string>{iss},
			         istream_iterator<string>{}};

//			cout << " " << tokens[0] << " "<< tokens[1] << " "<< tokens[2] << " "<< tokens[3] << " " << endl;
			if (tokens.size() < 4) {
				cerr << "ERROR: Expecting 4 columns in map file. Found " << tokens.size() <<
						" columns in row " << row+1 << " at " <<parms.mapFile << endl;
				exit(0);
			}
			// get chr
			data.chromosome[row] = tokens[0].at(0);
			data.basePairPos[row] = atoi(tokens[3].c_str());
			row++;
		}
		file.close();
	}
	else {
		cerr << "Unable to open file " << parms.mapFile << endl;
		exit(0);
	}
	return true;
}



/**
 * read phenotypes
 */
void readPheno(_data &data, _parms &parms) {

	long row = 0;

	ifstream pfile(parms.phenoFile.data());
	string line;

	// read pheno
	if (pfile.is_open()) {
		while ( pfile.good() && row < parms.n)
		{
			getline (pfile,line);
			data.pheno[row] = atof(line.c_str());
			row++;
		}
		pfile.close();
	}
	else {
		cerr << "Unable to open file " << parms.phenoFile << endl;
		exit(0);
	}

	if (row < parms.n) {
		cerr << "ERROR: Expecting "<<parms.n<<" rows in phenotype file. Found only " << row << " rows in " <<
				parms.phenoFile << endl;
		exit(0);
	}
	myscale(data.pheno, parms.n);

}


/**
 * read phenotypes
 */
void readBinEdges(_data &data, _parms &parms) {

	cout << "Reading bin edges file " << parms.binEdgesFile << endl;
	ifstream file(parms.binEdgesFile.data());
	string line;
	vector<double> temp;
	int col = 0;
	// read file
	if (file.is_open()) {
		getline (file,line);
		string word;
		stringstream stream(line,  ios_base::in);
		while( getline(stream, word, ',') ) {
			temp.push_back(atof(word.c_str()));
			col++;
		}
		file.close();
	}
	else {
		cerr << "Unable to open file " << parms.binEdgesFile << endl;
		exit(0);
	}
	parms.numBins = col-1;
	parms.binEdges = new double[col];
	for (int i = 0; i<col; i++) {
		if (temp[i] < 0 || temp[i] > 1) {
			cerr << "Invalid bin edge " << temp[i] << " in " << parms.binEdgesFile<<endl;
			exit(0);
		}
		if (i>0 && temp[i] < temp[i-1]) {
			cerr << "Invalid bin edge " << temp[i] << ": edges should be in ascending order." <<
					" File name: " << parms.binEdgesFile<<endl;
			exit(0);
		}
		parms.binEdges[i] = temp[i];
	}

}

/**
 * Normalize a vector: V_i = (V_i-mean)/stdev
 */
void myscale(double vec[], long n) {
	double sum = 0;
	for (long row = 0; row < n; ++row)
		sum += vec[row];

	double mean = sum/(double)n;
	double std = calcStd(vec, mean, n);
	for (long i = 0; i < n; ++i) {
		vec[i]= (vec[i]-mean)/std;
	}
}

//double calcMeanCol(SNP2** const snps, long n, long m, long col) {
//	double sum = 0;
//	for (long row = 0; row < n; ++row) {
//		sum += snps[row][col];
//	}
//	return sum/(double)n;
//}

//double calcStdCol(SNP2** const snps, double mean, long n, long m, long col) {
//	double sum = 0;
//	for (long row = 0; row < n; ++row) {
//		sum += pow(snps[row][col] - mean, 2);
//	}
//	return sqrt(sum/(double)n);
//}

double calcStd(double* vec, double mean, long n) {
	double sum = 0;
	for (long row = 0; row < n; ++row) {
		sum += pow(vec[row] - mean, 2);
	}
	return sqrt(sum/(double)n);
}

void printMat(SNP2 **mat, long n, long m) {
	for (long i = 0; i<n; i++){
		for (long j = 0; j<m; j++){
			printf("%d",mat[i][j]);
		}
		cout << endl;
	}
}

void printMat(double **mat, long n, long m) {
	for (long i = 0; i<n; i++){
		for (long j = 0; j<m; j++){
			printf("%f ",mat[i][j]);
		}
		cout << endl;
	}
}


/**
 * Do one iteration.
 * 1. calculate dot product for all pairs.
 * 			outputs:
 * 			(a)	origScoresA = v*r, origScoresB = 1/(u*r)
 *			(b) for each bin, sorted vectors of positive/negative scores1/scores2
 * 2. for each pair of bins, report all pairs that pass a threshold (|v*r x u*r| > t),
 * 	  where t is a different threshold for each pair of bins.
 *
 * Parameter pairs contains result (reported pairs).
 * snpsTocheck = 2 or m, depending on running mode.
 */
void oneIterationBins(long snpsTocheck, _parms &parms, _data &data,
		default_random_engine& generator, normal_distribution<double> normdist,
		double* r, double* ry1, double* ry2) {

	vector<vector<_score>> sortedScores_a; sortedScores_a.resize(parms.numBins);
	vector<vector<_score>> sortedScores_b; sortedScores_b.resize(parms.numBins);
	vector<_score> origScoresA(snpsTocheck); vector<_score> origScoresB(snpsTocheck);


	calcDotBins(sortedScores_a, sortedScores_b, origScoresA, origScoresB, snpsTocheck,
			parms, data, generator, normdist, r, ry1, ry2);

	// sort results
	for (int i = 0; i<parms.numBins; i++) {
		sort(sortedScores_a[i].begin(), sortedScores_a[i].end(), by_score_ascend()); // sum(r*sqrt(y)*snp)^2
		sort(sortedScores_b[i].begin(), sortedScores_b[i].end(), by_score_ascend()); // 1/sum(r*sign(y)sqrt(y)*snp)^2
	}

	for (int bin1 = 0; bin1<parms.numBins; bin1++) {
		for (int bin2 = bin1; bin2<parms.numBins; bin2++) {

			/* find i,j s.t. a^2*b^2 >= (t)^2 */
			double t_square = pow(parms.bins_t[bin1][bin2],2);
			findPairsBin(parms, data, sortedScores_a[bin1], sortedScores_b[bin2],
					origScoresA, origScoresB, snpsTocheck,t_square);
		}
	}
}




/***********************************
 *           BINS
 ***********************************/
/**
 * Only if pair is found twice (a^2 > t^2/b'^2 & a'^2 > t^2/b^2) then add to pairs
 */
void findPairsBin(_parms &parms, _data &data, vector<_score> &sortedScores_a, vector<_score> &sortedScores_b,
		vector<_score> const &origScores_a, vector<_score> const &origScores_b, long m, double t_square) {

	if (sortedScores_a.size() == 0 || sortedScores_b.size() == 0)
		return;
	long i2 = 0;
	for (long i1 =  sortedScores_a.size()-1; i1  >=0 &&
			(sortedScores_a[i1].score >= t_square * sortedScores_b[i2].score); i1--) {

		for (i2 = 0;i2 < (long)sortedScores_b.size()
				&& sortedScores_a[i1].score >= t_square * sortedScores_b[i2].score; i2++) {
			_snpid snp1 = sortedScores_a[i1].snpid;
			_snpid snp2 = sortedScores_b[i2].snpid;

			// Check that SNPs are different and not neighbors
			if (snp1 != snp2 &&	!areNeighbors(data, parms, snp1, snp2)) {
				double a2 = origScores_a[snp2].score;
				double t_div_b1 = t_square*origScores_b[snp1].score;
				if ( a2>=t_div_b1)  {

					// (a > t/b & b > t/a), add to set:
					doInsert(parms, data, snp1, snp2);
				}
			}
		}
		i2 = 0;
	}
}





/**
 * calculate "dot product" of r, pheno and each snps
 */
void calcDotBins(vector<vector<_score>> &sortedScores_a, vector<vector<_score>> &sortedScores_b,
		vector<_score> &origScores_a, vector<_score> &origScores_b,
		long snpsTocheck, _parms &parms, _data &data,
		default_random_engine& generator, normal_distribution<double> normdist,
		double* r, double* ry1, double* ry2) {

	// sample r
	randArr(parms.n,r, generator, normdist);
	for (long row = 0; row < parms.n; row++) {
		ry1[row] = r[row] * data.sqrt_abs_pheno[row];
		ry2[row] = ry1[row]*copysign(1,data.pheno[row]) ; // r*+-sqrt(|y|)
	}

	vector<int>::iterator it = data.snpIndices.begin();
	// calculate "dot product" of r, pheno and all snps
	for (long snp = 0; snp < snpsTocheck; snp++) {

		double sum1 = 0;
		double sum2 = 0;

		while(*it != parms.n) {
			// found a non zero index
			sum1+=ry1[*it]; // dot product <ry1,snp>
			sum2+=ry2[*it];
			it++;
		}
		it++;

		/*
		 * Save results
		 */
		_score score1;

		score1.score =  pow(sum1/(data.sqrt_p[snp]*parms.sqrt_n) ,2) ; // a
		score1.snpid = snp;
		origScores_a[snp] = score1;

		_score score2;
		score2.score = pow(data.sqrt_p[snp] / sum2 , 2); //  1/b
		score2.snpid = snp;
		origScores_b[snp] = score2;

		int bin_num = data.bin_assign[snp];
//		if (bin_num >= parms.numBins) {
//			bin_num = parms.numBins-1;
//		}
//		if (bin_num < 0 ) {
//			bin_num = 0;
//		}
		sortedScores_a[bin_num].push_back(score1);
		sortedScores_b[bin_num].push_back(score2);

	}
}

/**
 * Generate random array for projections
 */
void randArr(long n, double arr[], default_random_engine& generator,
		normal_distribution<double> distribution) {

	for (long i = 0; i < n; i++)  {
		arr[i] = (distribution(generator));
//		arr[i] = copysign(1,distribution(generator));
	}

}
//
///**
// * Do linear regression on all pairs that pass filtering stage
// */
//void validateAllPairs(_validated_pairs &vpairs, _pairs_set &pairs, double chi2cutoff,
//		SNP2** const snps, double* pheno,long n) {
//
//	_validated_pair vpair ;
//
//	// create buffers for lrt functions
//	Mat mat_3_3(3,3);
//	Mat mat_3_3b(3,3);
//	Mat mat_3_n(3,n);
//	Mat beta3(3,1);
//
//	Mat mat_4_4(4,4);
//	Mat mat_4_4b(4,4);
//	Mat mat_4_n(4,n);
//	Mat beta4(4,1);
//
//	Mat predicted_y(n,1);
//	// pheno
//	Mat Y(n,1);
//	for (unsigned int i = 0; i < n; i++) {
//		Y(i,0) = pheno[i];
//	}
//	SNP2 thrdCol[n];
//	SNP2 ones[n];
//	fill_n(ones, n, 1);
//
////	clock_t begin = clock();  // time
//
//	struct tm *current;
//	time_t now;
//
//	long i = 0;
//	for (_pairs_set::iterator it=pairs.begin(); it!=pairs.end(); ++it) {
//
//		vpair = validatePair(*it, snps, n, Y, ones,
//				mat_3_3, mat_3_3b, mat_3_n, beta3, mat_4_4, mat_4_4b, mat_4_n, beta4,
//				predicted_y, thrdCol);
//		if (i%1000000 == 0) {
//			time(&now);
//			current = localtime(&now);
//			printf("%2.i:%2.i:%2.i ", current->tm_hour, current->tm_min, current->tm_sec);
//			cout << "validation " << i << " of " << pairs.size() <<  endl;
//		}
//		i++;
//		if (vpair.LRTstat >= chi2cutoff) {
//			vpairs.push_back(vpair);
//		}
//	}
//}
//

/**
 * Do linear regression on all pairs that pass filtering stage
 */
void validateAllPairs(_data &data, _parms &parms) {

	_validated_pair vpair;
	long long i = 0;
	struct tm *current; time_t now;

	for (auto it = data.candidatePairs.begin(); it != data.candidatePairs.end(); ++it ) {
		i++;
		_pair pair = unhashPair(*it);

//		double t_score = t_test(data, parms, pair.first, pair.second);
//		if(t_score  >= parms.t_test_thresh || t_score  <= - parms.t_test_thresh) {

		vpair = validatePair(pair, parms, data);

		if (vpair.LRTstat >= parms.chi2cutoff) {
			data.validatedPairs.insert({hashPair(pair.first,pair.second),vpair}); // TODO maybe report here and not later
		}
//		}

		if (i%10000 ==0 ) {

			time(&now);
			current = localtime(&now);
			printf("%2.d:%2.d:%2.d: validation ", current->tm_hour, current->tm_min, current->tm_sec);

			cout  <<  i << " of " << data.candidatePairs.size() << ". "<< endl;
		}

	}

//	for (_snpid snp1 = 0; snp1 < parms.m; snp1++) { //BLOOM
//		for (_snpid snp2 = snp1+1; snp2 < parms.m; snp2++) {
//			if(pairs.contains(snp1*parms.m + snp2)) {
//
//				vpair = validatePair({snp1,snp2}, snps, n, Y, ones,
//								mat_3_3, mat_3_3b, mat_3_n, beta3, mat_4_4, mat_4_4b, mat_4_n, beta4,
//								predicted_y, thrdCol);
//				i++;
//				if (i%100000 ==0 )
//					cout << i << " " << snp1 << "," << snp2 <<  endl;
//				if (vpair.LRTstat >= parms.chi2cutoff) {
//					vpairs.push_back(vpair);
//				}
//			}
//		}
//	}

}



/**
 * multiple regrssion on one pair
 */
_validated_pair validatePair(_pair pair, _parms &parms, _data &data) {

	_validated_pair vpair;
	vpair.snpid1 = pair.first;
	vpair.snpid2 = pair.second;

	// calculate 3rd col
	for (unsigned int i = 0; i < parms.n; i++) {
		data.buffers.thrdCol[i] = (data.snps[vpair.snpid1][i] && data.snps[vpair.snpid2][i] ? 1 : 0);
	}

	// create matrices
	SNP2* Xnull_tr[3] = {data.buffers.ones, data.snps[vpair.snpid1], data.snps[vpair.snpid2]};
//	Mat Xnull(n,2,nullCols);

	SNP2* Xalt_tr[4] = {data.buffers.ones, data.snps[vpair.snpid1], data.snps[vpair.snpid2], data.buffers.thrdCol};
//	Mat Xalt(n,3,altCols);

	vpair.LRTstat = multipleRegression(Xnull_tr, Xalt_tr , *data.buffers.Y, *data.buffers.mat_3_3,
			*data.buffers.mat_3_3b, *data.buffers.mat_3_n, *data.buffers.beta3,
			*data.buffers.mat_4_4, *data.buffers.mat_4_4b, *data.buffers.mat_4_n, *data.buffers.beta4,
			*data.buffers.predicted_y);

	// TODO calculate coefficients and their p-values

	return vpair;

}

void genPheno(double pheno[], SNP2**snps, long n, long a, long b, double beta,
		default_random_engine& generator) {

	normal_distribution<double> distribution(0.0,1.0);

	for (long i = 0; i < n; i++) {
		// y = a*b*beta + err
		pheno[i] = snps[a][i] * snps[b][i] * beta +
				distribution(generator);
	}
	myscale(pheno, n);
	//cout << "y " << distribution(generator);
}

/**
 * choose snps so that score is larger than threshold. FOR DEBUG
 */
_pair chooseSnps(double pheno[], SNP2**snps, _parms &parms, _data &data,
		default_random_engine& generator, double &score) {
	_pair realPair(0,0);
	do {
		// choose real pair & generate pheno
		parms.test.realSNP1 = rand() % parms.m;
		parms.test.realSNP2 = rand() % parms.m;

		// Check snps intersect
		double intersect = 0;
		for(int i = 0; i<parms.n; i++)
			intersect += snps[parms.test.realSNP1][i] * snps[parms.test.realSNP2][i];

		while ((parms.test.realSNP1 == parms.test.realSNP2) | (intersect < parms.minIntersect)) {
			parms.test.realSNP2 = rand() % parms.m;
			for(int i = 0; i< parms.n; i++)
				intersect += snps[parms.test.realSNP1][i] * snps[parms.test.realSNP2][i];
		}
		cout << "Chose snps " << parms.test.realSNP1 << ", " << parms.test.realSNP2 << endl;
		genPheno(pheno, snps,parms.n, parms.test.realSNP1, parms.test.realSNP2, parms.test.beta,  generator);

		for (unsigned int i = 0; i < parms.n; i++) {
			data.buffers.Y->data[i][0] = pheno[i];
		}

		// regress real pair
		realPair.first = min(parms.test.realSNP1,parms.test.realSNP2);
		realPair.second = max(parms.test.realSNP1,parms.test.realSNP2);

		_validated_pair vp = validatePair(realPair, parms, data);
		cout << "Pair LRT statistic: " << vp.LRTstat << ". cutoff is " << parms.chi2cutoff << endl;
		score = vp.LRTstat;

	} while(score < parms.chi2cutoff);
	return realPair;
}



/**
 * read interacting snps. FOR DEBUG
 */
_pair readIntSnpsAndPheno(double pheno[], SNP2**snps, _parms &parms,  _data &data,
		default_random_engine& generator, double &score) {
	_pair realPair(0,0);
	do {
		// choose real pair & generate pheno
		parms.test.realSNP1 = rand() % parms.m;
		parms.test.realSNP2 = rand() % parms.m;

		// Check snps intersect
		double intersect = 0;
		for(int i = 0; i<parms.n; i++)
			intersect += snps[parms.test.realSNP1][i] * snps[parms.test.realSNP2][i];

		while ((parms.test.realSNP1 == parms.test.realSNP2) | (intersect < parms.minIntersect)) {
			parms.test.realSNP2 = rand() % parms.m;
			for(int i = 0; i< parms.n; i++)
				intersect += snps[parms.test.realSNP1][i] * snps[parms.test.realSNP2][i];
		}
		cout << "Chose snps " << parms.test.realSNP1 << ", " << parms.test.realSNP2 << endl;
		genPheno(pheno, snps,parms.n, parms.test.realSNP1, parms.test.realSNP2, parms.test.beta,  generator);
		for (unsigned int i = 0; i < parms.n; i++) {
			data.buffers.Y->data[i][0] = pheno[i];
		}

		// regress real pair
		realPair.first = min(parms.test.realSNP1,parms.test.realSNP2);
		realPair.second = max(parms.test.realSNP1,parms.test.realSNP2);

		_validated_pair vp = validatePair(realPair, parms, data);
		cout << "Pair LRT statistic: " << vp.LRTstat << ". cutoff is " << parms.chi2cutoff << endl;
		score = vp.LRTstat;

	} while(score < parms.chi2cutoff);
	return realPair;
}


/**
 * Sample genotype of interacting SNP pair, generate phenotype so that y=beta*x1*x2 + err.
 * Sampled genotype is written to snp matrix, overriding current value.
 */
_pair samplePairAndPheno(double pheno[], SNP2**snps,  double &score, _parms &parms,  _data &data,
		default_random_engine& generator) {

	parms.test.realSNP1 = 0;
	parms.test.realSNP2 = 1;
	_pair realPair(parms.test.realSNP1,parms.test.realSNP2);
	uniform_real_distribution<double> distribution(0.0,1.0);

	// sample snp genotypes with minimal intersection
	int intersect;
	do {
		for (int i = 0; i<parms.n; i++) {
			snps[parms.test.realSNP1][i] = (distribution(generator) > parms.test.p1? 0 : 1);
			snps[parms.test.realSNP2][i] = (distribution(generator) > parms.test.p2? 0 : 1);
		}
		intersect = 0;
		for(int i = 0; i<parms.n; i++)
			intersect += snps[parms.test.realSNP1][i] * snps[parms.test.realSNP2][i];
	} while (intersect < parms.minIntersect);

	// sample pheno
	genPheno(pheno, snps,parms.n, parms.test.realSNP1, parms.test.realSNP2, parms.test.beta,  generator);
	for (unsigned int i = 0; i < parms.n; i++) {
		data.buffers.Y->data[i][0] = pheno[i];
	}

	_validated_pair vp = validatePair(realPair, parms, data);
	cout << "Pair LRT statistic: " << vp.LRTstat ;
	cout.flush();
	cout << ". cutoff is " << parms.chi2cutoff << endl;
	cout.flush();
	score = vp.LRTstat;

	parms.test.realPairIntersect = intersect;
	return realPair;
}
/**
 * Run base model: find top sqrt(2m) pairs, test all pairs between them
 */
bool baseModel (SNP2** const snps, double* pheno,
		long n, long m, _parms &parms, _data &data, _pair realPair, double baseLRT_thresh) {

	/*
	 * do simple regression for all snps
	 */
	double* scores = new double[parms.m];
	simpleAll(scores, snps, pheno, n,  m);

	set<_score, by_score_ascend> best;
	findBest(best, scores,  m);
	delete [] scores;
	/*
	 * Check if real pair is in group "best" and it's score is > threshold
	 */
	bool found = baseTestRealPair(best, snps, pheno, n, parms, data, realPair, baseLRT_thresh);
	return found;

}

/**
 * Check if real pair is in group "best" and it's score is > threshold
 */
bool baseTestRealPair(set<_score, by_score_ascend> &best , SNP2** const snps, double* pheno,
		long n, _parms &parms, _data &data, _pair realPair, double baseLRT_thresh) {

	// set's keys are scores and not snpids, so have to go over it...
	bool found1 = false;
	bool found2 = false;
	set<_score, by_score_ascend>::iterator it;
	for (it=best.begin(); it!=best.end() && (!found1 || !found2); ++it) {
//		cout << "best " << it->snpid << ":"<< it->score << endl;
		if (it->snpid == realPair.first) {
			found1 = true;
//			cout << "snp " << realPair.first<<"!\n";
		}
		if (it->snpid == realPair.second) {
			found2 = true;
//			cout << "snp " << realPair.second << "!\n";
		}
	}
	if (!found1 || !found2)
		return false; // pair not found!

	_validated_pair vp = validatePair(realPair, parms, data);
	return (vp.LRTstat >= baseLRT_thresh); // score > threshold?

}

/**
 * Do simple regression on all snps
 * TODO: remove new
 */
void simpleAll(double scores[], SNP2** const snps, double* pheno, long n, long m) {

	// buffers
	Mat mat_2_2(2,2);
	Mat mat_2_2b(2,2);
	Mat mat_2_n(2,n);
	Mat beta(2,2);
	Mat predicted_y(n,1);

	SNP2 ones[n];
	fill_n(ones, n, 1);

	// pheno
	Mat Y(n,1);
	for (unsigned int i = 0; i < n; i++) {
		Y(i,0) = pheno[i];
	}

	// create matrix
	SNP2* X_tr[2] ;
	for (_snpid snpid = 0; snpid<m; snpid++) {
		X_tr[0] = ones;
		X_tr[1] = snps[snpid];
		scores[snpid] = simpleRegression(X_tr,Y,  mat_2_2, mat_2_2b, mat_2_n, beta,predicted_y);
//		cout << snpid << "-" << scores[snpid] << endl;
	}
}


/**
 * Do simple regression on 1 snp
 */
double doSimpleReg(_snpid snpid, SNP2** const snps, double* pheno, long n) {

	// create matrices
	Mat Y(n,1);
	for (unsigned int i = 0; i < n; i++) {
		Y(i,0) = pheno[i];
	}

	SNP2 ones[n];
	fill_n(ones, n, 1);

	SNP2* X_tr[2]  = {ones, snps[snpid]};
	return simpleRegression(X_tr,Y);
}


/**
 * find *highest* sqrt(2m) scores. Used for base model.
 */
void findBest(set<_score, by_score_ascend> &best, double scores[],  long m) {

	unsigned long numresults = round(sqrt(2*m));

	// first, find threshold
	for (long i = 0; i < m; i++) {
		_score s;
		if (best.size() < numresults) {
			// add current to fill "top" group
			s.score = scores[i];
			s.snpid = i;
			best.insert(s);
		}
		else {

			set<_score, by_score_ascend>::iterator it = best.begin(); // largest value
//			cout << it->score;
			if (scores[i] > it->score) {
				// insert new score instead
				best.erase(it);
				s.score = scores[i];
				s.snpid = i;
				best.insert(s);
			}
		}
	}
}

void cleanup(_parms &parms, _data &data){

	delete [] parms.bins_t[0];
	delete [] parms.bins_t;

	delete [] data.snps[0];
	delete [] data.snps;

	delete [] data.pheno;
	delete [] data.sqrt_abs_pheno;
	delete [] data.bin_assign;
	delete [] data.snpIndicesMap;
	delete [] data.sqrt_p;

	delete data.buffers.mat_3_3;
	delete data.buffers.mat_3_3b;
	delete data.buffers.mat_3_n;
	delete data.buffers.beta3;
	delete data.buffers.mat_4_4;
	delete data.buffers.mat_4_4b;
	delete data.buffers.mat_4_n;
	delete data.buffers.beta4;
	delete data.buffers.predicted_y;
	delete data.buffers.Y;

	delete [] data.buffers.thrdCol;
	delete [] data.buffers.ones;
	if (parms.binEdges)
		delete [] parms.binEdges;

}

void printDuration(time_t start, time_t end) {
	// calculate duration
	unsigned long long int seconds;
	unsigned int minutes, hours, secs_left, mins_left;
	seconds = difftime(end,start);
	hours = seconds / SECS_PER_HOUR;
	minutes = seconds / SECS_PER_MIN;
	mins_left = minutes % SECS_PER_MIN;
	secs_left = seconds % SECS_PER_MIN;
	printf("Duration: %02d:%02d:%02d\n", hours, mins_left, secs_left);
}

void shufflePheno(_data &data, _parms &parms, default_random_engine &generator) {

	// copy phenotype to vector
	vector<double> perm(parms.n);
	perm.assign(data.pheno, data.pheno + parms.n);
	shuffle (perm.begin(), perm.end(), generator);

	// copy back
	for (int i = 0; i < parms.n; ++i) {
		data.pheno[i] = perm[i];
	}
}

//
///**
// * calculate "dot product" of r, pheno and all snps
// * new
// */
//void calcDot(vector<_score> &scoreAneg, vector<_score> &scoreApos,
//		vector<_score> &scoreBneg, vector<_score> &scoreBpos,
//		vector<_score> &origScoresA, vector<_score> &origScoresB,
//		SNP2** const snps, double* pheno, long snpsTocheck, _parms &parms, double sqrt_p[],
//		default_random_engine& generator, normal_distribution<double> normdist,
//		double* r, double* ry1, double* ry2, vector<int> &snpIndices) {
//
//
//
//	// sample r // TODO: r~N(0,1)
//	randArr(parms.n,r, generator, normdist);
//	for (long row = 0; row < parms.n; row++) {
//		ry1[row] = r[row]*sqrt(abs(pheno[row]));
//				ry2[row] = ry1[row]*copysign(1,pheno[row]) ; // r*+-sqrt(|y|)
////		ry1[row] = r[row];
////		ry2[row] = r[row]*pheno[row] ; // r*+-sqrt(|y|)
//	}
//
//	vector<int>::iterator it = snpIndices.begin();
//	// calculate "dot product" of r, pheno and all snps
//	for (long snp = 0; snp < snpsTocheck; snp++) {
//
//		double sum1 = 0;
//		double sum2 = 0;
//
//		while(*it != parms.n) {
//			// found a non zero index
//			sum1+=ry1[*it]; // dot product <ry1,snp>
//			sum2+=ry2[*it];
//			it++;
//			if (it == snpIndices.end()) {
//				cerr << "reached end of snps too soon!" << endl;
//				exit(0);
//			}
//		}
//		it++;
//
////		double oldsum1 = 0;
////		double oldsum2 = 0;
////		for (long sample = 0; sample < parms.n; sample++) {
////			// snps[snp][sample] is in {0/1}
////			if (snps[snp][sample]) {
//////				cout << sample << " ";
////				oldsum1+=ry1[sample]; // dot product <ry1,snp>
////				oldsum2+=ry2[sample];
////			}
////		}
////		if (oldsum1 != sum1 || oldsum2 != sum2) {
////			cerr << "!!!";
////			exit(0);
////		}
//
//
//		/*
//		 * Save results
//		 */
//		_score score1;
//
//		score1.score =  sum1/(sqrt_p[snp]*parms.sqrt_n) ;
//		score1.snpid = snp;
//		origScoresA[snp] = score1;
//
//		_score score2;
//		score2.score = parms.t * sqrt_p[snp] / sum2;
//		score2.snpid = snp;
//		origScoresB[snp] = score2;
//
//		// TODO: check if saving scores in a vec, sorting them and then creating scoreApos/neg
//		// helps performance. (might reduce branch prediction fails), probably not...
//		if (score1.score >= 0)
//			scoreApos.push_back(score1);
//		else
//			scoreAneg.push_back(score1);
//		if (score2.score >= 0)
//			scoreBpos.push_back(score2);
//		else
//			scoreBneg.push_back(score2);
//
//	}
//}

//
///**
// * Do one iteration. report all pairs that pass threshold t using
// * Parameter pairs contains result.
// * snpsTocheck = 2 or m, depending on running mode.
// */
//void oneIteration(_pairs_set& pairs, long snpsTocheck, _parms &parms, _data &data,
//		default_random_engine& generator, normal_distribution<double> normdist,
//		double* r, double* ry1, double* ry2) {
//
////	struct tm *current;
////	time_t now;
////	time(&now);
////	current = localtime(&now);
////	printf("%2.d:%2.d:%2.d: Started one iteration\n", current->tm_hour, current->tm_min, current->tm_sec);
//
//
//	vector<_score> scoreAneg; scoreAneg.reserve(snpsTocheck);
//	vector<_score> scoreApos; scoreApos.reserve(snpsTocheck); // TODO maybe get as parameter
//	vector<_score> scoreBneg; scoreBneg.reserve(snpsTocheck);
//	vector<_score> scoreBpos; scoreBpos.reserve(snpsTocheck);
//	vector<_score> origScoresA(snpsTocheck); vector<_score> origScoresB(snpsTocheck);
//
////	time(&now);
////	current = localtime(&now);
////	printf("%2.d:%2.d:%2.d: After allocations\n", current->tm_hour, current->tm_min, current->tm_sec);
//
////	int timedot1=time(NULL);
//	calcDot(scoreAneg, scoreApos, scoreBneg, scoreBpos, origScoresA, origScoresB, data.snps, data.pheno, snpsTocheck,
//			parms, data.sqrt_p, generator, normdist, r, ry1, ry2, data.snpIndices);
////	int timedot2=time(NULL);
////	cout << "Calc dot " << timedot2-timedot1 << endl;
//
////	time(&now);
////	current = localtime(&now);
////	printf("%2.d:%2.d:%2.d: After calcDot\n", current->tm_hour, current->tm_min, current->tm_sec);
//
//	sort(scoreAneg.begin(), scoreAneg.end(), by_score_ascend()); // sum(r*sqrt(y)*snp)
//	sort(scoreApos.begin(), scoreApos.end(), by_score_ascend()); // sum(r*sqrt(y)*snp)
//	sort(scoreBneg.begin(), scoreBneg.end(), by_score_ascend()); // t/sum(r*sign(y)sqrt(y)*snp)
//	sort(scoreBpos.begin(), scoreBpos.end(), by_score_ascend()); // t/sum(r*sign(y)sqrt(y)*snp)
//
////	time(&now);
////	current = localtime(&now);
////	printf("%2.d:%2.d:%2.d: After sort\n", current->tm_hour, current->tm_min, current->tm_sec);
//
////	int timefind1=time(NULL);
//
//	/* find i,j s.t. sum1*sum2 >= t*n  or sum1*sum2 <= -t*n  */
//	findPairsA(pairs, parms, data, scoreApos, scoreBpos, origScoresA, origScoresB, snpsTocheck);
//	findPairsB(pairs, parms, data, scoreAneg, scoreBneg, origScoresA, origScoresB, snpsTocheck);
//	findPairsC(pairs, parms, data, scoreApos, scoreBneg, origScoresA, origScoresB, snpsTocheck);
//	findPairsD(pairs, parms, data, scoreAneg, scoreBpos, origScoresA, origScoresB, snpsTocheck);
//
////	int timefind2=time(NULL);
////	cout << "find " << timefind2 - timefind1 << endl;
////	time(&now);
////	current = localtime(&now);
////	printf("%2.d:%2.d:%2.d: After find pairs\n\n", current->tm_hour, current->tm_min, current->tm_sec);
//
//}
//
//
////inline void doInsert(_pairs_set &passedPairs, _snpid snp1, _snpid snp2){
////	_snpid minid = min(snp1, snp2);
////	_snpid maxid = max(snp1, snp2);
////	passedPairs.insert({minid, maxid});
//////	uint64_t key = minid;
//////	key = key << 32 | maxid;
//////	passedPairs.insert({key, newpair});
////}
//
///**
// * A: a+*b+ > t => a+ > t/b+
// * Only if pair is found twice (a > t/b & b > t/a) then add to pairs
// */
//void findPairsA(_pairs_set &passedPairs, _parms &parms, _data &data,
//		vector<_score> &scoreApos, vector<_score> &scoreBpos,
//		vector<_score> const &origScoresA, vector<_score> const &origScoresB, long m) {
//
//	if (scoreApos.size() == 0 || scoreBpos.size() == 0)
//		return;
//	long i2 = 0;
//	for (long i1 =  scoreApos.size()-1; i1  >=0 &&  (scoreApos[i1].score >= scoreBpos[i2].score); i1--) {
//
//		for (i2 = 0;i2 < (long)scoreBpos.size() && scoreApos[i1].score >= scoreBpos[i2].score; i2++) {
//
//			_snpid snp1 = scoreApos[i1].snpid;
//			_snpid snp2 = scoreBpos[i2].snpid;
//			double s1 = origScoresA[snp2].score;
//			double s2 = origScoresB[snp1].score;
//			if (snp1 != snp2 &&
//					( (s1>=0 && s2>=0 && s1>=s2)  ||  (s1<=0 && s2<=0 && s1<=s2)) ) { // TODO: check LD block
//				// (a > t/b & b > t/a)
//				/* add to set */
//				doInsert(passedPairs, parms, data, snp1, snp2);
//			}
//		}
//		i2 = 0;
//	}
//}
//
///**
// * B: a-*b- > t => a- < t/b-
// * Only if pair is found twice (a > t/b & b > t/a) then add to pairs
// */
//void findPairsB(_pairs_set &passedPairs, _parms &parms, _data &data,
//		vector<_score> &scoreAneg, vector<_score> &scoreBneg,
//		vector<_score> const &origScoresA, vector<_score> const &origScoresB, long m) {
//
//	if (scoreAneg.size() == 0 || scoreBneg.size() == 0) {
//		return;
//	}
//	// Go over snps in scoreApos
//	long i2 = scoreBneg.size()-1;
//	for (long i1 = 0; i1 < (long)scoreAneg.size() && (scoreAneg[i1].score <= scoreBneg[i2].score); i1++) {
//		for (i2 = scoreBneg.size()-1; i2 >= 0 && scoreAneg[i1].score <= scoreBneg[i2].score; i2--) {
//
//			_snpid snp1 = scoreAneg[i1].snpid;
//			_snpid snp2 = scoreBneg[i2].snpid;
//
//			double s1 = origScoresA[snp2].score;
//			double s2 = origScoresB[snp1].score;
//
//			if (snp1 != snp2 &&
//					( (s1>=0 && s2>=0 && s1>=s2)  ||  (s1<=0 && s2<=0 && s1<=s2)) ){ // TODO: check LD block
//				//(a > t/b & b > t/a)
//				/* add to set */
//				doInsert(passedPairs, parms, data, snp1, snp2);
//			}
//		}
//		i2 = scoreBneg.size()-1;
//	}
//}
//
///**
// * C: a+*b- < -t => a+ > -t/b-
// * Only if pair is found twice (a < -t/b & b < -t/a) then add to pairs
// */
//void findPairsC(_pairs_set &passedPairs, _parms &parms, _data &data,
//		vector<_score> &scoreApos, vector<_score> &scoreBneg,
//		vector<_score> const &origScoresA, vector<_score> const &origScoresB, long m) {
//
//	if (scoreApos.size() == 0 || scoreBneg.size() == 0)
//		return;
//
//	// Go over snps in scoreApos
//	long i2 = scoreBneg.size()-1;
//	for (long i1 = scoreApos.size()-1; i1 >= 0 && (scoreApos[i1].score >= -scoreBneg[i2].score); i1--) {
//		for (i2 = scoreBneg.size()-1; i2 >= 0 && scoreApos[i1].score >= -scoreBneg[i2].score; i2--) {
//
//			_snpid snp1 = scoreApos[i1].snpid;
//			_snpid snp2 = scoreBneg[i2].snpid;
//
//			double s1 = origScoresA[snp2].score;
//			double s2 = origScoresB[snp1].score;
//			if (snp1 != snp2 &&
//					( (s1>=0 && s2<=0 && s1 >= -s2)  ||  (s1<=0 && s2>=0 && s1 <= -s2)) ) {
//				// (a < -t/b & b < -t/a)
//
//				/* add to set */
//				doInsert(passedPairs, parms, data, snp1, snp2);
//			}
//		}
//		i2 = scoreBneg.size()-1;
//	}
//}
//
///**
// * D: a-*b+ < -t => a- < -t/b+
// * Only if pair is found twice (a > t/b & b > t/a) then add to pairs
// */
//void findPairsD(_pairs_set &passedPairs, _parms &parms, _data &data,
//		vector<_score> &scoreAneg, vector<_score> &scoreBpos,
//		vector<_score> const &origScoresA, vector<_score> const &origScoresB, long m) {
//
//	if (scoreAneg.size() == 0 || scoreBpos.size() == 0)
//		return;
//
//	// Go over snps in scoreApos
//	long i2 = 0;
//	for (long i1 = 0; i1 < (long)scoreAneg.size() && (scoreAneg[i1].score <= -scoreBpos[i2].score); i1++) {
//		for (i2 = 0; i2 < (long)scoreBpos.size() && scoreAneg[i1].score <= -scoreBpos[i2].score; i2++) {
//
//			_snpid snp1 = scoreAneg[i1].snpid;
//			_snpid snp2 = scoreBpos[i2].snpid;
//
//			double s1 = origScoresA[snp2].score;
//			double s2 = origScoresB[snp1].score;
//
//			if (snp1 != snp2 &&
//					( (s1>=0 && s2<=0 && s1 >= -s2)  ||  (s1<=0 && s2>=0 && s1 <= -s2)) ) {
//				// (a < -t/b & b < -t/a)
//
//				/* add to set */
//				doInsert(passedPairs, parms, data, snp1, snp2);
//			}
//		}
//		i2 = 0;
//	}
//}
//



///**
// * a+*b+ > t => a+ > t/b+
// * Only if pair is found twice (a > t/b & b > t/a) then add to pairs
// */
//void findPairsABin(_parms &parms, _data &data, vector<_score> &scoreApos, vector<_score> &scoreBpos,
//		vector<_score> const &origScoresA, vector<_score> const &origScoresB, long m, double t) {
//
//	if (scoreApos.size() == 0 || scoreBpos.size() == 0)
//		return;
//	long i2 = 0;
//	for (long i1 =  scoreApos.size()-1; i1  >=0 && (scoreApos[i1].score >= t*scoreBpos[i2].score); i1--) {
//
//		for (i2 = 0;i2 < (long)scoreBpos.size() && scoreApos[i1].score >= t*scoreBpos[i2].score; i2++) {
//
//			_snpid snp1 = scoreApos[i1].snpid;
//			_snpid snp2 = scoreBpos[i2].snpid;
//			double Sa = origScoresA[snp2].score;
//			double Sb = t*origScoresB[snp1].score;
//			if (snp1 != snp2 &&
//					( (Sa>=0 && Sb>=0 && Sa>=Sb)  ||  (Sa<=0 && Sb<=0 && Sa<=Sb)) ) {
//				// (a > t/b & b > t/a)
//				/* add to set */
//				doInsert(parms, data, snp1, snp2);
//			}
//		}
//		i2 = 0;
//	}
//}
//
//
//
///**
// * B: a-*b- > t => a- < t/b-
// * Only if pair is found twice (a > t/b & b > t/a) then add to pairs
// */
//void findPairsBBin(_parms &parms, _data &data, vector<_score> &scoreAneg, vector<_score> &scoreBneg,
//		vector<_score> const &origScoresA, vector<_score> const &origScoresB, long m, double t) {
//
//	if (scoreAneg.size() == 0 || scoreBneg.size() == 0) {
//		return;
//	}
//	// Go over snps in scoreApos
//	long i2 = scoreBneg.size()-1;
//	for (long i1 = 0; i1 < (long)scoreAneg.size()
//			&& (scoreAneg[i1].score <= t*scoreBneg[i2].score); i1++) {
//		for (i2 = scoreBneg.size()-1; i2 >= 0 &&
//				t*scoreAneg[i1].score <= scoreBneg[i2].score; i2--) {
//
//			_snpid snp1 = scoreAneg[i1].snpid;
//			_snpid snp2 = scoreBneg[i2].snpid;
//
//			double s1 = origScoresA[snp2].score;
//			double s2 = t*origScoresB[snp1].score;
//
//			if (snp1 != snp2 &&
//					( (s1>=0 && s2>=0 && s1>=s2)  ||  (s1<=0 && s2<=0 && s1<=s2)) ){ // TODO: check LD block
//				//(a > t/b & b > t/a)
//				/* add to set */
//				doInsert(parms, data, snp1, snp2);
//			}
//		}
//		i2 = scoreBneg.size()-1;
//	}
//}
//
///**
// * C: a+*b- < -t => a+ > -t/b-
// * Only if pair is found twice (a < -t/b & b < -t/a) then add to pairs
// */
//void findPairsCBin(_parms &parms, _data &data, vector<_score> &scoreApos, vector<_score> &scoreBneg,
//		vector<_score> const &origScoresA, vector<_score> const &origScoresB, long m, double t) {
//
//	if (scoreApos.size() == 0 || scoreBneg.size() == 0)
//		return;
//
//	// Go over snps in scoreApos
//	long i2 = scoreBneg.size()-1;
//	for (long i1 = scoreApos.size()-1; i1 >= 0 &&
//			(scoreApos[i1].score >= -t*scoreBneg[i2].score); i1--) {
//		for (i2 = scoreBneg.size()-1; i2 >= 0 &&
//				scoreApos[i1].score >= -t*scoreBneg[i2].score; i2--) {
//
//			_snpid snp1 = scoreApos[i1].snpid;
//			_snpid snp2 = scoreBneg[i2].snpid;
//
//			double s1 = origScoresA[snp2].score;
//			double s2 = t*origScoresB[snp1].score;
//			if (snp1 != snp2 &&
//					( (s1>=0 && s2<=0 && s1 >= -s2)  ||  (s1<=0 && s2>=0 && s1 <= -s2)) ) {
//				// (a < -t/b & b < -t/a)
//
//				/* add to set */
//				doInsert(parms, data, snp1, snp2);
//			}
//		}
//		i2 = scoreBneg.size()-1;
//	}
//}
//
///**
// * D: a-*b+ < -t => a- < -t/b+
// * Only if pair is found twice (a > t/b & b > t/a) then add to pairs
// */
//void findPairsDBin(_parms &parms, _data &data, vector<_score> &scoreAneg, vector<_score> &scoreBpos,
//		vector<_score> const &origScoresA, vector<_score> const &origScoresB, long m, double t) {
//
//	if (scoreAneg.size() == 0 || scoreBpos.size() == 0)
//		return;
//
//	// Go over snps in scoreApos
//	long i2 = 0;
//	for (long i1 = 0; i1 < (long)scoreAneg.size()
//			&& (scoreAneg[i1].score <= -t*scoreBpos[i2].score); i1++) {
//		for (i2 = 0; i2 < (long)scoreBpos.size() &&
//				scoreAneg[i1].score <= -t*scoreBpos[i2].score; i2++) {
//
//			_snpid snp1 = scoreAneg[i1].snpid;
//			_snpid snp2 = scoreBpos[i2].snpid;
//
//			double s1 = origScoresA[snp2].score;
//			double s2 = t*origScoresB[snp1].score;
//
//			if (snp1 != snp2 &&
//					( (s1>=0 && s2<=0 && s1 >= -s2)  ||  (s1<=0 && s2>=0 && s1 <= -s2)) ) {
//				// (a < -t/b & b < -t/a)
//
//				/* add to set */
//				doInsert(parms, data, snp1, snp2);
//			}
//		}
//		i2 = 0;
//	}
//}
//
//
