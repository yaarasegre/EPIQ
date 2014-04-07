/*
 * typedefs.h
 *
 *  Created on: Sep 3, 2013
 *      Author: yaara
 */

#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
//#include "bloom_filter.hpp" // BLOOM

//#include <sparsehash/sparse_hash_set>
//using google::sparse_hash_set;      // namespace where class lives by default
//using ext::hash;  // or __gnu_cxx::hash, or maybe tr1::hash, depending on your OS

using namespace std;

//typedefs
typedef unsigned char SNP2; // binary snps, 0/1/2
typedef uint32_t _snpid;
typedef pair<_snpid,_snpid> _pair;
typedef uint64_t _hashed_pair;
typedef unordered_set<uint64_t> _pairs_set;

// struct for pairs
struct _validated_pair {

	_snpid snpid1;
	_snpid snpid2;
	double coeff1;
	double coeff2;
	double coeff3;
	double pval1;
	double pval2;
	double pval3;
	double LRTstat;
};
typedef  unordered_map<uint64_t, _validated_pair> _validated_pairs;


//typedef bloom_filter _pairs_set;  // BLOOM

// Do annoying stuff to make the hash work
//struct PairHasher
//{
//	inline std::size_t operator()(const _pair& key) const
//	{
//		using std::size_t;
//		using std::hash;
//		using std::string;
//
//		uint64_t both = key.first;
//
//		return (hash<uint64_t>()(both<<32 | key.first));
////		return (hash<_snpid>()(key.first) ^ (hash<_snpid>()(key.second) << 1));
//	}
//};
//typedef unordered_set<_pair, PairHasher> _pairs_set;
//typedef unordered_map<uint64_t,_pair> _pairs_set;


#endif /* TYPEDEFS_H_ */
