#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h> 

#include "compositionVectors.h"

using namespace std;


//Define Iterators for each Level in the Map
typedef map<int, map<int, map<string, int> > >::iterator sequence_it;
typedef map<int, map<string, int> >::iterator kmers_it;
typedef map<string, int>::iterator freq_it;
typedef map<string, double>::iterator score_it;

/*
* Get all scores for each k-mer set.
*/
map<string, double> compositionVector_getScoreVectors(map<int, map<string, int> > KMERS_MAP, int SeqLength, int K) {

	map<string, double> KMER_SCORE;

	map<string, int> kmer_map = KMERS_MAP[K];

	//LEN: for finding probabilities
	int LEN = (SeqLength - K + 1);
	//LENS: using in prob0
	double LENS = (double)(LEN - 2) / ((double)(LEN - 1)*(double)(LEN - 1));

	for (freq_it i = kmer_map.begin(); i != kmer_map.end(); i++) {

		int frequency = i->second;
		string kmer = i->first;

		double prob = (double)frequency / LEN;

		double prob0_freq = (double)KMERS_MAP[K - 1][kmer.substr(0, K - 1)] * (double)KMERS_MAP[K - 1][kmer.substr(1, K - 1)] / (double)KMERS_MAP[K - 2][kmer.substr(1, K - 2)];
		double prob0 = prob0_freq*LENS;

		double A = (prob - prob0) / prob0;

		KMER_SCORE[kmer] = A;

	}
	return KMER_SCORE;
}

double compositionVector_getDistance(map<string, double> seq1, map<string, double> seq2) {

	double up = 0, downa = 0, downb = 0;
	for (score_it i = seq1.begin(); i != seq1.end(); i++) {

		string seq1_key = i->first;
		double seq1_value = i->second;

		if (seq2[seq1_key]) {
			up += seq1_value*seq2[seq1_key];
			downa += seq1_value*seq1_value;
			downb += seq2[seq1_key] * seq2[seq1_key];
		}
	}
	double C = up / sqrt(downa*downb);

	return (1 - C) / 2;
}

double ** compositionVector_getDistanceMatrix(vector<map<string, double> >Vectors, int size) {

	//Initialize Matrix to Zeros
	double** Matrix = new double*[size];
	for (int i = 0; i < size; ++i) {
		Matrix[i] = new double[size];
		for (int j = 0; j < size; j++) {
			Matrix[i][j] = 0;
		}
	}

	int i = 0, j = 0;

	for (vector<map<string, double> >::iterator seq1 = Vectors.begin(); seq1 != Vectors.end(); ++seq1) {
		map<string, double> seq1_map = (*seq1);
		j = i + 1;
		vector<map<string, double> >::iterator seq2 = Vectors.begin();
		for (advance(seq2, j); seq2 != Vectors.end(); ++seq2) {
			if (i == j) {
				Matrix[i][j] = 0;
				j++;
				continue;
			}
			map<string, double> seq2_map = (*seq2);
			Matrix[i][j] = Matrix[j][i] = compositionVector_getDistance(seq1_map, seq2_map);
			j++;
		}
		i++;
	}

	return Matrix;
}

/**
* Composition Vector Algorithm
* For each sequence get the k-mer you want.
* Compute the probability of appearence of each k-word
* Compute the probability of k-strings(MM)
* Compute score A for each k-word
* Get all scores for each k-word in a map
*/
double ** compositionVector(map<int, map<int, map<string, int> > > Map, int seqL[], int K) {

	int j = 0;

	vector< map<string, double> > maps;
	map<string, double> temp;
	for (sequence_it i = Map.begin(); i != Map.end(); i++) {
		temp.clear();

		//kmer map;
		map<int, map<string, int> > KMERS_MAP = i->second;
		temp = compositionVector_getScoreVectors(KMERS_MAP, seqL[j], K);
		j++;
		maps.push_back(temp);
	}

	//NOW I HAVE ALL SCORES OF A KMER FOR ALL SEQUENCES
	//DERIVE THE DISTANCE MATRIX FOR ALL SEQUENCES
	double **DistanceMatrix = compositionVector_getDistanceMatrix(maps, Map.size());
	return DistanceMatrix;
}