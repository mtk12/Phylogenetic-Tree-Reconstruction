#pragma once
#ifndef COMPOSITIONVECTORS_H_
#define COMPOSITIONVECTORS_H_
using namespace std;

map<string, double> compositionVector_getScoreVectors(map<int, map<string, int> > KMERS_MAP, int SeqLength, int K);
double compositionVector_getDistance(map<string, double> seq1, map<string, double> seq2);
double ** compositionVector_getDistanceMatrix(vector<map<string, double> >Vectors, int size);
double ** compositionVector(map<int, map<int, map<string, int> > > Map, int seqL[], int K);

#endif 

