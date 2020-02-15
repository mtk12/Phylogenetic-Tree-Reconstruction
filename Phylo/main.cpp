#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <map>
#include <math.h> 
#include <iomanip>
#include <algorithm>

#include "compositionVectors.h"
#include "neighborJoining.h"
#include "GraphClass.h"


using namespace std;

vector< vector<string> > getSequencesFastaFile(string filename) 
{

	ifstream FASTA_FILE;

	FASTA_FILE.open(filename.c_str());

	vector <string>  FASTA_SEQUENCE_SET;
	string  FASTA_SEQUENCE = "";
	string line;

	string names = "";

	while (getline(FASTA_FILE, line)) {

		if (line.empty()) continue;

		if (line[0] == '>') {//line starts with '>' , save the current in vector and set new sequence
			names.append(line);
			if (!FASTA_SEQUENCE.empty()) {
				FASTA_SEQUENCE_SET.reserve(1);
				FASTA_SEQUENCE_SET.push_back(FASTA_SEQUENCE);
				FASTA_SEQUENCE = "";
			}

			continue;
		}
		FASTA_SEQUENCE.append(line);
	}

	//save the last sequence
	if (!FASTA_SEQUENCE.empty()) {
		FASTA_SEQUENCE_SET.reserve(1);
		FASTA_SEQUENCE_SET.push_back(FASTA_SEQUENCE);
		FASTA_SEQUENCE = "";
	}
	FASTA_SEQUENCE_SET.reserve(1);
	FASTA_SEQUENCE_SET.push_back(names);

	//Close file
	FASTA_FILE.close();

	vector< vector<string> > sets;
	sets.push_back(FASTA_SEQUENCE_SET);

	return sets;
}

vector< vector<string> > getSequencesFromFile(string filename) 
{
	return getSequencesFastaFile(filename);
}

vector<string> getSequencesNames(string SequencesNames) 
{
	vector<string> NameSet;
	std::string delimiter = ">";
	size_t pos = 0;
	std::string token;
	while ((pos = SequencesNames.find(delimiter)) != std::string::npos) {
		token = SequencesNames.substr(0, pos);
		SequencesNames.erase(0, pos + delimiter.length());
		if (token.empty()) { continue; }
		NameSet.reserve(1);
		NameSet.push_back(token);
	}
	NameSet.reserve(1);
	NameSet.push_back(SequencesNames);//push back the last name
	return NameSet;
}

map<int, map<string, int> > getKMers(string Sequence) 
{
	map<int, map<string, int> > Maps;

	int k = 1;
	//int kmer;
	while (k < 11) { // up to 10-mer
		//kmer = k;
		map<string, int> kmer_map;
		for (int i = k; i < Sequence.size() - k + 1; i++) {
			if (kmer_map[Sequence.substr(i - k, k)]) {//non zero -- increase
				kmer_map[Sequence.substr(i - k, k)] += 1;
			}
			else {//zero -- add
				kmer_map[Sequence.substr(i - k, k)] = 1;
			}
		}
		Maps[k] = kmer_map;
		k++;
	}
	return Maps;
}

int* getSequencesLength(vector<string> Set) 
{
	int * seqL = new int[Set.size()];
	int j = 0;
	for (vector<string>::iterator i = Set.begin(); i != Set.end(); ++i) {
		seqL[j] = (*i).length();
		j++;
	}
	return seqL;
}

double** normalization(double** Matrix, int size) {

	//MIN VALUE IS ALWAYS ZERO BECAUSE OF DIAGONAL
	double MAX = 0;

	for (int i = 0; i < size; i++) {
		for (int j = i + 1; j < size; j++) {
			if (MAX < Matrix[i][j]) {
				MAX = Matrix[i][j];
			}
		}
	}

	for (int i = 0; i < size; i++) {
		for (int j = i + 1; j < size; j++) {
			Matrix[i][j] = Matrix[i][j] / MAX;
			Matrix[j][i] = Matrix[i][j];
		}
	}

	return Matrix;
}

void writeToFile(string filename, double ** Matrix, int SEQUENCES_SIZE)
{
	ofstream of(filename.c_str(),ios::out);
	for (int i = 0; i< SEQUENCES_SIZE; i++) {
		for (int j = 0; j< SEQUENCES_SIZE; j++) {
			of << Matrix[i][j]<<"\t";
		}
		of << endl;
	}
	of.close();
	cout << "Results on: " << filename << " file." << endl;

}

void calculatePhyloTrees(string filename, int option, int kmer, bool writeMatrices) 
{
	string f_out_name;

	string Netwick = "Tree.nwk";
	ofstream ov(Netwick.c_str());

	vector< vector< string > > Sets = getSequencesFromFile(filename);

	int iter = 1;
	for (vector< vector <string> >::iterator Set = Sets.begin(); Set != Sets.end(); ++Set) 
	{
		cout << "Sequence Set " << iter << " /" << Sets.size() << endl << endl;
		string SequencesNames = (*Set).back();//get the last element of the Sequence Set
		(*Set).pop_back();//It's the Sequences' names!
		vector<string> SequencesNamesV = getSequencesNames(SequencesNames);

		int SEQUENCES_SIZE = (*Set).size();

		/*
		* BEFORE READING ALL SEQUENCES WE SHOULD
		* PREPROCESS THEM FOR 'N' OCCURENCES.
		*/
		

		map<int, map<int, map<string, int> > > kmers_freq;

		/*
		* READ SEQUENCES AND PARSE THEM
		*/
		int j = 1;
		cout << "Reading sequences.." << endl;
		for (vector<string>::iterator i = (*Set).begin(); i != (*Set).end(); ++i) 
		{
			cout << "Sequence " << j << " / " << SEQUENCES_SIZE <<"..."<< endl;
			kmers_freq[j] = getKMers((*i));
			j++;
		}
		int* seqL = getSequencesLength((*Set));
		double ** Matrix =compositionVector(kmers_freq, seqL, kmer);
		f_out_name = "CompositionVectors.txt";
		
		/* If true write the DistanceMatrix on file */
		
		/* NORMALIZE MATRIX TO [0, 1] RANGE */
		Matrix = normalization(Matrix, SEQUENCES_SIZE);

		if (writeMatrices) {
			ofstream of(f_out_name.c_str());
			for (int i = 0; i< SEQUENCES_SIZE; i++) {
				for (int j = 0; j< SEQUENCES_SIZE; j++) {
					of << *((*Matrix + j) + i) << "\t";
				}
				of << endl;
			}
			of.close();
			cout << "Results on: " << f_out_name << " file." << endl;
		}
		/*GET THE NEIGHBOR JOINING FORM*/
		Graph tree = neighborJoining(Matrix, SEQUENCES_SIZE, SequencesNamesV);
		tree.print();
		ov << tree.getNetwickTreeFormat() << ";" << endl;
		iter++;
	}
	ov.close();
	cout << "Netwick Tree Format on: " << Netwick << endl;
}

int main() 
{
	string filename;
	//filename = "sequence.fasta";
	filename = "sequence (2).fasta";

	cout << "\tFinding Distance Matrix by Composition Vectors..." << endl;

	int k = 4;

	calculatePhyloTrees(filename, 1, k, true);

	cout << "\n\nTree Generated" << endl;
	system("pause");
	string completeUrl = "http://etetoolkit.org/treeview/";
	system(string("start " + completeUrl).c_str());
	return 0;
}