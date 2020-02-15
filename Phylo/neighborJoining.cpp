#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <map>
#include <cmath> 
#include <iomanip>

#include "neighborJoining.h"
#include "GraphClass.h"
using namespace std;
/*
* 1.GET THE DISTANCE MATRIX AND CALCULATE MATRIX Q
* 2.
*/

double* calculateAllDistances(double** DistanceMatrix, int size) {


	//Initialize Matrix to Zeros
	double* Table = new double[size];
	for (int i = 0; i < size; ++i) {
		Table[i] = 0;
	}
	for (int i = 0; i<size; i++) {
		for (int j = 0; j<size; j++) {
			Table[i] += DistanceMatrix[i][j];
		}
	}
	return Table;
}

double** neighborJoining_getQMatrix(double ** DistanceMatrix, double *AllDistances, int size) {

	double** Matrix = new double*[size];
	for (int i = 0; i < size; ++i) {
		Matrix[i] = new double[size];
		for (int j = 0; j < size; j++) {
			Matrix[i][j] = ((double)size - 2.0)*DistanceMatrix[i][j] - AllDistances[i] - AllDistances[j];
		}
	}
	return Matrix;
}

Graph initializeTree(int size, std::vector<std::string> SequencesNames) {

	Graph graph;
	graph.addVertex("centraal");
	for (std::vector<std::string>::iterator i = SequencesNames.begin(); i != SequencesNames.end(); ++i) {
		graph.addVertex((*i));
		graph.addEdge("centraal", (*i), 1);
	}
	return graph;
}

std::string  synthesizeName(std::string name1, std::string name2, double distance1, double distance2) {

	std::string dist1 = static_cast<std::ostringstream*>(&(std::ostringstream() << distance1))->str();
	std::string dist2 = static_cast<std::ostringstream*>(&(std::ostringstream() << distance2))->str();

	std::string name = "(";

	//add first
	name.append(name1);
	name.append(":");
	name.append(dist1);
	name.append(",");

	//add second
	name.append(name2);
	name.append(":");
	name.append(dist2);
	name.append(")");
	return name;
}
int findPosInDistance(int i, int min1, int min2, int size) {

	int biggerMin = (min1>min2) ? min1 : min2;
	int smallerMin = (min1<min2) ? min1 : min2;

	if (i < smallerMin) return i; // if smaller than smaller min return it!
	if (i >= smallerMin) i = i + 1;
	if (i >= biggerMin) i = i + 1;//skip two rows/columns
	return i;
}

Graph neighborJoining(double ** DistanceMatrix, int size, vector<string> SequencesNames) {

	//Create the graph for Neighbor Joining algorithm.
	//First initialize graph as it's all nodes connected to a central node.
	Graph tree = initializeTree(size, SequencesNames);

	while (size != 1) {

		//Get Sum of Distances for every Row
		double* AllDistances = calculateAllDistances(DistanceMatrix, size);

		//Get Q Matrix
		double** QMatrix = neighborJoining_getQMatrix(DistanceMatrix, AllDistances, size);

		double minQ = 1000000.0;
		int min1 = 0, min2 = 0;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (i == j) continue;
				if (QMatrix[i][j] < minQ) {
					minQ = QMatrix[i][j];
					min1 = i;
					min2 = j;
				}
			}
		}

		/*
		Negative branch lengths
		As the neighbor-joining algorithm seeks to represent the data in the form of an additive tree,
		it can assign a negative length to the branch. Here the interpretation of branch lengths as an
		estimated number of substititions gets into difficulties. When this occurs it is adviced to set
		the branch length to zero and transfer the difference to the adjacent branch length so that the
		total distance between an adjacent pair of terminal nodes remains unaffected. This does not
		alter the overall topology of the tree (Kuhner and Felsenstein, 1994).
		*/
		double NodeDistance1 = 0.5*DistanceMatrix[min1][min2] + 0.5*(1.0 / (double)size - 2)*(AllDistances[min1] - AllDistances[min2]);
		double NodeDistance2 = DistanceMatrix[min1][min2] - NodeDistance1;
		if (NodeDistance1 < 0) {
			NodeDistance2 += std::abs(NodeDistance1);
			NodeDistance1 = 0;
		}
		else if (NodeDistance2 < 0) {
			NodeDistance1 += std::abs(NodeDistance2);
			NodeDistance2 = 0;
		}

		//add Nodes in the Graph and fix Connections
		std::string newnodename = synthesizeName(SequencesNames.at(min1), SequencesNames.at(min2), NodeDistance1, NodeDistance2);
		tree.addVertex(newnodename);
		tree.addEdge("centraal", newnodename, 1);
		tree.destroyEdge("centraal", SequencesNames.at(min1));
		tree.destroyEdge("centraal", SequencesNames.at(min2));
		tree.addEdge(newnodename, SequencesNames.at(min1), NodeDistance1);
		tree.addEdge(newnodename, SequencesNames.at(min2), NodeDistance2);

		//Update Elements for next Iteration
		double** distTemp = new double*[size - 1];
		int it, jt;
		for (int i = 0; i < size - 1; i++) {
			distTemp[i] = new double[size - 1];
			for (int j = 0; j < size - 1; j++) {
				distTemp[i][j] = 0;
			}
		}
		for (int i = 0; i < size - 1; i++) {
			if (i == size - 2) {
				continue;
			}
			it = findPosInDistance(i, min1, min2, size);
			for (int j = 0; j < size - 1; j++) {
				if (j == size - 2) {
					distTemp[i][j] = 0.5*(DistanceMatrix[it][min1] + DistanceMatrix[it][min2] - DistanceMatrix[min1][min2]);
					distTemp[j][i] = 0.5*(DistanceMatrix[it][min1] + DistanceMatrix[it][min2] - DistanceMatrix[min1][min2]);;
				}
				else {
					jt = findPosInDistance(j, min1, min2, size);
					distTemp[i][j] = DistanceMatrix[it][jt];
				}
			}
		}

		delete[] DistanceMatrix;
		DistanceMatrix = distTemp;

		//Remove the distant element first 
		if (min2 > min1) {
			SequencesNames.erase(SequencesNames.begin() + min2);
			SequencesNames.erase(SequencesNames.begin() + min1);
		}
		else {
			SequencesNames.erase(SequencesNames.begin() + min2);
			SequencesNames.erase(SequencesNames.begin() + min1);
		}
		SequencesNames.push_back(newnodename);
		size = size - 1;

	}

	return tree;
}