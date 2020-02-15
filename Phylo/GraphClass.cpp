#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <map>
#include <math.h> 
#include <time.h>
#include <iomanip>

#include "GraphClass.h"
using namespace std;

void Graph::addVertex(string name) 
{
	int id;
	if (!vertices.empty()) { id = vertices.size(); }
	else { id = 0; }
	Vertex newVertex(name, id);
	vertices.push_back(newVertex);
}

Vertex Graph::findVertex(string name)
{

	for (vector<Vertex>::iterator i = vertices.begin(); i != vertices.end(); ++i) {
		if (name.compare((*i).getName()) == 0) {
			return (*i);
		}
	}
}

void Graph::addEdge(string name1, string name2, double dist) 
{

	Vertex node1 = findVertex(name1);
	Vertex node2 = findVertex(name2);

	Edge edge(node1.getId(), node2.getId(), dist);
	edges.push_back(edge);
}

bool Graph::destroyEdge(string name1, string name2) 
{

	Vertex node1 = findVertex(name1);
	Vertex node2 = findVertex(name2);

	for (vector<Edge>::iterator i = edges.begin(); i != edges.end(); ++i) {
		if (((*i).getOrigin() == node1.getId() && (*i).getDestination() == node2.getId()) ||
			((*i).getOrigin() == node2.getId() && (*i).getDestination() == node1.getId())) {
			//delete the edge
			edges.erase(i);
			return true;
		}
	}
	return false;
}

void Graph::print() 
{
	ofstream of("Graph.txt");
	for (vector<Vertex>::iterator i = vertices.begin(); i != vertices.end(); ++i) 
	{
		of << (*i).getName() << std::endl;
	}
	for (vector<Edge>::iterator i = edges.begin(); i != edges.end(); ++i)
	{
		of << (*i).getOrigin() << ", " << (*i).getDestination() << " => " << (*i).getDistance() << std::endl;
	}
	of.close();
}

string Graph::getNetwickTreeFormat() 
{
	return vertices.back().getName();
}