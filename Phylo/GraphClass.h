#pragma once
#ifndef GRAPHCLASS_H_
#define GRAPHCLASS_H_
using namespace std;

class Edge
{
private:
	int origin;
	int destination;
	double distance;

public:
	Edge(int org, int dest, double dist)
	{
		origin = org;
		destination = dest;
		distance = dist;
	}
	int getOrigin() { return origin; }
	int getDestination() { return destination; }
	double getDistance() { return distance; }

};

class Vertex
{
private:
	string name;
	int id;
public:
	Vertex(string n, int i)
	{
		name = n;
		id = i;
	}
	string getName() { return name; }
	int getId() { return id; }
};


class Graph
{
private:
	vector<Vertex> vertices;
	vector<Edge> edges;
public:
	Graph() {};
	Vertex findVertex(string name);
	void addVertex(string name);
	void addEdge(string name1, string name2, double distance);
	bool destroyEdge(string name1, string name2);

	void print();
	string getNetwickTreeFormat();
};

#endif 

