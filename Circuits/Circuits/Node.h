#pragma once
#include <iostream>
#include <complex>       
using namespace std;
class Node {
private:
	bool simple;
	int num_connections;
	double voltage;
	Node* connections[10];
public:
	Node();
	//return true if simple &&  false if non_simple
	bool ifsimple();
	void increment_connections(Node* Node2);
	void setvoltage(double x);
	double getvoltage();
	
};