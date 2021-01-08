#pragma once
#include <iostream>
#include <complex>       
#include "Node.h"
#include "isrc.h"
using namespace std;

#include"vsrc.h"
class Branch {
private:
	int num = 0;
	Node* node1; Node* node2;
	complex <double>z;
	complex <double> current;
	Vsrc *volt;
	Isrc* csource;
public:
	Branch(Node*, Node*);
	void setnode1(Node* N1);
	void setNode2(Node* N2);
	void setz(complex<double> x);
	void setcurrent(complex<double> c);
	void setvolt(Vsrc *c);
	void setcurentsource(Isrc*);
	Node* getnode1();
	Node* getnode2();
	complex <double>getz();
	Vsrc *getvolt();
	Isrc* getcsource();
	complex<double> getcurrent();
	////////////// sabry

	void Calculatecurrent();

	///////////////
	int getnum() {
		return num;
	}
};

