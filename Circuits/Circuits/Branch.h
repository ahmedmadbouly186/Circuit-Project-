#pragma once
#include <iostream>
#include <complex>       
#include "Node.h"
using namespace std;

#include"vsrc.h"
class Branch {
private:
	int num = 0;
	Node* node1; Node* node2;
	complex <double>z; complex <double> current;
	Vsrc volt;
public:
	Branch(Node*, Node*);
	void setnode1(Node* N1);
	void setNode2(Node* N2);
	void setz(complex<double> x);
	void setcurrent(complex<double> c);
	void setvolt(Vsrc c);
	Node* getnode1();
	Node* getnode2();
	complex <double>getz();
	Vsrc getvolt();
	complex<double> getcurrent();
	////////////// sabry

	void Calculatecurrent();

	///////////////
	int getnum() {
		return num;
	}
};

