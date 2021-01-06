#pragma once
#include <iostream>
#include <complex>       
#include "Node.h"
using namespace std;

class Branch {
private:
	Node* node1; Node* node2;
	complex <double>z; complex <double> current;
public:
	Branch(Node*, Node*);
	void setnode1(Node* N1);
	void setNode2(Node* N2);
	void setz(complex<double> x);
	void setcurrent(complex<double> c);
	Node* getnode1();
	Node* getnode2();
	complex <double>getz();
	complex<double> getcurrent();

};

