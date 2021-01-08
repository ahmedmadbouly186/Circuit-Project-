#pragma once
#include<complex>
#include "Node.h"
#include "Component.h"
class Isrc : public Component
{
private:
	Node *n1;
	Node *n2;
	string name;
	double phase;
	double mag;
	complex<double> x;
public:
	Isrc(string, Node*, Node*, double, double);
	void setnode1(Node*);
	void setnode2(Node*);
	double getmag();
	double getphase();
	complex<double> getcomplex();
	~Isrc();
	void zerosetters() {
		phase = mag = 0;
		x.real(0);
		x.imag(0);

	}
};

