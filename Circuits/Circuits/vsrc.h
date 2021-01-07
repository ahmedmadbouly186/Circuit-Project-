#pragma once
#include"Node.h"
#include "Component.h"
class Vsrc :public Component
{
private:
	Node* n1;
	Node* n2;
	string name;
	double phase;
	double mag;
	complex<double> x;

public:
	Vsrc(string, Node*, Node*,double,double);
	void setnode1(Node*);
	void setnode2(Node*);
	void zerosetters() {
		phase = mag = 0;
		x.real(0);
		x.imag(0);
		
	}
	Node* getnode1();
	Node* getnode2();
	string getname() {
		return name;
	}
	double getmag();
	double getphase();
	complex<double> getcomplex();
	~Vsrc();
};

