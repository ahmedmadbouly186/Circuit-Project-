#pragma once
#include"Component.h"

class Resistance :public Component
{
private:
	complex<double> z;
public:
	Resistance(double,string, Node*, Node*);
	void set_Impedance(double);
	complex <double> get_Impedance();
	~Resistance()
	{

	}
};

