#pragma once
#include "Component.h"
class Indepentsourcr : public Component
{
private:

	double phase;
	double mag;
	complex<double> x;

public:
	Indepentsourcr(string, int, int, double, double);

	void zerosetters() {
		phase = mag = 0;
		x.real(0);
		x.imag(0);

	}
	double getmag();
	double getphase();
	complex<double> getcomplex();
	~Indepentsourcr();
};

