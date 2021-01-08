#include "Indepentsourcr.h"
Indepentsourcr::Indepentsourcr(string c, int a, int b, double x, double y) :Component(c, a, b) {

	phase = x;
	mag = y;
}
double Indepentsourcr::getmag() {
	return mag;
}

double Indepentsourcr::getphase() {
	return phase;
}
complex<double> Indepentsourcr::getcomplex() {
	x.real(mag * cos(phase));
	x.imag(mag * sin(phase));
	return x;
}
Indepentsourcr::~Indepentsourcr() {}
