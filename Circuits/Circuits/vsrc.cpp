#include "vsrc.h"
Vsrc::Vsrc(string c, Node* a, Node* b,double x,double y):Component(c,a,b) {
	name = c;
	n1 = a;
	n2 = b;
	phase = x;
	mag = y;
}
double Vsrc::getmag() {
	return mag;
}
void Vsrc::setnode1(Node* a) {
	n1 = a;
}
void Vsrc::setnode2(Node* a) {
	n2 = a;
}
double Vsrc::getphase() {
	return phase;
}
complex<double> Vsrc::getcomplex() {
	x.real(mag * cos(phase));
	x.imag(mag * sin(phase));
}
Vsrc::~Vsrc() {}
Node* Vsrc::getnode1() {
	return n1;
}
Node* Vsrc::getnode2() {
	return n2;
}