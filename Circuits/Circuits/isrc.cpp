#include "isrc.h"
Isrc::Isrc(string c, Node* a, Node* b, double x, double y) :Component(c, a, b) {
	name = c;
	n1 = a;
	n2 = b;
	phase = x;
	mag = y;
}
double Isrc::getmag() {
	return mag;
}
void Isrc::setnode1(Node* a) {
	n1 = a;
}
void Isrc::setnode2(Node* a) {
	n2 = a;
}
double Isrc::getphase() {
	return phase;
}
complex<double> Isrc::getcomplex() {
	x.real(mag * cos(phase));
	x.imag(mag * sin(phase));
}
Isrc::~Isrc() {}