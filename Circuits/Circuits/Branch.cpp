#include "Branch.h"
#include <iostream>
#include <complex>       
#include "Node.h"
Branch::Branch(Node* N1, Node* N2)
{
	node1 = N1;
	node2 = N2;
}
void Branch:: setnode1(Node* N1)
{
	node1 = N1;
}
void Branch::setNode2(Node* N2)
{
	node2 = N2;
}
void Branch::setz(complex<double> x)
{
	z = x;
}
void Branch::setcurrent(complex<double> c)
{
	current = c;
}
Node* Branch::getnode1()
{
	return node1;
}
Node* Branch::getnode2()
{
	return node2;
}
complex <double>Branch::getz()
{
	return z;
}
complex<double> Branch::getcurrent()
{
	return current;
}