#include "Dependent_sources.h"
Dependent_sources::Dependent_sources(Node* x1, Node* x2, string n, Node* d1, Node* d2, double coff):Component(x1,x2,n)
{
	node_d1 = d1;
	node_d2 = d2;
	cofficient = coff;
}
void Dependent_sources::set_Noded1(Node* x)
{
	node_d1 = x;
}
void Dependent_sources::set_Noded2(Node* x)
{
	node_d2 = x;
}
void Dependent_sources::set_cofficient(double x)
{
	cofficient = x;
}
Node* Dependent_sources::get_Noded1()
{
	return node_d1;
}
Node* Dependent_sources::get_Noded2()
{
	return node_d2;
}
double Dependent_sources::get_cofficient()
{
	return cofficient;
}
