#pragma once
#include "Component.h"
class Dependent_sources: public Component
{
private:
	Node* node_d1;
	Node* node_d2;
	double cofficient;
public:
	Dependent_sources(Node* x1, Node* x2, string n, Node* d1, Node* d2, double coff);
	void set_Noded1(Node* x);
	void set_Noded2(Node* x);
	void set_cofficient(double x);
	Node* get_Noded1();
	Node* get_Noded2();
	double get_cofficient();
};