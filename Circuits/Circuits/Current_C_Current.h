#pragma once
#include "Dependent_sources.h"
class Current_C_Current :public Dependent_sources
{
private:
	Component* name_d;
public:
	Current_C_Current(Node* x1, Node* x2, string n, Node* d1, Node* d2, Component* n_d,double coff);
	void set_name_d(Component* x);
	Component* get_name_d();

};

