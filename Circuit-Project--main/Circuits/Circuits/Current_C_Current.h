#pragma once
#include "Dependent_sources.h"
class Current_C_Current :public Dependent_sources
{
private:
	Component* name_d;
public:
	Current_C_Current(int x1, int x2, string n, int d1, int d2, Component* n_d,double coff);
	void set_name_d(Component* x);
	Component* get_name_d();

};

