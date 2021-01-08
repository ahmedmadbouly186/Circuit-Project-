#include "Current_C_Current.h"
Current_C_Current::Current_C_Current(Node* x1, Node* x2, string n, Node* d1, Node* d2, Component* n_d, double coff):Dependent_sources(x1, x2, n, d1, d2, coff)
{
	name_d = n_d;
}
void Current_C_Current::set_name_d(Component* x)
{
	name_d = x;
}
Component* Current_C_Current::get_name_d()
{
	return name_d;
}
