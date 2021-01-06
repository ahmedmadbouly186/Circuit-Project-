#include "Component.h"
Component::Component(string a, int n1, int n2)
{
	a = name;
	node1 = n1;
	node2 = n2;
}
int Component::get_node1()
{
	return node1;
}
int Component::get_node2()
{
	return node2;
}
string Component::get_name()
{
	return name;
}