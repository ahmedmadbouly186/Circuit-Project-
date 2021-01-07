#include "Component.h"
Component::Component(string a, Node* n1, Node* n2)
{
	a = name;
	node1 = n1;
	node2 = n2;
}
Node* Component::get_node1()
{
	return node1;
}
Node* Component::get_node2()
{
	return node2;
}
string Component::get_name()
{
	return name;
}