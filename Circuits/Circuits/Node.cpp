#include "Node.h"
#include <iostream>
Node::Node() {
	int num_connections;
	bool simple = true;
	for (int i = 0; i < 10; i++)
	{
		connections[i] = NULL;
	}
}
bool Node::ifsimple() 
{
	return simple;
}
void Node::increment_connections(Node* Node2)
{
	bool find = false;
	for (int i = 0; i < 10; i++)
	{
		if (connections[i] == Node2)
		{
			find = true;
		}
	}
	if (find == false) {
		connections[num_connections++] = Node2;
	}
	if (num_connections > 2)
	{
		simple = false;
	}
}
void Node::setvoltage(double x)
{
	voltage = x;
}
double Node::getvoltage()
{
	return voltage;
}