#include "Node.h"
#include <iostream>
Node::Node(int n) {
	rank = n;
	 num_connections = 0;
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
	//for (int i = 0; i < 10; i++)
	//{
	//	if (connections[i] == Node2)
	//	{
	//		find = true;
	//	}
	//}
	if (find == false) {
		connections[num_connections++] = Node2;
	}
	if (num_connections > 2)
	{
		simple = false;
	}
}
void Node::setvoltage(complex <double> x)
{
	voltage = x;
}
complex <double> Node::getvoltage()
{
	return voltage;
}

int Node::getrank()
{
	return rank;
}
Node** Node::getNodeconnections()
{
	return connections;
}