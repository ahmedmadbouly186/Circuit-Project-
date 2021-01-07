#pragma once
#include<iostream>
#include<string>
#include<complex>
#include"Node.h"
using namespace std;
using namespace std::complex_literals;


class Component
{
private:
	Node* node1;
	Node* node2;
	string name;
public:
	Component(string , Node*, Node*);
	Node* get_node1();
	Node* get_node2();
	string get_name();
	virtual ~Component()
	{

	}
};

