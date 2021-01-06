#pragma once
#include<iostream>
#include<string>
#include<complex>
using namespace std;
using namespace std::complex_literals;


class Component
{
private:
	int node1;
	int node2;
	string name;
public:
	Component(string ,int ,int);
	int get_node1();
	int get_node2();
	string get_name();
	~Component()
	{

	}
};

