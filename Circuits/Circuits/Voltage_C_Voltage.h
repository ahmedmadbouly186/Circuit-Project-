#pragma once
#include "Dependent_sources.h"
class Voltage_C_Voltage: Dependent_sources
{
public:
	Voltage_C_Voltage(Node* x1, Node* x2, string n, Node* d1, Node* d2, double coff);
};

