#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>

using namespace std;

class RungeKut
{
public:
	RungeKut();
	double Dx_Dt(double v);
	double Dv_Dt(double x);
	void PrintData(string name, int N, vector<double>& x, vector<double>& v, vector<double>& t);
	void Calculation(double tau, unsigned int numberOfIteration, vector<double>& x, vector<double>& v,
		vector<double>& t);
	void PrintErrorRK(int numberOfIteration);
};

