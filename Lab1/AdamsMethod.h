#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

class AdamsMethod
{
public:
	AdamsMethod();
	void CalculationAdams(int N, double tau, vector<double>& x, vector<double>& v, vector<double>& t);
	double dvdt(double x);
	double dxdt(double v);
	void print_data(string name, int N, vector<double>& x, vector<double>& v, vector<double>& t);
	void PrintErrorAdams(int numberOfIteration);
};

