#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#define x0 0
#define t0 0
#define v0 (-2 * sqrt(9.1) * 0.01) 
using namespace std;

class EulerMethod
{
	void PrintDataER(int numberOfIteration, std::vector<double>& x, std::vector<double>& v, std::vector<double>& t);
	void PrintErrorEuler(int numberOfIteration);
public:
	EulerMethod();

	void CalculationEuler(int N, double tau, vector<double>& x, vector<double>& v, vector<double>& t);
	void EulerForRunge(int N, double tau1, double tau2);

	double dvdt(double x);
	double dxdt(double v);
	void print_data(string name, int N, vector<double>& x, vector<double>& v, vector<double>& t);
};

