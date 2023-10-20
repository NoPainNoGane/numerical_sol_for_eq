#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#define a0 1.
#define b0 1.
#define a1 0.
#define b1 1.
#define k0 1.

class KRS
{
public:
	KRS();

private:
	int N, L;
	double a, b;
	double h;
	double tau;
	std::vector<std::vector<double>> u;
	std::vector<std::vector<double>> u_sol;

	void CalculationKRS(double a);
	double function(double x, double t);
	double Solution(double x, double t) {
		return t * pow(x, 2) + x;
	};
	double phi(double x) { return x; };
	double psi_0(double t) { return 1; };
	double psi_1(double t) { return 2 * t + 1; };

	void PrintData(const std::vector<std::vector<double>>& u, double h, double tau);
};

