#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include "EulerMethod.h"
#include "AdamsMethod.h"
#include "RungeKut.h"

#include "KRS.h"


int main()
{
	setlocale(LC_ALL, "Russian");

	vector<double> x;
	vector<double> v;
	vector<double> t;
	int N;

	cout << "Введите число узлов: ";
	cin >> N;

	double a = 0, b = 10; // где а = t[0]
	double tau = double(b - a) / double(N);

	EulerMethod EM;
	EM.CalculationEuler(N, tau, x, v, t);
	EM.EulerForRunge(N, tau, tau / 2);

	AdamsMethod AM;
	AM.CalculationAdams(N, tau, x, v, t);

	RungeKut RK;
	RK.Calculation(tau, N, x, v, t);

	KRS KM;

}


