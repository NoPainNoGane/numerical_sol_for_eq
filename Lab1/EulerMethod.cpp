#include "EulerMethod.h"

EulerMethod::EulerMethod()
{

}

double EulerMethod::dvdt(double x)
{
	return (3 * x * x + 2 * x);
}

double EulerMethod::dxdt(double v)
{
	return -v;
}

void EulerMethod::print_data(string name, int N, vector<double>& x, vector<double>& v, vector<double>& t)
{
	string dir = "E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/";
	ofstream foutX(dir + name + "X" + to_string(N) + ".txt");
	cout.precision(10);
	foutX.precision(10);

	for (int i = 0; i < N + 1; i++)
		foutX << t[i] << " " << x[i] << endl;

	foutX.close();

	ofstream foutV(dir + name + "V" + to_string(N) + ".txt");
	cout.precision(10);
	foutV.precision(10);

	for (int i = 0; i < N + 1; i++)
		foutV << t[i] << " " << v[i] << endl;

	foutV.close();

	x.clear();
	v.clear();
	t.clear();
}

void EulerMethod::CalculationEuler(int N, double tau, vector<double>& x, vector<double>& v, vector<double>& t)
{
	x.resize(N+1);
	v.resize(N+1);
	t.resize(N+1);

	x[0] = x0;
	v[0] = v0;
	t[0] = t0;

	for (size_t i = 1; i < N + 1; i++) {
		t[i] = tau * i;
		v[i] = v[i - 1] + tau * dvdt(x[i - 1]);
		x[i] = x[i - 1] + tau * dxdt(v[i - 1]);
	}

	print_data("Euler-", N, x, v, t);
	PrintErrorEuler(N);
}

void EulerMethod::EulerForRunge(int N, double tau1, double tau2)
{

	std::vector<double> x1;
	std::vector<double> v1;
	std::vector<double> x2;
	std::vector<double> v2;
	std::vector<double> t1;

	std::vector<double> x;
	std::vector<double> v;


	x.resize(N + 1);
	v.resize(N + 1);
	x1.resize(N + 1);
	v1.resize(N + 1);
	t1.resize(N + 1);
	x2.resize(2 * N + 1);
	v2.resize(2 * N + 1);

	x1[0] = x2[0] = 0; v1[0] = v2[0] = -2 * sqrt(9.1) * 0.01; t1[0] = 0;

	double sigma = (tau1 / 2) / ((tau1 / 2) - tau1);

	for (size_t i = 1; i <= N; i++) {

		t1[i] = tau1 * i;
		v1[i] = v1[i - 1] + tau1 * dvdt(x1[i - 1]);
		x1[i] = x1[i - 1] + tau1 * dxdt(v1[i - 1]);
	}

	for (size_t i = 1; i <= 2 * N; i++) {

		v2[i] = v2[i - 1] + tau2 * dvdt(x2[i - 1]);
		x2[i] = x2[i - 1] + tau2 * dxdt(v2[i - 1]);
	}

	for (size_t i = 0; i <= N; i++)
	{
		x[i] = sigma * x1[i] + (1. - sigma) * x2[i * 2];
		v[i] = sigma * v1[i] + (1. - sigma) * v2[i * 2];

	}

	print_data("EulerRunge", N, x, v, t1);
	//PrintDataER(N, x, v, t1);

}

void EulerMethod::PrintDataER(int numberOfIteration, std::vector<double>& x, std::vector<double>& v, std::vector<double>& t)
{
	std::ofstream fout("EulerR-" + std::to_string(numberOfIteration) + ".txt");
	std::cout.precision(10);
	fout.precision(10);

	for (int i = 0; i <= numberOfIteration; i++) {

		fout << x[i] << " " << t[i] << " " << v[i] << std::endl;
	}
	fout.close();

	x.clear();
	v.clear();
	t.clear();
}

void EulerMethod::PrintErrorEuler(int numberOfIteration) {
	std::vector<double> x;
	std::vector<double> v;
	std::vector<double> t;

	std::vector<double> xM;
	std::vector<double> vM;
	std::vector<double> tM;

	x.resize(numberOfIteration + 1);
	v.resize(numberOfIteration + 1);
	t.resize(numberOfIteration + 1);

	xM.resize(numberOfIteration + 1);
	vM.resize(numberOfIteration + 1);
	tM.resize(numberOfIteration + 1);

	std::ifstream finX;
	std::ifstream finV;
	std::ifstream finMx;
	std::ifstream finMv;
	finX.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/Euler-X" + to_string(numberOfIteration) + ".txt");
	finV.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/Euler-V" + to_string(numberOfIteration) + ".txt");
	finMx.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/maple_" + to_string(numberOfIteration) + "x.txt");
	finMv.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/maple_" + to_string(numberOfIteration) + "v.txt");
	for (int i = 0; i <= numberOfIteration; i++) {

		finX >> t[i] >> x[i];
		finV >> t[i] >> v[i];
		finMx >> tM[i] >> xM[i];
		finMv >> tM[i] >> vM[i];
	}

	finX.close();
	finV.close();
	finMx.close();
	finMv.close();

	string dir = "E:/USATU/3 курс/6 семестр/ТРС/Lab1/ErrorCalculations/";
	std::ofstream foutX(dir + "EulerError_X-" + std::to_string(numberOfIteration) + ".txt");
	std::cout.precision(10);
	foutX.precision(10);

	double delta = 0;
	for (int i = 0; i <= numberOfIteration; i++) {
		delta = abs(x[i] - xM[i]);
		foutX << t[i] << " " << delta << std::endl;
	}
	foutX.close();

	std::ofstream foutV(dir + "EulerError_V-" + std::to_string(numberOfIteration) + ".txt");
	std::cout.precision(10);
	foutV.precision(10);

	delta = 0;
	for (int i = 0; i <= numberOfIteration; i++) {
		delta = abs(v[i] - vM[i]);
		foutV << t[i] << " " << delta << std::endl;
	}
	foutV.close();

	x.clear();
	v.clear();
	t.clear();

	xM.clear();
	vM.clear();
	tM.clear();

}