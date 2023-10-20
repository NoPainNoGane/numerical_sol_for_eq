#include "AdamsMethod.h"

AdamsMethod::AdamsMethod()
{

}

void AdamsMethod::CalculationAdams(int N, double tau, vector<double>& x, vector<double>& v, vector<double>& t)
{
	x.resize(N+1);
	v.resize(N+1);
	t.resize(N+1);

	x[0] = 0;
	v[0] = -2 * sqrt(9.1) * 0.01;
	t[0] = 0;

	t[1] = t[0] + tau;
	v[1] = v[0] + tau * dvdt(x[0]);
	x[1] = x[0] + tau * dxdt(v[0]);

	for (size_t i = 2; i <= N; i++) {

		t[i] = i * tau;
		x[i] = x[i - 1] + 1.5 * tau * dxdt(v[i - 1]) - 0.5 * tau * dxdt(v[i - 2]);
		v[i] = v[i - 1] + 1.5 * tau * dvdt(x[i - 1]) - 0.5 * tau * dvdt(x[i - 2]);
	}

	print_data("Adams-", N, x, v, t);
	PrintErrorAdams(N);

}

double AdamsMethod::dvdt(double x)
{
	return (3 * x * x + 2 * x);
}

double AdamsMethod::dxdt(double v)
{
	return -v;
}

void AdamsMethod::print_data(string name, int N, vector<double>& x, vector<double>& v, vector<double>& t)
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

void AdamsMethod::PrintErrorAdams(int numberOfIteration)
{
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
	finX.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/Adams-X" + to_string(numberOfIteration) + ".txt");
	finV.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/Adams-V" + to_string(numberOfIteration) + ".txt");
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
	std::ofstream foutX(dir + "AdamsError_X-" + std::to_string(numberOfIteration) + ".txt");
	std::cout.precision(10);
	foutX.precision(10);

	double delta = 0;
	for (int i = 0; i <= numberOfIteration; i++) {
		delta = abs(x[i] - xM[i]);
		foutX << t[i] << " " << delta << std::endl;
	}
	foutX.close();

	std::ofstream foutV(dir + "AdamsError_V-" + std::to_string(numberOfIteration) + ".txt");
	std::cout.precision(10);
	foutV.precision(10);

	delta = 0;
	for (int i = 0; i <= numberOfIteration; i++) {
		delta = abs(v[i] - vM[i]);
		foutV << t[i] << " " << delta << std::endl;
	}
	foutV.close();


	/*std::vector<double> xE;
	std::vector<double> vE;
	std::vector<double> tE;

	xE.resize(numberOfIteration + 1);
	vE.resize(numberOfIteration + 1);
	tE.resize(numberOfIteration + 1);

	std::ifstream finXEuler;
	std::ifstream finVEuler;
	finXEuler.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/Euler-X" + to_string(numberOfIteration) + ".txt");
	finVEuler.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/Euler-V" + to_string(numberOfIteration) + ".txt");
	for (int i = 0; i <= numberOfIteration; i++) {

		finXEuler >> tE[i] >> xE[i];
		finVEuler >> tE[i] >> vE[i];

	}

	finXEuler.close();
	finVEuler.close();

	std::ofstream foutXEuler(dir + "Adams-Euler-Error_X-" + std::to_string(numberOfIteration) + ".txt");
	std::cout.precision(10);
	foutXEuler.precision(10);

	double delta = 0;
	for (int i = 0; i <= numberOfIteration; i++) {
		delta = abs(x[i] - xE[i]);
		foutXEuler << t[i] << " " << delta << std::endl;
	}
	foutXEuler.close();

	std::ofstream foutVEuler(dir + "Adams-Euler-Error_V-" + std::to_string(numberOfIteration) + ".txt");
	std::cout.precision(10);
	foutVEuler.precision(10);

	delta = 0;
	for (int i = 0; i <= numberOfIteration; i++) {
		delta = abs(v[i] - vE[i]);
		foutVEuler << t[i] << " " << delta << std::endl;
	}
	foutVEuler.close();*/

	x.clear();
	v.clear();
	t.clear();

	/*xE.clear();
	vE.clear();
	tE.clear();*/

	xM.clear();
	vM.clear();
	tM.clear();

}