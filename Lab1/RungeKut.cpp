#include "RungeKut.h"


RungeKut::RungeKut()
{
}

double RungeKut::Dx_Dt(double v)
{
	return -v;
}

double RungeKut::Dv_Dt(double x)
{
	return (3 * x * x + 2 * x);
}

void RungeKut::PrintData(string name, int N, vector<double>& x, vector<double>& v, vector<double>& t)
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

void RungeKut::Calculation(double tau, unsigned int numberOfIteration, vector<double>& x, vector<double>& v, vector<double>& t)
{
	x.resize(numberOfIteration + 1);
	v.resize(numberOfIteration + 1);
	t.resize(numberOfIteration + 1);

	x[0] = 0.;
	v[0] = -2 * sqrt(9.1) * 0.01;
	t[0] = 0.;

	double k1_x, k2_x, k3_x, k4_x;
	double k1_v, k2_v, k3_v, k4_v;

	for (size_t i = 0; i < numberOfIteration; ++i)
	{
		k1_x = Dx_Dt(v[i]);
		k1_v = Dv_Dt(x[i]);

		k2_x = Dx_Dt(v[i] + tau * k1_v / 2.);
		k2_v = Dv_Dt(x[i] + tau * k1_x / 2.);

		k3_x = Dx_Dt(v[i] + tau * k2_v / 2.);
		k3_v = Dv_Dt(x[i] + tau * k2_x / 2.);

		k4_x = Dx_Dt(v[i] + tau * k3_v);
		k4_v = Dv_Dt(x[i] + tau * k3_x);

		x[i + 1] = x[i] + (k1_x + 2. * k2_x + 2. * k3_x + k4_x) * (tau / 6.);
		v[i + 1] = v[i] + (k1_v + 2. * k2_v + 2. * k3_v + k4_v) * (tau / 6.);
		t[i] = tau * i;

	}

	PrintData("RungeKut-", numberOfIteration, x, v, t);
	PrintErrorRK(numberOfIteration);
}

void RungeKut::PrintErrorRK(int numberOfIteration)
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
	finX.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/RungeKut-X" + to_string(numberOfIteration) + ".txt");
	finV.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/RungeKut-V" + to_string(numberOfIteration) + ".txt");
	finMx.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/maple_" + to_string(numberOfIteration) + "x.txt");
	finMv.open("E:/USATU/3 курс/6 семестр/ТРС/Lab1/MethodCalculations/maple_" + to_string(numberOfIteration) + "x.txt");
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
	std::ofstream foutX(dir + "RungeCutError_X-" + std::to_string(numberOfIteration) + ".txt");
	std::cout.precision(10);
	foutX.precision(10);

	double delta = 0;
	for (int i = 0; i <= numberOfIteration; i++) {
		delta = abs(x[i] - xM[i]);
		foutX << t[i] << " " << delta << std::endl;
	}
	foutX.close();

	std::ofstream foutV(dir + "RungeCutError_V-" + std::to_string(numberOfIteration) + ".txt");
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