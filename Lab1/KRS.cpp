#include "KRS.h"

KRS::KRS() {

	double T = 1.;
	N = 40;
	L = 2 * k0 * pow(N, 2) * T;
	h = 1. / N;
	tau = T / L;
	u_sol.resize(L + 1);

	for (size_t i = 0; i < L + 1; i++) {
		u_sol[i].resize(N + 1);
		for (size_t j = 0; j < N + 1; j++) {
			u_sol[i][j] = Solution(j * h, i * tau);
		}
	}

	std::cout << "N = " << N << ", L = " << L << ", tau = " << tau << ", h = " << h << std::endl;
	CalculationKRS(k0);
	PrintData(u, h, tau);
}

void KRS::CalculationKRS(double a)
{
	u.clear();
	u.resize(L + 1);
	for (size_t i = 0; i < L + 1; i++) {
		u[i].resize(N + 1);
	}

	for (size_t i = 0; i < N + 1; i++) {
		u[0][i] = phi(i * h);
	}

	for (size_t n = 0; n < L; n++) {
		for (size_t i = 1; i < N; i++)
		{
			u[n + 1][i] = a / pow(h, 2) * tau * (u[n][i - 1] - 2 * u[n][i] + u[n][i + 1]) + u[n][i] + function(i * h, tau * (n + 1));
		}
		u[n + 1][0] = (h * psi_0(tau * (n + 1)) - b0 * u[n + 1][1]) / (a0 * h - b0);
		u[n + 1][N] = (h * psi_1(tau * (n + 1)) + b1 * u[n + 1][N - 1]) / (a1 * h + b1);
	}
	
}

double KRS::function(double x, double t)
{
    return pow(x, 2) - 2 * t;
}

void KRS::PrintData(const std::vector<std::vector<double>>& u, double h, double tau)
{
	std::ofstream in("KRS_" + std::to_string(N) + "_" + std::to_string(L) + ".txt");
	for (int n = 0; n < u.size(); n++)
	{
		for (int i = 0; i < u[n].size(); i++)
		{
			in << i * h << "\t" << n * tau << "\t" << u[n][i] << "\n";
		}
	}
	in.close();
}
