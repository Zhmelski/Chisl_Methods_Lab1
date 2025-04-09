#include "Main_Header.h"

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	cout << setw(35) << "Метод Гаусса для варианта 14";
	cout << "\n*******************************************\n";

	vector<vector<double>> A14 = {
		{2.21, 3.65, 1.69, 6.99},
		{8.30, 2.62, 4.10, 1.90},
		{3.92, 8.45, 7.78, 2.46},
		{3.77, 7.21, 8.04, 2.28}
	};
	vector<double> b14 = { -8.35, -10.65, 12.21, 15.45 };

	vector<vector<double>> A14_orig = A14;
	vector<double> b14_orig = b14;

	vector<double> x14 = GaussMethod(A14, b14);

	cout << "Решение по методу Гаусса:\n";
	for (int i = 0; i < x14.size(); ++i)
	{
		cout << "x" << i + 1 << ": " << x14[i] << endl;
	}
	
	vector<double> residual14 = CalculateResidual(A14_orig, b14_orig, x14);

	cout << "Вектор невязки:\n";
	for (double val : residual14) cout << val << endl;

	cout << "Норма вектора невязки: " << NormVector(residual14) << endl;

	vector<double> b14_new(b14_orig.size());
	for (int i = 0; i < b14_new.size(); ++i)
		for (int j = 0; j < x14.size(); ++j)
			b14_new[i] += A14_orig[i][j] * x14[j];

	vector <double> x14_2 = GaussMethod(A14_orig, b14_new);
	if (x14_2.empty()) return 1;

	double error14 = 0.0;
	double xNorm14 = 0.0;
	for (int i = 0; i < x14.size(); ++i)
	{
		error14 = max(error14, abs(x14_2[i] - x14[i]));
		xNorm14 = max(xNorm14, abs(x14[i]));
	}
	double relativeError14 = error14 / xNorm14;

	cout << "Относительная погрешность: " << relativeError14;

	cout << "\n*******************************************\n";
	cout << setw(35) << "Метод Гаусса для варианта 21";
	cout << "\n*******************************************\n";

	double lambda1 = 1.0, lambda2 = 1e3, lambda3 = 1e6;

	vector<vector<double>> A21 = {
		{2 * lambda1 + 4 * lambda2, 2 * (lambda1 - lambda2), 2 * (lambda1 - lambda2)},
		{2 * (lambda1 - lambda2), 2 * lambda1 + lambda2 + 3 * lambda3, 2 * lambda1 + lambda2 - 3 * lambda3},
		{2 * (lambda1 - lambda2), 2 * lambda1 + lambda2 - 3 * lambda3, 2 * lambda1 + lambda2 + 3 * lambda3}
	};

	vector<double> b21 = {
		-4 * lambda1 - 2 * lambda2,
		-4 * lambda1 + lambda2 + 9 * lambda3,
		-4 * lambda1 + lambda2 - 9 * lambda3
	};

	vector<vector<double>> A21_orig = A21;
	vector<double> b21_orig = b21;

	vector<double> x21 = GaussMethod(A21, b21);

	cout << "Решение по методу Гаусса:\n";
	for (int i = 0; i < x21.size(); ++i)
	{
		cout << "x" << i + 1 << ": " << x21[i] << endl;
	}

	vector<double> residual21 = CalculateResidual(A21_orig, b21_orig, x21);

	cout << "Вектор невязки:\n";
	for (double val : residual21) cout << val << endl;

	cout << "Норма вектора невязки: " << NormVector(residual21) << endl;

	vector<double> b21_new(b21_orig.size());
	for (int i = 0; i < b21_new.size(); ++i)
		for (int j = 0; j < x21.size(); ++j)
			b21_new[i] += A21_orig[i][j] * x21[j];

	vector <double> x21_2 = GaussMethod(A21_orig, b21_new);
	if (x21_2.empty()) return 1;

	double error21 = 0.0;
	double xNorm21 = 0.0;
	for (int i = 0; i < x21.size(); ++i)
	{
		error21 = max(error21, abs(x21_2[i] - x21[i]));
		xNorm21 = max(xNorm21, abs(x21[i]));
	}
	double relativeError21 = error21 / xNorm21;

	cout << "Относительная погрешность: " << relativeError21;

	cout << "\n*******************************************\n";
	cout << "\n===========================================\n";
	cout << setw(46) << "Метод LDLT-Факторизации для варианта 21";
	cout << "\n*******************************************\n";

	vector<double> x21_LDLT = LDLT_Factorization(A21_orig, b21_orig);

	cout << "Решение по методу LDLT-Факторизации:\n";
	for (int i = 0; i < x21_LDLT.size(); ++i)
	{
		cout << "x" << i + 1 << ": " << x21_LDLT[i] << endl;
	}

	auto Ax = multiply(A21_orig, x21_LDLT);
	vector<double> F(x21_LDLT.size());
	for (int i = 0; i < x21_LDLT.size(); ++i)
		F[i] = Ax[i] - b21_orig[i];

	cout << "Норма невязки: " << NormVector(F) << endl;

	return 0;
}