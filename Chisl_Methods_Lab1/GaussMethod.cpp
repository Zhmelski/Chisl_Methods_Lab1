#include "Main_Header.h"

vector<double> GaussMethod(vector<vector<double>> A, vector<double> b)
{
	int n = A.size();
	vector<int> IOR(n);
	for (int i = 0; i < n; ++i) IOR[i] = i;

	//Прямой ход
	for (int k = 0; k < n; ++k)
	{
		//Выбор главного элемента (максимального по модулю)
		int p = k;
		double maxVal = 0.0;
		for (int i = k; i < n; ++i)
		{
			int row = IOR[i];
			if (abs(A[row][k] > maxVal))
			{
				maxVal = abs(A[row][k]);
				p = i;
			}
		}
		if (maxVal < 1e-12)
		{
			cerr << "Матрица вырождена!" << endl;
			exit(1);
		}
		if (p != k) swap(IOR[k], IOR[p]);

		int pivotRow = IOR[k];
		double pivot = A[pivotRow][k];
		for (int j = k; j < n; ++j) A[pivotRow][j] /= pivot; // Нормировка строки
		b[pivotRow] /= pivot;

		//Исключение переменных
		for (int i = k + 1; i < n; ++i)
		{
			int row = IOR[i];
			double factor = A[row][k];
			for (int j = k; j < n; ++j)
				A[row][j] -= factor * A[pivotRow][j];
			b[row] -= factor * b[pivotRow];
		}
	}

	//Обратный ход
	vector<double> x(n);
	for (int k = n - 1; k >= 0; --k)
	{
		int row = IOR[k];
		x[k] = b[row];
		for (int j = k + 1; j < n; ++j)
			x[k] -= A[row][j] * x[j];
	}

	return x;
}

//Вычисления вектора невязки (умножение матрицы на вектор)
vector<double> CalculateResidual(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& x)
{
	int n = A.size();
	vector<double> residual(n, 0.0);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			residual[i] += A[i][j] * x[j];
		}
		residual[i] -= b[i];
	}

	return residual;
}

//Вычисление нормы вектора нувязки
double NormVector(const vector<double>& v)
{
	double norm = 0.0;
	for (double val : v)
		norm = max(norm, abs(val));
	return norm;
}