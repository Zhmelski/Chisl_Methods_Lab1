#include "Main_Header.h"

vector<double> LDLT_Factorization(const vector<vector<double>>& A, const vector<double>& b)
{
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0.0));
    vector<double> D(n, 0.0), y(n), z(n), x(n);

    // LDL^T разложение
    for (int j = 0; j < n; ++j)
    {
        D[j] = A[j][j];
        for (int k = 0; k < j; ++k)
            D[j] -= L[j][k] * L[j][k] * D[k];

        for (int i = j + 1; i < n; ++i)
        {
            L[i][j] = A[i][j];
            for (int k = 0; k < j; ++k)
                L[i][j] -= L[i][k] * L[j][k] * D[k];
            L[i][j] /= D[j];
        }

        L[j][j] = 1.0;
    }

    // Решение Ly = b
    for (int i = 0; i < n; ++i)
    {
        y[i] = b[i];
        for (int k = 0; k < i; ++k)
            y[i] -= L[i][k] * y[k];
    }

    // Решение Dz = y
    for (int i = 0; i < n; ++i)
        z[i] = y[i] / D[i];

    // Решение Lᵗx = z
    for (int i = n - 1; i >= 0; --i)
    {
        x[i] = z[i];
        for (int k = i + 1; k < n; ++k)
            x[i] -= L[k][i] * x[k];
    }

    return x;
}

vector<double> multiply(const vector<vector<double>>& A, const vector<double>& x)
{
    int n = A.size();
    vector<double> res(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < x.size(); ++j)
            res[i] += A[i][j] * x[j];
    return res;
}