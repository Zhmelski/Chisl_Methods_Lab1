#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <Windows.h>
using namespace std;

vector<double> GaussMethod(vector<vector<double>> A, vector<double> b);
vector<double> CalculateResidual(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& x);
double NormVector(const vector<double>& v);

vector<double> LDLT_Factorization(const vector<vector<double>>& A, const vector<double>& b);
vector<double> multiply(const vector<vector<double>>& A, const vector<double>& x);