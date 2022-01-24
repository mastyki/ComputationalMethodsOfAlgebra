#include<iostream>
#include<iomanip>
#include <ctime>
using namespace std;

void printMatrix(double* matrix, int n)
{
	//cout << "Matrix" << endl;
	cout << endl;
	for (int i = 0; i < n; ++i)
	{
		cout << setw(30) << setprecision(25) << matrix[i]<< endl;
	}
	cout << endl;
}

void inputMatrix(double* matrix, int n) {
	cout << "Введите матрицу" << endl;
	for (int i = 0; i < n; ++i)
	{
		cin >> matrix[i];
	}
}

void randomMatrix(double* matrix, int n)
{
	
	for (int i = 0; i < n; ++i)
	{
		matrix[i] = -100 + rand() % (200 + 1);
	}
}

// rand() % (b - a + 1) + a
void randomC(double* c, double* a, double* b, int n)
{
	srand(time(0));
	c[0] = rand() % int(abs(b[0]) + 6 * 2 - (abs(b[0]) + 6 + 1)) + abs(b[0]) + 6;
	for (int i = 1; i < n - 1; ++i)
	{
		c[i] = rand() % int(abs(b[i]) + abs(a[i]) + 6 * 2 - (abs(b[i])+ abs(a[i]) + 6 + 1)) + abs(b[i]) + abs(a[i]) + 6;
	}
	c[n-1] = rand() % int(abs(a[n-1]) + 6 * 2 - (abs(a[n-1]) + 6 + 1)) + abs(a[n-1]) + 6;
}
void defineF(double* f, double* c, double* a, double* b,double* exactSol, int n)
{
	f[0] = c[0] * exactSol[0] + b[0] * exactSol[1];
	f[n - 1] = a[n - 2] * exactSol[n - 2] + c[n - 1] * exactSol[n - 1];
	for (int i = 1,  j = 0; i < n - 1; ++i, ++j )
	{ 
		f[i] = a[j] * exactSol[i-1] + c[i] * exactSol[i] + b[i] * exactSol[i + 1];
	}

}
double findAccuracy(double* x, int n) {
	double max = 0;
	double accuracy;
	for (int i = 0; i < n; i++)
	{
		if (abs(x[i]-(i+1)) > max)
		{
			max = abs(x[i] - (i + 1));
		}
	}
	accuracy = max / n;
	//cout << setprecision(40)<< max << endl;
	return accuracy;
}
void main()
{
	srand(time(0));
	int n = 10;
	double* a = new double[n - 1];// ниже главной диагонали
	double* b = new double[n - 1];//выше главной диагонали
	double* c = new double[n];//на главной диагонали
	double* f = new double[n];//вектор правой части
	double* exactSol = new double[n];//точное решение
	double* x = new double[n];//столбец решений
	double gama;//вспомогательный коэффициент знаменателя
	double* beta = new double[n];//на главной диагонали
	double* alpha = new double[n-1];//на главной диагонали


	randomMatrix(a, n - 1);
	randomMatrix(b, n - 1);
	randomC(c, a, b, n);
	for (size_t i = 0; i < n; i++)
	{
		exactSol[i] = i + 1;
	}
	defineF(f, c, a,b, exactSol, n);

	cout << "a:\t";
	printMatrix(a, n - 1);
	cout << "b:\t";
	printMatrix(b, n - 1);
	cout << "c:\t";
	printMatrix(c, n);
	cout << "f:\t";
	printMatrix(f, n);
	cout << "exactSolution:\t";
	printMatrix(exactSol, n);

// прямая прогонка

	gama = c[0];
	alpha[0] = -b[0] / gama;
	beta[0] = f[0] / gama;

	for (size_t i = 1; i < n-1; i++)
	{
		gama = c[i] + a[i-1] * alpha[i - 1];
		alpha[i] = -b[i] / gama;
		beta[i] = (f[i] - a[i-1] * beta[i - 1]) / gama;
	}

	gama = c[n - 1] + a[n - 2] * alpha[n - 2];
	beta[n - 1] = (f[n - 1] - a[n - 2] * beta[n - 2]) / gama;

	// обратная прогонка

	x[n - 1] = beta[n - 1];
	for (int i = n-2; i >=  0; i--)
	{
		x[i] = alpha[i] * x[i + 1] + beta[i];
	}
	
	cout << "found Solution:\t";
	printMatrix(x, n);
	
	cout  << "Inaccuracy:" << setprecision(20) << findAccuracy(x, n) << endl;

}




/*a[0] = 5;
a[1] = 1;
b[0] = -1;
b[1] = 2;
c[0] = 2;
c[1] = 4;
c[2] = -3;
f[0] = 3;
f[1] = 6;
f[2] = 2;*/
/*
	cout << "a";
	for (size_t i = 0; i < n-1; i++)
{
	cin >> a[i];
}
cout << "b";
for (size_t i = 0; i < n - 1; i++)
{
	cin >> b[i];
}
cout << "c";
for (size_t i = 0; i < n ; i++)
{
	cin >> c[i];
}
cout << "f";
for (size_t i = 0; i < n; i++)
{
	cin >> f[i];
}*/
