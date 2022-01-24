#include<iostream>
#include<iomanip>
#include<math.h>
using namespace std;

double* matrixMultiplication(double** first, double* second, int dim);
double* generateF(double** matrix, int dim);
void printMatrix(double** matrix, int dim)
{
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			cout << setw(10) << matrix[i][j];
		}

		cout << endl;
	}
	cout << endl;
}
void inputMatrix(double** matrix, int dim)
{
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			cin >> matrix[i][j];
		}

	}
}
void inputMatrix(double* matrix, int dim)
{
	for (int i = 0; i < dim; ++i)
	{
		cin >> matrix[i];
	}
}

void printMatrix(double* matrix, int dim);
double** generateMatrix(int n)
{
	srand(time(0));
	double sum = 0; // required to generate doubleerval
	double** matrix = new double* [n];
	int a=-100, b =100;
	for (int i = 0; i < n; i++)
	{
		matrix[i] = new double[n];
	}
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			matrix[i][j] = 0;
		}
	}
	for (int j = 0; j < n; ++j)
	{
		for(int i = j + 1; i < n; ++i)
		{
			matrix[j][i] = matrix[i][j] =( a + rand() % (b - a + 1));
		}
	}
	for (int i = 0; i < n; i++)
	{
		sum = 0;
		for (int j = 0; j < n; j++)
		{
			sum += abs(matrix[i][j]);
		}
		a = sum + 6;
		b = sum + 60;
		matrix[i][i] = rand() % (b - a + 1) + a;

	}
	return matrix;
}

double* generateF(double** matrix, int dim)
{
	double* vector  = new double[dim];

	for (int i = 0; i < dim; i++)
	{
		vector[i] = i + 1;
	}
	cout << "Exact solution: " << endl;
	printMatrix(vector, dim);
	double* result = matrixMultiplication(matrix, vector, dim);

	return result;	
}
double* generateInitialAppr(double** matrix,double*b, int dim)
{
	
	double* result = new double[dim];
	for (int i = 0; i < dim; i++)
	{
		result[i] = b[i];
	}

	return result;
}
void matrixMultiplication(double** first, double** second, int dim)
{
	double** result = new double* [dim];
	for (int i = 0; i < dim; ++i)
	{
		result[i] = new double[dim];
	}

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			result[i][j] = 0;
			for (int k = 0; k < dim; k++)
			{
				result[i][j] += (first[i][k] * second[k][j]);
			}

		}
	}
}
double* matrixMultiplication(double** first, double* second, int dim)
{
	double* result = new double [dim];
	
	for (int i = 0; i < dim; i++)
	{
		result[i] = 0;
		for (int j = 0; j < dim; j++)
		{
			result[i] += (first[i][j] * second[j]);

		}
	}
	return result;
}
double scolarMultiplication(double* first, double*second, int dim)
{
	double result = 0;
	for (int i = 0; i < dim; i++)
	{
		result += first[i] * second[i];
	}
	return result;
}
double* calculateR(double** matrix, double* x, double* f, int dim)
{
	double* r = new double[dim];
	double* temp = matrixMultiplication(matrix, x, dim);
	for (int i = 0; i < dim; i++)
	{
		r[i] = temp[i] - f[i];
	}

	return r;
}
double* calculateX(double** matrix, double* x, double* r, int dim)
{
	double* result = new double[dim];
	
	for (int i = 0; i < dim; i++)
	{
		result[i] = x[i] - (scolarMultiplication(r, r, dim) / scolarMultiplication(matrixMultiplication(matrix, r, dim), r, dim) ) * r[i];
	}
	return result;

}
double* calculateX(double** matrix, double* x,double* f, double w, int dim)
{
	double* result = new double[dim];
	double temp = 0;
	
	for (int i = 0; i < dim; i++)
	{
		temp = 0;
		for (int j = 0; j <i ; j++)
		{
			temp += matrix[i][j] * result[j];
		}
		for (int j = i+1; j < dim; j++)
		{
			temp += matrix[i][j] * x[j];
		}
		result[i] = (1 - w) * x[i] + w * (f[i] - temp) / matrix[i][i];
		
	}
	return result;
}
double findMaxNorm(double** matrix, double* f, double*x, int dim)
{
	 double* temp = new double[dim];
	 temp = matrixMultiplication(matrix, x, dim);
	 for (int i = 0; i < dim; i++)
	 {
		 temp[i] -= f[i];
	 }
	 double res = 0;
	for (int i = 0; i < dim; i++)
	{
		if (abs(temp[i]) > res) 
		{
			res = abs(temp[i]);
		}
	}
	return res;
}
double findPogr(double* a, double* b, int dim)
{
	double res = 0;
	for (int i = 0; i < dim; i++)
	{
		if(abs(a[i]-b[i]) > res)
		{
			res = abs(a[i] - b[i]);
		}
	}
	return res;
}
void printMatrix(double* matrix, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		cout << matrix[i] << endl;
	}
	cout << endl;
}


void main()
{
	setlocale(LC_ALL, "Russian");
	
	int dim = 10;
	double* exactSolution = new double[dim];
	for (int i = 0; i < dim; i++)
	{
		exactSolution[i] = i + 1;
	}
	double** matrix = generateMatrix(dim);
	cout << "matrix" << endl;
	printMatrix(matrix, dim);
    double* f = generateF(matrix, dim);
	cout << "f" << endl;
	printMatrix(f, dim);

	double* x0 = generateInitialAppr(matrix, f, dim);
	cout << "Initial approximation: " << endl;
	printMatrix(x0, dim);
	
	double* r;
	r = calculateR(matrix, x0, f, dim);
	cout << "r" << endl;
	printMatrix(r, dim);
	double* x;
	x = calculateX(matrix, x0, r, dim);
	
	double eps = 1e-7;
	cout << "Accuracy: " << eps << endl;
	int max = 1;
	while(findMaxNorm(matrix, f, x, dim) > eps && max < 5000)
	{

		x = calculateX(matrix, x, r, dim);;
		r = calculateR(matrix, x, f, dim);
		
		max++;

	}
	cout << "Solution" << endl;
	printMatrix(x, dim);
	cout << "Num of iterations: " << max << endl;
	cout << "Neviazka: " << findMaxNorm(matrix, f, x, dim) << endl;
	cout << "Inaccuracy: " << findPogr(x, exactSolution, dim) << endl;


	double w[] = {0.2,0.5, 0.8, 1, 1.3, 1.5, 1.8};
	for (int i = 0; i < 7; i++)
	{
		cout << w[i] << endl;
		x = calculateX(matrix, x0, f, w[i], dim);
		

		int max = 1;
		while (findMaxNorm(matrix, f, x, dim) > eps && max < 5000)
		{

			x = calculateX(matrix, x, f, w[i], dim);;
			max++;

		}
		cout << "Solution" << endl;
		printMatrix(x, dim);
		cout << "Num of iterations: " << max << endl;
		cout << "Neviazka: " << findMaxNorm(matrix, f, x, dim) << endl;
		cout << "Inaccuracy: " << findPogr(x, exactSolution, dim) << endl;
	}
	

}

