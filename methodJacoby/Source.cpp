#include<iostream>
#include<iomanip>
#include <ctime>
using namespace std;

double** matrixOfEigenvectors;
const double EPS = pow(10, -7);
const int NUM_OF_ITERATIONS_LIMIT = 1000;

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
void printMatrix(double* matrix, int dim)
{
	for (int i = 0; i < dim; ++i)
	{
		cout << setw(10) << matrix[i];
	}
	cout << endl;
}
double** generateMatrix(int n)
{
	srand(time(0));
	double** matrix = new double* [n];
	for (int i = 0; i < n; i++)
	{
		matrix[i] = new double[n];
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			
			matrix[i][j] = matrix[j][i] = rand()%(100+100+1)-100;
		}
	}
	return matrix;
}
double** inputMatrix(int n)
{
	double** matrix = new double* [n];
	for (int i = 0; i < n; i++)
	{
		matrix[i] = new double[n];
	}

	cout << "Input matrix:\n";
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{

			cin >> matrix[i][j];
		}
	}
	return matrix;
}
bool checkIfEnought(double** matrix, int dim) {
	double sum = 0;
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim, j != i; ++j)
		{
			sum += pow(matrix[i][j], 2);
		}
	}
	return sum < EPS;

}
int sign(double a) {
	if (a == 0)  return 0;
	if (a > 0)  return 1;
	else return -1;
}
double** rotateMatrix(double** matrix, int dim, int row, int col) {
	double angle, cos, sin;
	if (matrix[row][row] != matrix[col][col])
	{
		angle = 2 * matrix[row][col] / (matrix[row][row] - matrix[col][col]);
		cos = sqrt((1 + 1 / sqrt(1 + pow(angle, 2))) / 2);
		sin = sign(angle) * sqrt((1 - 1 / sqrt(1 + pow(angle, 2))) / 2);
	}
	else
	{
		cos = sin = 1 / sqrt(2);
	}
	
	for (int i = 0; i < dim; ++i)
	{
		double temp1 = matrix[row][i];
		double temp2 = matrix[col][i];
		matrix[row][i] = matrix[row][i] * cos + matrix[col][i] * (sin);
		matrix[col][i] = temp1* (-sin) + temp2 * (cos);
	}
	
	for (int i = 0; i < dim; ++i)
	{
		double temp1 = matrix[i][row];
		double temp2 = matrix[i][col];
		matrix[i][col] = matrix[i][row] * (-sin) + matrix[i][col] * (cos);
		matrix[i][row] = temp1 * cos + temp2 * (sin);
	}
	for (int i = 0; i < dim; ++i)
	{
		double temp1 = matrixOfEigenvectors[i][row];
		double temp2 = matrixOfEigenvectors[i][col];
		matrixOfEigenvectors[i][col] = matrixOfEigenvectors[i][row] * (-sin) + matrixOfEigenvectors[i][col] * (cos);
		matrixOfEigenvectors[i][row] = temp1 * cos + temp2 * (sin);
	}
	return matrix;
}
void getResults(double** a, double** matrix, int dim) {
	
	double* r = new double[dim];

	for (int i = 0; i < dim; i++)
	{
		
		for (int t = 0; t < dim; t++)
		{
			r[t] = 0;
		}
		
		cout << "Eigenvalue: " << matrix[i][i] << " Eigenvector: ( " ;
	
		for (int j = 0; j < dim; j++)
		{
			r[j] = 0;
			for (int k = 0; k < dim; k++)
			{
				
				if (j == k) {
					r[j] += (a[j][k] - matrix[i][i]) * matrixOfEigenvectors[k][i];
				}
				else {
					r[j] += a[j][k] * matrixOfEigenvectors[k][i];
				}
			}
			cout << matrixOfEigenvectors[j][i] << " ";
		}
		cout << ")" << endl;
		cout << "Vector r: ";
		for (int t = 0; t < dim; t++)
		{
			cout << r[t] << "\t";
		}	
		cout << endl;
	}	
}
double** nextIterationOfJacobyMethod(double** matrix, int dim) {

	int OptimalElementRow;
	int OptimalElementCol;
	double currentMaxElement = 0;
	for (int i = 0; i < dim; i++)
	{
		for (int j = i + 1; j < dim; j++) {
			if (abs(currentMaxElement) < abs(matrix[i][j])) {
				currentMaxElement = matrix[i][j];
				OptimalElementRow = i;
				OptimalElementCol = j;
			}
		}
	}
	matrix = rotateMatrix(matrix, dim, OptimalElementRow, OptimalElementCol);
	return matrix;
}
void JacobyMethod(double** a)
{
	int dim = 10;
	double** matrix = new double* [dim];
	matrixOfEigenvectors = new double* [dim];
	for (int i = 0; i < dim; i++)
	{
		matrixOfEigenvectors[i] = new double[dim];
		matrix[i] = new double[dim];
	}
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			matrixOfEigenvectors[i][j] = 0;
			if (i == j) matrixOfEigenvectors[i][j] = 1;
		}
	}




	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			matrix[i][j] = a[i][j];
		}
	}
	printMatrix(matrix, dim);
	int numOfIterations = 0;
	while (!checkIfEnought(matrix, dim))
	{
		matrix = nextIterationOfJacobyMethod(matrix, dim);
		numOfIterations++;
	}
	cout << numOfIterations << endl;
	getResults(a, matrix, dim);
};
double findMax(double* vector, int dim) {
	double temp = 0;
	for (int i = 0; i < dim; i++)
	{
		if (abs(vector[i]) > abs(temp)) {
			temp = vector[i];
		}
	}
	return temp;
}
bool StopCriterion(double** matrix, double* eigenvector, double eigenvalue, int dim)
{
	double sum = 0;
	double temp = 0;
	for (int i = 0; i < dim; i++)
	{
		temp = 0;
		for (int j = 0; j < dim; j++)
		{
			if (i == j)
			{
				temp += (matrix[i][j] - eigenvalue) * eigenvector[j];
			}
			else {
				temp += matrix[i][j] * eigenvector[j];
			}
		}
		sum += pow(temp, 2);
	}
	return sqrt(sum) < EPS;

}
double* NextIterationOfPowerMethod(double** matrix, double* y, int dim)
{
	double* result = new double[dim];
	for (int i = 0; i < dim; i++)
	{
		result[i] = 0;
		for (int j = 0; j < dim; j++)
		{
			result[i] += matrix[i][j] * y[j];
		}

	}
	return result;
}
void Normalize(double* vector, int dim) {
	double norm = 0;
	for (int i = 0; i < dim; i++)  {
		norm += pow(vector[i], 2);
	}
	norm = sqrt(norm);
	for (int i = 0; i < dim; i++) {
		vector[i] /= norm;
	}
	
}
void PowerIterationMethod(double** matrix)
{
	
	int dim = 10;
	double* yPrevious = new double[dim];
	double* yNext = new double[dim];
	double* eigenValues = new double [dim];
	double maxEigenvalue = 0;
	for (int i = 0; i < dim; i++)
	{
		yPrevious[i] = yNext[i] = 1;
	}
	cout << "Power Iteration Method" << endl;
	cout << "Initial Matrix" << endl;
	printMatrix(matrix, dim);
	cout << "Initial vector" << endl;
	printMatrix(yPrevious, dim);

	int numOfIterations = 0;
	

	while ((!StopCriterion(matrix, yNext, maxEigenvalue, dim)) && numOfIterations < NUM_OF_ITERATIONS_LIMIT)
	{
		for (int i = 0; i < dim; i++)
		{
			yPrevious[i] = yNext[i];
		
		}
		Normalize(yPrevious, dim);
		
		yNext = NextIterationOfPowerMethod(matrix, yPrevious, dim);
		for (int i = 0; i < dim; i++)
		{
			eigenValues[i] = yNext[i]/yPrevious[i];
		}
		maxEigenvalue = findMax(eigenValues, dim);
		numOfIterations++;

	} 
	if (numOfIterations == NUM_OF_ITERATIONS_LIMIT)
	{
		cout << "Power metod differs with max num of iterations = "<< NUM_OF_ITERATIONS_LIMIT << endl << "Values on the last iteration: " << endl;
	}
	cout << "EigenVector" << endl;
	printMatrix(yNext, dim);
	cout << "EigenValue" << endl;
	printMatrix(eigenValues, dim);
	cout << "Max EigenValue: " << maxEigenvalue << endl;
	
	

	double temp, sum = 0;
	cout << "Innacuracy vector: " ;
	for (int i = 0; i < dim; i++)
	{
		temp = 0;
		for (int j = 0; j < dim; j++)
		{ 
			if (i == j)
			{
				temp += (matrix[i][j] - maxEigenvalue) * yNext[j];
			}
			else {
				temp += matrix[i][j] * yNext[j];
			}
		}
		sum += pow(temp, 2);
		cout << temp << "\t";
	}
	cout <<endl<< "Norm of Innacuracy vector:" << sqrt(sum)<< endl;
	cout << "Num Of Iterations: " << numOfIterations << endl;
	/*
	1 1 3
	1 5 1
	3 1 1

	*/

}

void main()
{
	int dim = 10;
	double** a = generateMatrix(dim);
	JacobyMethod(a);
	PowerIterationMethod(a);
}










/*
1 -3 -1
-3 1 1
-1 1 5
	*/
