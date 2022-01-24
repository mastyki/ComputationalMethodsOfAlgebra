#include <iostream>
#include <iomanip>

using namespace std;

int  i, j, k;



void calculatePogr(double* x, int dim)
{
	double* y = new double[dim];
	double m, n;
	double pogr;
	for (int i = 0; i < dim; i++)
	{
		y[i] =( i + 1);
	}
	cout << "y" << y[0] << endl;
	cout << "x" << x[0] << endl;
	n = y[0];
	n	-=x[0];
	/*for (int i = 0; i < dim; i++)
	{
		if (abs(y[i] - x[i]) > n)
		{
			cout << "y" << y[i] << endl;
			cout << "x" << x[i] << endl;
			n = (y[i] - x[i]);

		}
		cout << n << endl;
	}*/
	
	cout << y[dim-1] << endl;
	pogr = x[dim-1]/ y[dim - 1] - y[dim- 1]/y[dim-1];
	cout << "pogr   " << pogr;

}
void matrixMultiplication(double** first, double** second, int dim)
{
	double** result = new double* [dim];
	for (int i = 0; i < dim; ++i)
	{
		result[i] = new double[dim];
	}



	cout << "Fir" << endl;
	for (int i = 0; i < dim; ++i)
	{

		for (int j = 0; j < dim ; ++j)
		{
			cout << setw(15) << first[i][j];
		}

		cout << endl;
	}

	cout << "SEc" << endl;
	for (int i = 0; i < dim; ++i)
	{

		for (int j = 0; j < dim ; ++j)
		{
			cout << setw(15) << second[i][j];
		}

		cout << endl;
	}

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			result[i][j] = 0;
			for (int k = 0; k < dim; k++)
			{
				result[i][j] += (first[i][k] * second[k][j]);
				cout << result[i][j] << endl;
			}
		
		}
	}

	
	cout << "Matrix hsxj" << endl;
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim ; ++j)
		{
			cout << setw(15) << result[i][j];
		}

		cout << endl;
	}
}
void findReverse(double** a, int dim)
{
	//прямой ход
	int j_max;
	double t;


	double** matrix = new double*[dim];
	for (int i = 0; i < dim; ++i)
	{
		matrix[i] = new double[dim * 2];
	}
	double** result = new double* [dim];
	for (int i = 0; i < dim; ++i)
	{
		result[i] = new double[dim];
	}
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			matrix[i][j] = a[i][j];
		}
	}
	for (int i = 0; i < dim; ++i)
	{
		for (int j = dim; j < dim*2; ++j)
		{
			if (i == j - dim)
			{
				matrix[i][j] = 1;
			}
			else
			{
				matrix[i][j] = 0;
			}
		}
	}
	cout << "Matrix" << endl;
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim * 2; ++j)
		{
			cout << setw(10) << matrix[i][j];
		}

		cout << endl;
	}

	for (int k = 0; k < dim; k++)
	{
		j_max = k;
		for (j = k; j < dim; j++)
			if (fabs(matrix[j_max][k]) < fabs(matrix[j][k]))
				j_max = j;

		//переставляем строки
		for (i = 0; i < dim * 2; i++)
		{
			t = matrix[k][i];
			matrix[k][i] = matrix[j_max][i];
			matrix[j_max][i] = t;
		}
		cout << "Matrix" << endl;
		for (int i = 0; i < dim; ++i)
		{

			for (int j = 0; j < dim * 2; ++j)
			{
				cout << setw(15) << matrix[i][j];
			}

			cout << endl;
		}
		//выполнение вычислений
	for (j = dim*2 -1; j >= k; j--)
		matrix[k][j] = matrix[k][j] / matrix[k][k];//делим на коэффиц. a[1][1],a[2][2] и т.д.
	for (i = k + 1; i < dim; i++)
		for (j = dim * 2 - 1; j >= k; j--)
			matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];//вычитаем из 2 ур-ия 1-ое умнож. на коэфф. и т.д. */

		
	}
	cout << "Matrix" << endl;
	for (int i = 0; i < dim; ++i)
	{

		for (int j = 0; j < dim * 2; ++j)
		{
			cout << setw(15) << matrix[i][j];
		}

		cout << endl;
	}


	//обратный ход
	for (int i = dim - 1; i >= 0; --i)
	{
		for (int j = 0; j < dim; ++j)
		{
			result[i][j] = matrix[i][dim + j ];// -1


			for (int k = dim - 1; k > i; --k)
			{
				result[i][j] -= matrix[i][k] * result[k][j];
			}
		}
	}
	


	//for (int i = 0; i < dim; i++) {
	//	for (int j = 0; j < dim; j++) {
	//		b[i][j] = 0;
	//		for (int k = 0; k < dim; k++) {
	//			//cout << a[i][k] << endl;
	//			b[i][j] = a[i][k] * result[k][j];
	//		}
	//	}
	//}
	//cout << "Matrix" << endl;
	//for (int i = 0; i < dim; ++i)
	//{
	//	for (int j = 0; j < dim ; ++j)
	//	{
	//		cout << setw(15) << b[i][j];
	//	}

	//	cout << endl;
	//}
	matrixMultiplication(result, a, dim);

}

void multiply(double**matrix, int dim)
{
	double vector[] = { 1,2,3,4,5,6,7,8,9,10 };
	double* result = new double[dim];
	for (int i = 0; i < dim; i++)
	{

		for (int j = 0; j < dim; j++)
		{
			result[i] = 0;
			for (int k = 0; k < dim; k++)// columns_1 = lines_2
			{
				result[i] += matrix[i][k] * vector[k];
			}
		}
	}

	for (int i = 0; i < dim; ++i)
	{
		matrix[i][dim] = result[i];
	}
}
void printMatrix(double** matrix, int dim)
{
	cout << "Matrix" << endl;
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim+1; ++j)
		{
			cout << setw(10) << matrix[i][j] ;
		}

		cout << endl;
	}
}
void randomMatrix(double** matrix, int dim)
{
	srand(time(0));
	cout << "Matrix" << endl;
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim+1; ++j)
		{
			matrix[i][j] = 1 + rand() % 10;
		}
	//	matrix[i][j + 1] = 1 + rand() % 10;
	}
}
int main()

{
	int dim;
	cout << "Poryadok:" << endl;
	cin >> dim;

	double* x = new double[dim];
	double** matrix = new double* [dim];

	for (int i = 0; i < dim; ++i)
	{
		matrix[i] = new double[dim + 1];
	}

	randomMatrix(matrix, dim);
	//findReverse(matrix, dim);

	multiply(matrix, dim);

	/*cout << "Koef" << endl;
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim + 1; ++j)
		{
			cout << "[" << i << "]" << "[" << j << "]";
			cin >> matrix[i][j];
		}
	}*/
	printMatrix(matrix, dim);
	//прямой ход
	int j_max;
	double t;
for (int k = 0; k < dim; k++)
{
	j_max = k;
	for (j = k; j < dim; j++)
		if (fabs(matrix[j_max][k]) < fabs(matrix[j][k]))
			j_max = j;

	//переставляем строки
	for (i = 0; i < dim + 1; i++)
	{
		t = matrix[k][i];
		matrix[k][i] = matrix[j_max][i];
		matrix[j_max][i] = t;
	}
	//если максимальный элемент нулевой
	int count;
	if (matrix[j_max][k] == 0)
	{
		count = 0;
		for (i = 0; i < dim; i++)
		{
			if (matrix[k][i] == 0)
				count++;
			printf("count of nulls: %d\n", count);
			if (count == dim && (matrix[k][dim + 1] == 0))//если коэфф. в строке = 0 и своб. член = 0
			{
				printf("too many solutions\n");
				return(0);
			}
			else if (count == dim && matrix[k][dim + 1] != 0)//если коэфф. в строке = 0 и своб. член не равен 0
			{
				printf("no solutions\n");
				return(0);
			}
		}
	}

	//выполнение вычислений
	for (j = dim; j >= k; j--)
		matrix[k][j] = matrix[k][j] / matrix[k][k];//делим на коэффиц a[1][1],a[2][2] и т.д.
	for (i = k + 1; i < dim; i++)
		for (j = dim; j >= k; j--)
			matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];//вычитаем из 2 ур-ия 1-ое умнож. на коэфф. и т.д. */
}

printMatrix(matrix, dim);


//обратный ход
for (int i = dim-1; i >= 0; --i)
{
	x[i] = matrix[i][dim];

	for (int j = dim - 1; j > i; --j)
	{
		x[i] -= matrix[i][j] * x[j];
	}
}

//вывод решения
cout << "roots" << endl;
for (i = 0; i < dim; i++)
	cout << setprecision(20) << x[i] << endl;
double max = 0;
for (i = 0; i < dim; i++)
{
	if (x[i] - i - 1 > max) {
		max = x[i] - i - 1;
	}
}

cout  << max;
//calculatePogr(x, dim);

return 0;
}

