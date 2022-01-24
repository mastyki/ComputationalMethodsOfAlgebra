#include <iostream>
#include <exception>
#include "Matrix.h"

using namespace std;

matrix::matrix()
{
	m = 0;
	n = 0;
	p = 0;
}
matrix::matrix(int m, int n)
{
	this->m = m;
	this->n = n;
	p = new double* [m];
	for (int i = 0; i < m; i++)
		p[i] = new double[n];
}
matrix::matrix(const matrix& m)
{
	this->m = m.m;
	this->n = m.n;
	p = new double* [this->m];
	for (int i = 0; i < this->m; i++)
	{
		p[i] = new  double[n];
	}
	for (int i = 0; i < this->m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			p[i][j] = m.p[i][j];
		}
	}
}
matrix::~matrix()
{
	for (int i = 0; i < m; i++)
	{
		delete[] p[i];
	}
	delete[] p;
}
int matrix::getRow()
{
	return m;
}
int matrix::getCol()
{
	return n;
}
double& matrix::operator ()(int i, int j)
{
	return p[i][j];
}
matrix operator +(matrix& m1, matrix& m2)
{
	if (m1.getCol() != m2.getCol() || m1.getRow() != m2.getRow())
	{
		throw exception("Error!!!   The number of rows and columns does not match.");
	}
	matrix m3(m1);
	for (int i = 0; i < m3.getRow(); i++)
	{
		for (int j = 0; j < m3.getCol(); j++)
		{
			m3(i, j) += m2(i, j);
		}
	}
	return m3;
}
matrix operator -(matrix& m1, matrix& m2)
{
	if (m1.getCol() != m2.getCol() || m1.getRow() != m2.getRow())
	{
		throw exception("Error!!!   The number of rows and columns does not match.");
	}
	matrix m3(m1);
	for (int i = 0; i < m3.getRow(); i++)
	{
		for (int j = 0; j < m3.getCol(); j++)
		{
			m3(i, j) -= m2(i, j);
		}
	}
	return m3;
}
matrix& matrix::operator = (const matrix& m)
{
	if (this == &m)
	{
		return *this;
	}
	this->~matrix();
	this->m = m.m;
	this->n = m.n;
	p = new double* [this->m];
	for (int i = 0; i < this->m; i++)
	{
		p[i] = new  double[n];
	}
	for (int i = 0; i < this->m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			p[i][j] = m.p[i][j];
		}
	}
	return *this;
}
matrix& operator +=(matrix& m1, matrix& m2)
{
	if (m1.getCol() != m2.getCol() || m1.getRow() != m2.getRow())
	{
		throw exception("Error!!!   The number of rows and columns does not match.");
	}
	m1 = m1 + m2;
	return m1;
}
matrix& operator -=(matrix& m1, matrix& m2)
{
	if (m1.getCol() != m2.getCol() || m1.getRow() != m2.getRow())
	{
		throw exception("Error!!!   The number of rows and columns does not match.");
	}
	m1 = m1 - m2;
	return m1;
}
matrix operator *(double n, matrix& m)
{
	matrix m1(m);
	for (int i = 0; i < m1.getRow(); i++)
	{
		for (int j = 0; j < m1.getCol(); j++)
		{
			m1(i, j) *= n;
		}
	}
	return m1;
}
matrix operator *(matrix& m, double n)
{
	matrix m1(m);
	for (int i = 0; i < m1.getRow(); i++)
	{
		for (int j = 0; j < m1.getCol(); j++)
		{
			m1(i, j) *= n;
		}
	}
	return m1;
}
matrix operator /(matrix& m, double n)
{
	matrix m1(m);
	for (int i = 0; i < m1.getRow(); i++)
	{
		for (int j = 0; j < m1.getCol(); j++)
		{
			m1(i, j) /= n;
		}
	}
	return m1;
}
matrix operator *(matrix& m1, matrix& m2)
{
	if (m1.getCol() != m2.getRow())
	{
		throw exception("Error!!!   The number of rows and columns does not match.");
	}
	matrix m3(m1.getRow(), m2.getCol());

	for (int i = 0; i < m3.getRow(); i++)
	{
		for (int j = 0; j < m3.getCol(); j++)
		{
			m3(i, j) = 0;
			for (int k = 0; k < m1.getCol(); k++)
			{
				m3(i, j) += m1(i, k) * m2(k, j);
			}
		}
	}

	return m3;
}
matrix& operator *=(matrix& m1, matrix& m2)
{
	if (m1.getCol() != m2.getRow())
	{
		throw exception("Error!!!   The number of rows and columns does not match.");
	}
	m1 = m1 * m2;
	return m1;
}
bool operator ==(matrix& m1, matrix& m2)
{
	if (m1.getCol() != m2.getCol() || m1.getRow() != m2.getRow())
	{
		return false;
	}
	else
	{
		for (int i = 0; i < m1.getRow(); i++)
		{
			for (int j = 0; j < m1.getCol(); j++)
			{
				if (m1(i, j) == m2(i, j))
					continue;
				else
					return false;
			}
		}
		return true;
	}
}
bool operator !=(matrix& m1, matrix& m2)
{
	return !(m1 == m2);
}
matrix& operator *=(matrix& m1, double n)
{
	m1 = m1 * n;
	return m1;
}
matrix& operator /=(matrix& m1, double n)
{
	m1 = m1 / n;
	return m1;
}
istream& operator >>(istream& in, matrix& m)
{
	for (int i = 0; i < m.getRow(); i++)
	{
		for (int j = 0; j < m.getCol(); j++)
		{
			in >> m(i, j);
		}
	}
	return in;
}
ostream& operator <<(ostream& out, matrix& m)
{
	for (int i = 0; i < m.getRow(); i++)
	{
		for (int j = 0; j < m.getCol(); j++)
		{
			out << m(i, j) << "\t";
		}
		out << endl;
	}
	return out;
}