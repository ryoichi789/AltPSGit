#include <valarray>
#include <limits>//For CentOS6.2
#include "matrix.h"
#include "exception.h"

using namespace Cba;
using std::valarray;
using std::numeric_limits;

//
Matrix::Matrix(size_t nrows, size_t ncolumns)
	: m_nrows(nrows), m_ncolumns(ncolumns)
{
	m_elems = new valarray<double>(m_nrows*m_ncolumns);
	for (int i = 0; i < m_nrows*m_ncolumns; i++) (*m_elems)[i] = 0;
}

//
Matrix::Matrix(const Matrix& m)
	: m_nrows(m.m_nrows), m_ncolumns(m.m_ncolumns)
{
	m_elems = new valarray<double>(m.m_nrows*m.m_ncolumns);
	*m_elems = *(m.m_elems);
}

// operator=:
Matrix& Matrix::operator=(const Matrix& m)
{
	if (this != &m) {
		if (m_nrows != m.m_nrows || m_ncolumns != m.m_ncolumns) {
			delete m_elems;
			m_nrows = m.m_nrows;
			m_ncolumns = m.m_ncolumns;
			m_elems = new valarray<double>(m.m_nrows*m.m_ncolumns);
		}
		*m_elems = *(m.m_elems);
	}
	return *this;
}

//
Matrix::~Matrix()
{
	delete m_elems;
}

// clear:
void Matrix::clear()
{ for (int i = 0; i < m_nrows*m_ncolumns; i++) (*m_elems)[i] = 0.0; }

//
Matrix Matrix::transpose() const
{
	Matrix tm(m_ncolumns, m_nrows);
	for (int i = 0; i < m_nrows; i++)
		for (int j = 0; j < m_ncolumns; j++)
			tm.elem(j, i) = elem(i, j);
	return tm;
}

//////////////////////////////
// helper functions
//////////////////////////////

// operator+:
Matrix Cba::operator+(const Matrix& a, const Matrix& b) throw(Exception)
{
	if (a.nrows() != b.nrows())
		throw("two matrices have different numbers of rows");
	if (a.ncolumns() != b.ncolumns())
		throw("two matrices have different numbers of columns");

	Matrix c(a.nrows(), a.ncolumns());
	*(c.m_elems) = *(a.m_elems) + *(b.m_elems);
	return c;
}

// operator-:
Matrix Cba::operator-(const Matrix& a, const Matrix& b) throw(Exception)
{
	if (a.nrows() != b.nrows())
		throw("two matrices have different numbers of rows");
	if (a.ncolumns() != b.ncolumns())
		throw("two matrices have different numbers of columns");

	Matrix c(a.nrows(), a.ncolumns());
	*(c.m_elems) = *(a.m_elems) - *(b.m_elems);
	return c;
}

// operator*:
Matrix Cba::operator*(const Matrix& a, const Matrix& b) throw(Exception)
{
	if (a.ncolumns() != b.nrows())
		throw("cannot multiply the two matrices");

	Matrix c(a.nrows(), b.ncolumns());
	for (int i = 0; i < a.nrows(); i++)
		for (int j = 0; j < b.ncolumns(); j++)
			for (int k = 0; k < a.ncolumns(); k++)
				c.elem(i,j) += a.elem(i,k) * b.elem(k,j);
	return c;
}
