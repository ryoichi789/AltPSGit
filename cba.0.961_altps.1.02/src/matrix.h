#ifndef CBA_MATRIX_H
#define CBA_MATRIX_H

#include <valarray>
#include "exception.h"

namespace Cba {

// Matrix represents matrix of real numbers
class Matrix {
public:
	Matrix(size_t nrows = 3, size_t ncolumns = 3);
	Matrix(const Matrix& m);
	Matrix& operator=(const Matrix& m);
	~Matrix();
	void clear();

	int nrows() const;
	int ncolumns() const;
	int nelems() const;
	double& elem(size_t i, size_t j);
	const double& elem(size_t i, size_t j) const;
		// NOTE: elem()'s do no check index range error
	Matrix transpose() const; // returns transposition of the matrix

private:
	size_t m_nrows;
	size_t m_ncolumns;
	std::valarray<double>* m_elems;
	// nxm matrix is represented by valarray<double>(n*m);
	// (i,j)-element is (i*m+j)-th element of the valarray;
	// that is, elements 0,1,2,...,m,...,n of the valarray is
	// mapped to the matrix such as:
	// 0,1,2,...,m-1
	// m,...,2m-1
	// ...
	// (n-1)*m,...,n*m-1

	friend Matrix operator+(const Matrix& a, const Matrix& b) throw(Exception);
	friend Matrix operator-(const Matrix& a, const Matrix& b) throw(Exception);
	friend Matrix operator*(const Matrix& a, const Matrix& b) throw(Exception);
};

inline int Matrix::nrows() const { return m_nrows; }
inline int Matrix::ncolumns() const { return m_ncolumns; }
inline int Matrix::nelems() const { return m_nrows*m_ncolumns; }
inline double& Matrix::elem(size_t i, size_t j)
{ return (*m_elems)[i*m_ncolumns + j]; }
inline const double& Matrix::elem(size_t i, size_t j) const
{ return (*m_elems)[i*m_ncolumns + j]; }

// helper functions
Matrix operator+(const Matrix& a, const Matrix& b) throw(Exception);
Matrix operator-(const Matrix& a, const Matrix& b) throw(Exception);
Matrix operator*(const Matrix& a, const Matrix& b) throw(Exception);
}
#endif
