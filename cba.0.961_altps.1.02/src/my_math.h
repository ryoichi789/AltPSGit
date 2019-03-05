#ifndef CBA_MY_MATH_H
#define CBA_MY_MATH_H

#include <valarray>
#include "matrix.h"
#include "exception.h"

namespace Cba {

class MyMath {
public:
	static Matrix get_eigen(const Matrix& mat) throw(Exception);
		// assumes mat is a symmetric nxn matrix;
		// returns nx(n+1) matrix:
		//	first column contains eigenvalues,
		//	the other columns contain eigenvectors;

private:
	// to implement get_eigen()
	static const double EIGVAL_EPS;
		// relative error tolerance for eigenvalues
	static void eigen_house(Matrix& a, Matrix& work);
		// a = nxn matrix (copy of this matrix), work = nx6 matrix
	static void eigen_bisec(std::valarray<double>& eigval, Matrix& work);
	static void eigen_invitr(Matrix& a, std::valarray<double>& eigval, Matrix& eigvec,
			Matrix& work, std::valarray<int>& lex);
	static void eigen_ludecomp(int n, double small,
			std::valarray<double>& di, std::valarray<double>& bl,
			std::valarray<double>& bu, std::valarray<double>& bv,
			std::valarray<double>& cm, std::valarray<int>& lex);
	static void eigen_initvec(Matrix& eigvec);
};

}
#endif
