#include <valarray>
#include <cmath>
#include <limits>
#include "my_math.h"
#include "matrix.h"
#include "exception.h"

using namespace Cba;
using std::valarray;
using std::numeric_limits;

const double MyMath::EIGVAL_EPS = 1.0e-12;

//////////////////////////////////////////////////////////
// functions for calculating eigenvalues and eigenvectors
//////////////////////////////////////////////////////////

//
Matrix MyMath::get_eigen(const Matrix& mat) throw(Exception)
{
	if (mat.nrows() != mat.ncolumns())
		throw("matrix is not square symmetric");

	int nrows = mat.nrows();
	Matrix mtmp(mat);
	Matrix w(nrows, 6);
	eigen_house(mtmp, w);

	valarray<double> eigval(nrows);
	eigen_bisec(eigval, w);

	Matrix work(nrows, 8);
	for (int i = 0; i < nrows; i++)
		for (int j = 0; j < 6; j++)
			work.elem(i,j) = w.elem(i,j);
	Matrix eigvec(nrows, nrows);
	valarray<int> lex(nrows);
	eigen_invitr(mtmp, eigval, eigvec, work, lex);

	Matrix em(nrows, nrows+1);
	for (int i = 0; i < nrows; i++)
		em.elem(i, 0) = eigval[i];
	for (int j = 0; j < nrows; j++)
		for (int i = 0; i < nrows; i++)
			em.elem(i, j+1) = eigvec.elem(i, j);
	return em;
}

//
void MyMath::eigen_house(Matrix& a, Matrix& work)
{
	int ndim = a.nrows();
	valarray<double> al(ndim);
	valarray<double> be(ndim);
	valarray<double> co(ndim);
	valarray<double> w(ndim);
	valarray<double> p(ndim);
	valarray<double> q(ndim);

	for (int k = 0; k < (ndim-2); k++) {
		double t = 0.0;
		for (int i = k+1; i < ndim; i++)
			t += a.elem(i,k) * a.elem(i,k);
		double s = sqrt(t);
		if (a.elem(k+1,k) < 0.0) s = -s;
		al[k] = a.elem(k,k);
		be[k] = -s;
		if (t == 0.0) continue;
		co[k] = 1.0 / (t + a.elem(k+1,k)*s);
		w[k+1] = a.elem(k+1,k) + s;
		for (int i = k+2; i < ndim; i++) w[i] = a.elem(i,k);
		a.elem(k+1,k) = w[k+1];
		for (int i = k+1; i < ndim; i++) {
			t=0.0;
			for (int j = k+1; j <= i; j++)
				t += (a.elem(i,j)*w[j]);
			for (int j = i+1; j < ndim; j++)
				t += (a.elem(j,i)*w[j]);
			p[i] = co[k]*t;
		}
		t = 0.0;
		for (int i = k+1; i < ndim; i++)
			t += p[i]*w[i];
		s = 0.5*(co[k]*t);
		for (int i = k+1; i < ndim; i++) q[i] = p[i] - s*w[i];
		for (int j = k+1; j < ndim; j++)
			for (int i = j; i < ndim; i++)
				a.elem(i,j) -= w[i]*q[j] + q[i]*w[j];
	}
	al[ndim-2] = a.elem(ndim-2,ndim-2);
	al[ndim-1] = a.elem(ndim-1,ndim-1);
	be[ndim-2] = a.elem(ndim-1,ndim-2);

	for(int i = 0; i < ndim; i++) work.elem(i,0) = al[i];
	for(int i = 0; i < ndim; i++) work.elem(i,1) = be[i];
	for(int i = 0; i < ndim; i++) work.elem(i,2) = co[i];
	for(int i = 0; i < ndim; i++) work.elem(i,3) = w[i];
	for(int i = 0; i < ndim; i++) work.elem(i,4) = p[i];
	for(int i = 0; i < ndim; i++) work.elem(i,5) = q[i];
}

//
void MyMath::eigen_bisec(valarray<double>& eigval, Matrix& work)
{
	int n = work.nrows();
	valarray<double> al(n), be(n), b2(n);
	for (int i = 0; i < n; i++) al[i] = work.elem(i,0);
	for (int i = 0; i < n; i++) be[i] = work.elem(i,1);
	for (int i = 0; i < n; i++) b2[i] = work.elem(i,3);
	be[n-1] = 0.0;
	double range = fabs(al[0]) + fabs(be[0]);
	for (int i = 1; i < n; i++) {
		double range1 = fabs(be[i-1]) + fabs(al[i]) + fabs(be[i]);
		if (range1 > range) range = range1;
	}
	//if (nls < 0) range = -range; // nls is fixed at 1 in this implementation
	b2[0] = 0.0;
	for (int i = 1; i < n; i++) b2[i] = be[i-1] * be[i-1];
	double eps = EIGVAL_EPS;
	double epsa = fabs(range) * eps;
	double small = epsa * eps;
	for (int i = 0; i < n; i++) eigval[i] = -range;
	double b = range;
	for (int k = 0; k < n; k++) {
		double a = eigval[k];
	loop:
		double c = 0.5 * (a+b);
		if ((fabs(b-a) > epsa) && (c != a) && (c != b)) {
			int nneg = 0;
			double  g = 1.0;
			int ipass = 0;
			for (int i = 0; i < n; i++) {
				if (ipass == 0) g = c - al[i] - b2[i]/g;
				else if (ipass == 1) ipass = 2;
				else { g = c - al[i]; ipass = 0;}
				if (ipass==0) {
					if (g <= 0.0) ++nneg;
					if (fabs(g) <= fabs(b2[i]*small)) ipass = 1;
				}
			}
			// if (nls < 0) nneg = n - nneg; // nls is fixed at 1 in this implementation
			if (nneg <= k) b = c;
			else {
				a = c;
				for (int i = k; i < ((nneg < n) ? nneg : n); i++) eigval[i] = c;
			}
			goto loop;
		}
	}

	for (int i = 0; i < n; i++) work.elem(i,0) = al[i];
	for (int i = 0; i < n; i++) work.elem(i,1) = be[i];
	for (int i = 0; i < n; i++) work.elem(i,3) = b2[i];
}

//
void MyMath::eigen_invitr(Matrix& a, valarray<double>& eigval, Matrix& eigvec,
	Matrix& work, valarray<int>& lex)
{
	int n = a.nrows();
	valarray<double> al(n), be(n), co(n);
	valarray<double> di(n), bl(n), bu(n), bv(n), cm(n);
	for (int i = 0; i < n; i++) al[i] = work.elem(i,0);
	for (int i = 0; i < n; i++) be[i] = work.elem(i,1);
	for (int i = 0; i < n; i++) co[i] = work.elem(i,2);
	for (int i = 0; i < n; i++) di[i] = work.elem(i,3);
	for (int i = 0; i < n; i++) bl[i] = work.elem(i,4);
	for (int i = 0; i < n; i++) bu[i] = work.elem(i,5);
  
	int km, lend;
	double s, t;
	double epsmac = numeric_limits<double>::epsilon();
	double eps = epsmac *
		(fabs(eigval[0]) > fabs(eigval[n-1]) ?
			fabs(eigval[0]) : fabs(eigval[n-1]));
	double small = eps;
	eigen_initvec(eigvec); // initial guess of eigenvectors
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			di[i] = eigval[k] - al[i];
			bl[i] = -be[i];
			bu[i] = -be[i];
		}
		eigen_ludecomp(n, small, di, bl, bu, bv, cm, lex); // LU decomposition
		if (k == 0) km = k;
		else if (fabs(eigval[k]-eigval[km]) > eps) km = k;
		else {
			for (int i = km; i <= k-1; i++) {
				t = 0.0;
				for (int j = 0; j < n; j++)
					t += eigvec.elem(j,i) * eigvec.elem(j,k);
				s = t;
				for (int j = 0; j < n; j++)
					eigvec.elem(j,k) -= s * eigvec.elem(j,i);
			}
		}
		lend = k - km + 2;
		if (lend > 5) lend = 5;
		for (int l = 0; l < lend; l++) {
			if ((l != 0) || (k != km)) {
				for(int i = 0; i < n-1; i++) { // forward substitution
					if (lex[i] == 0) eigvec.elem(i+1,k) -= cm[i+1] * eigvec.elem(i,k);
					else {
						s = eigvec.elem(i,k);
						eigvec.elem(i,k) = eigvec.elem(i+1,k);
						eigvec.elem(i+1,k) = s - cm[i+1] * eigvec.elem(i,k);
					}
				}
			}
			for (int i = n-1; i >= 0; --i) { // backward substitution
				s = eigvec.elem(i,k);
				if (i < n-1) s -= bu[i] * eigvec.elem(i+1,k);
				if (i < n-2) s -= bv[i] * eigvec.elem(i+2,k);
				eigvec.elem(i,k) = s/di[i];
			}

			// normalize vector to avoid overflow
			t = 0.0;
			for (int j = 0; j < n; j++)
				t += eigvec.elem(j,k) * eigvec.elem(j,k);
			s = 0.0;
			if (t != 0.0) s = sqrt(1.0/t);
			for (int j = 0; j < n; j++) eigvec.elem(j,k) *= s;
		}
	}

	// transformation of eigenvectors to original space
	for (int k = 0; k < n; k++) {
		for (int i = n-3; i >= 0; i--) {
			t = 0.0;
			for (int j = i+1; j < n; j++)
				t += a.elem(j,i) * eigvec.elem(j,k);
			s = t * co[i];
			for (int j = i+1; j < n; j++)
				eigvec.elem(j,k) -= s * a.elem(j,i);
		}
	}

	// orthogonalize eigen vectors if degenerate
	km = 0;
	for (int k = 1; k < n; k++) {
		if (fabs(eigval[k]-eigval[km]) >= eps) km = k;
		else {
			for (int i = km; i <= k-1; i++) {
				t=0.0;
				for (int j = 0; j < n; j++)
					t += eigvec.elem(j,i) * eigvec.elem(j,k);
				s = t;
				for (int j = 0; j < n; j++)
					eigvec.elem(j,k) -= s * eigvec.elem(j,i);
			}
			t = 0.0;
			for (int j = 0; j < n; j++)
				t += eigvec.elem(j,k) * eigvec.elem(j,k);
			s = sqrt(1.0/t);
			for (int j = 0; j < n; j++)
				eigvec.elem(j,k) *= s;
		}
	}
}

//
void MyMath::eigen_ludecomp(int n, double small,
	valarray<double>& di, valarray<double>& bl, 
	valarray<double>& bu, valarray<double>& bv,
	valarray<double>& cm, valarray<int>& lex)
{
	for (int i = 0; i < n-1; i++) {
		if (fabs(di[i]) >= fabs(bl[i])) {
			lex[i] = 0;
			if (fabs(di[i]) < small) di[i] = small;
			cm[i+1] = bl[i]/di[i];
			di[i+1] = di[i+1] - cm[i+1]*bu[i];
			bv[i] = 0.0;
		} else {
			lex[i] = 1;
			cm[i+1] = di[i]/bl[i];
			di[i] = bl[i];
			double s = bu[i];
			bu[i] = di[i+1];
			bv[i] = bu[i+1];
			di[i+1] = s - cm[i+1]*bu[i];
			bu[i+1] = -cm[i+1]*bv[i];
		}
	}
	if (fabs(di[n-1]) < small) di[n-1] = small;
}

//
void MyMath::eigen_initvec(Matrix& eigvec)
{
	int n = eigvec.nrows();
	int ml = 1664501;
	int ia = 1229;
	int ic = 351750;
	double rnorm = 1.0/(double)ml;
	int ir = 1;
	int nr = n * n;
	valarray<double> r(nr);
	for (int i = 0; i < nr; i++) {
		ir = (ia*ir + ic) % ml;
		r[i] = (double)ir;
	}

	int k = 0;
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			eigvec.elem(i,j) = rnorm * r[k++];
}
