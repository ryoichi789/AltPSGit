#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "superposer.h"
#include "position.h"
#include "matrix.h"
#include "my_math.h"

using namespace Cba;
using std::vector;

// superpose:
double Superposer::superpose(vector<Position*>& a, vector<Position*>& b)
{
	m_translation1.clear();
	m_translation2.clear();
	m_rotation.clear();
	int npoints = std::min(a.size(), b.size());
	if (npoints > 0) {
		Position ca = calc_center(npoints, a);
		Position cb = calc_center(npoints, b);
		// translation vectors
		m_translation1 -= ca;
		m_translation2 += cb;
		// move ceters of a and b to origin
		for (int i = 0; i < a.size(); i++) *a[i] -= ca;
		for (int i = 0; i < b.size(); i++) *b[i] -= cb;
		if (npoints > 1) {
			// calculate rotation matrix
			calc_rotation(npoints, a, b);
			// rotate a
			for (int i = 0; i < a.size(); i++) a[i]->rotate_by(m_rotation);
		}
		// translate a and b
		for (int i = 0; i < a.size(); i++) *a[i] += m_translation2;
		for (int i = 0; i < b.size(); i++) *b[i] += m_translation2;
	}
	return Superposer::root_mean_square_deviation(a, b);
}

// superpose:
double Superposer::superpose(vector<Position>& a, vector<Position>& b)
{
	vector<Position*> pa;
	for (vector<Position>::iterator i = a.begin(); i != a.end(); ++i)
		pa.push_back(&(*i));
	vector<Position*> pb;
	for (vector<Position>::iterator i = b.begin(); i != b.end(); ++i)
		pb.push_back(&(*i));
	return superpose(pa, pb);
}

// calc_center:
Position Superposer::calc_center(size_t npoints, vector<Position*>& ps) const
{
	Position c;
	for (size_t i = 0; i < npoints; i++) c += *ps[i];
	c /= npoints;
	return c;
}

// calc_rotation:
// using Kabsch's algorithm (1978) for best superposition
void Superposer::calc_rotation(int npoints,
	vector<Position*>& p1, vector<Position*>& p2)
{
	Matrix r; // rij = Sum_over_k{ p2(k,i)*p1(k,j) }
	for (int k = 0; k < npoints; k++) {
		r.elem(0,0) += p2[k]->x() * p1[k]->x();
		r.elem(0,1) += p2[k]->x() * p1[k]->y();
		r.elem(0,2) += p2[k]->x() * p1[k]->z();
		r.elem(1,0) += p2[k]->y() * p1[k]->x();
		r.elem(1,1) += p2[k]->y() * p1[k]->y();
		r.elem(1,2) += p2[k]->y() * p1[k]->z();
		r.elem(2,0) += p2[k]->z() * p1[k]->x();
		r.elem(2,1) += p2[k]->z() * p1[k]->y();
		r.elem(2,2) += p2[k]->z() * p1[k]->z();
	}
	Matrix trr = r.transpose() * r;

	// for tRR, calculates eigenvalues mu(k) and mutually orthogonal eigenvectors
	Matrix eigen = MyMath::get_eigen(trr);
	// eigenvalues
	double mu[THREE_D];
	for (int i = 0; i < THREE_D; i++) mu[i] = eigen.elem(i, 0);
	// eigenvectors
	double a[THREE_D][THREE_D];
	for (int k = 0; k < THREE_D; k++)
		for (int i = 0; i < THREE_D; i++)
			a[k][i] = eigen.elem(i, k+1);
			// a[k][i] is the i-th component of k-th vector
			// (to follow Kabsch(1976,1978)'s notation)
	make_right_handed(a);
		// sets a[2] = a[0] x a[1] to be sure to have a right-handed system

	// determines Ra[k], and normalizes the first two vectors
	// to obtain b[0], b[1], b[2]  (b[k]=Ra[k]/mu[k], k=1,2)
	double b[THREE_D][THREE_D];
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < THREE_D; j++) {
			b[i][j] = 0.0;
			for (int k = 0; k < THREE_D; k++)
				b[i][j] += r.elem(j,k) * a[i][k];
		}
		// normalize b[i]
		double w = 0.0;
		for (int j = 0; j < THREE_D; j++) w += b[i][j] * b[i][j];
		double t = std::sqrt(1.0/w);
		for (int j = 0; j < THREE_D; j++) b[i][j] *= t;
	}
	calc_vector_product(b[0], b[1], b[2]); // b[2] = b[0] x b[1]
	set_rotation(a, b); // rotation matrix rot[i][j] = b[k][i]*a[k][j]
}

// make_right_handed:
// sets a[2] = a[0] x a[1] to be sure to have a right-handed system
void Superposer::make_right_handed(double a[THREE_D][THREE_D])
{
	double v1[THREE_D], v2[THREE_D], v3[THREE_D];
	for (int i = 0; i < THREE_D; i++) v1[i] = a[0][i];
	for (int i = 0; i < THREE_D; i++) v2[i] = a[1][i];
	calc_vector_product(v1, v2, v3);
	for (int i = 0; i < THREE_D; i++) a[2][i] = v3[i];
}

// calc_vector_product:
void Superposer::calc_vector_product(
	const double v1[THREE_D], const double v2[THREE_D], double v3[THREE_D])
{
	v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
	v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
	v3[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

// set_rotation:
void Superposer::set_rotation(
	const double a[THREE_D][THREE_D], const double b[THREE_D][THREE_D])
{
	Matrix rot(THREE_D, THREE_D);
	for (int i = 0; i < THREE_D; i++) {
		for (int j = 0; j < THREE_D; j++) {
			for (int k = 0; k < THREE_D; k++)
				rot.elem(i,j) += b[k][i] * a[k][j];
		}
	}
	m_rotation = rot;
}

// transform:
void Superposer::transform(Position& p) const
{
	p += m_translation1;
	p.rotate_by(m_rotation);
	p += m_translation2;
}

// mean_square_deviation:
double Superposer::mean_square_deviation(
	const vector<Position*>& a, const vector<Position*>& b)
{
	size_t n = std::min(a.size(), b.size());
	if (n == 0) return 0;
	double msd = 0;
	for (size_t i = 0; i < n; i++)
		msd += a[i]->square_distance_from(*b[i]);
	return (msd/n);
}

// root_mean_square_deviation:
double Superposer::root_mean_square_deviation(
	const vector<Position*>& a, const vector<Position*>& b)
{
	return std::sqrt(mean_square_deviation(a, b));
}

// print:
void Superposer::print(std::ostream& os) const
{
	using std::endl;
	using std::fixed;
	using std::setprecision;
	using std::setw;

	os << std::fixed << std::setprecision(3) << std::setw(8) << m_translation1.x() << ' ';
	os << std::fixed << std::setprecision(3) << std::setw(8) << m_translation1.y() << ' ';
	os << std::fixed << std::setprecision(3) << std::setw(8) << m_translation1.z() << endl;

	for (int i = 0; i < THREE_D; i++) {
		for (int j = 0; j < THREE_D; j++) {
			os << std::fixed << std::setprecision(3) << std::setw(8)
				<< m_rotation.elem(i,j);
			if (j == (THREE_D-1)) os << endl;
			else os << ' ';
		}
	}

	os << std::fixed << std::setprecision(3) << std::setw(8) << m_translation2.x() << ' ';
	os << std::fixed << std::setprecision(3) << std::setw(8) << m_translation2.y() << ' ';
	os << std::fixed << std::setprecision(3) << std::setw(8) << m_translation2.z() << endl;
}
