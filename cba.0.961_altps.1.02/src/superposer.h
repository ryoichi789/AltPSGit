#ifndef CBA_SUPERPOSER_H
#define CBA_SUPERPOSER_H
/*
 * superposer.h
 *	class for superposing a series of points A (typically, a molecule)
 *	onto another series of points B (another molecule);
 *
 *	The basic procedure is as follows:
 *
 *	(1) The class first receives the same number of points (two series of
 *	points), A = {a1,a2,...,an} and B = {b1,b2,...,bn};
 *
 *	(2) It then calculate 'superposition matrices' (translation 1,
 *	rotation, and translation 2) to minimize RMSD between the two
 *	series of points after superposition (A is superposed onto B);
 *	translation 1 translates A to the origin (A->A'); 
 *	rotation rotates A' (A'->A'');
 *	translation 2 then translates A'' to the center of B (A''->A''');
 *
 *	(3) Using the superposition matrices, it can transform coordinates
 *	of any given set of points.
 *
 *	Typically, a set of atoms (a molecule) are superposed onto another
 *	set of atoms (another molecule) using the superposition matrices
 *	derived from equivalent atoms from the two sets.
 *
 *	Superposition matrices are derived using Kabsch's algorithm (1978).
 */
#include <vector>
#include <iostream>
#include "position.h"
#include "matrix.h"

namespace Cba {

class Superposer {
public:
	Superposer() {};
	double superpose(std::vector<Position*>& a, std::vector<Position*>& b);
		// superposes a onto b;
		// returns r.m.s.d after superposition;
	double superpose(std::vector<Position>& a, std::vector<Position>& b);
	void transform(Position& p) const;
		// transforms coordinates of p
		//	using the superposition matrix given by superpose();
	static double mean_square_deviation(
		const std::vector<Position*>& a, const std::vector<Position*>& b);
	static double root_mean_square_deviation(
		const std::vector<Position*>& a, const std::vector<Position*>& b);
	void print(std::ostream& os) const; // prints superposition matrix

private:
	// superposition matrices
	Position m_translation1; // translation 1
	Matrix m_rotation; // rotation
	Position m_translation2; // translation 2

	static const int THREE_D = 3;
	Position calc_center(size_t n_points, std::vector<Position*>& ps) const;
	void calc_rotation(int n_points,
		std::vector<Position*>& p1, std::vector<Position*>& p2);
	void make_right_handed(double a[THREE_D][THREE_D]);
	void calc_vector_product(const double v1[THREE_D],
		const double v2[THREE_D], double v3[THREE_D]); // v3 = v1 x v2
	void set_rotation(
		const double a[THREE_D][THREE_D], const double b[THREE_D][THREE_D]);
};

}
#endif
