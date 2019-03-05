#ifndef CBA_ASA_CALCULATOR
#define CBA_ASA_CALCULATOR

#include <vector>
#include <deque>
#include <functional>
#include "atom.h"
#include "position.h"

namespace Cba {

class AsaCalculator {
public:
	AsaCalculator() {}
	~AsaCalculator() {}

	static const double SOLVENT_RADIUS;
	static const double ATOM_SLICING_STEP;
	static const double PI;

	void calc(std::vector<Atom*>& atoms,
		double rsolv = SOLVENT_RADIUS, double dz = ATOM_SLICING_STEP);
		// rsolv: probe radius
		// dz: used in approximating a sphere by a series of circles

private:
	std::vector<double> m_rs; // vdw radius + probe radius of each atom
	double m_rsolv; // radius of solvent (probe)
	double m_rmax; // max atom radius
	double m_zstep;
	std::deque<bool> m_igna;

	struct Arc { // representing an arc in a circle
		double ti; // starting position (in radian)
		double te; // ending position (in radian)
		Arc() : ti(-999.99), te(-999.99) {}
	};

	struct CompArcs : std::binary_function<Arc, Arc, bool> {
		bool operator() (const Arc& a1, const Arc& a2) const;
	};

	// functions used by calc()
	void set_params(const std::vector<Atom*>& atoms, double rsolv, double dz);
	void ignore_atoms_for(size_t iatom, const std::vector<Atom*>& atoms);
	void get_zs(size_t iatom, std::vector<double>& zs,
		const std::vector<Atom*>& atoms);
	double get_accarc_in_circle(
		int iatom, double z, double r, std::vector<Atom*>& atoms);
	Arc cc_overlap(
		double x1, double y1, double r1, double x2, double y2, double r2);
	void combine_arcs(const std::vector<Arc>& arcs0, std::vector<Arc>& arcs);
};

}
#endif
