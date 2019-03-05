#ifndef CBA_MOL_FIN_H
#define CBA_MOL_FIN_H

#include <vector>
#include <string>
//#include <list>
#include <map>
#include <iostream>
#include <valarray>
#include "matrix.h"
#include "atom.h"
#include "molecule.h"
#include "bit_string.h"

namespace Cba {

class MolFin {
public:
	MolFin();
	~MolFin();

	const std::string& id() const;
	void set_id(const std::string& id);
	void make(Molecule& mol);
	int natoms() const;
	int size() const;
	int at(size_t i) const;
	void clear();

	double rmsd(const MolFin& mf) const;
	double dist(const MolFin& mf) const;
	double tc(const MolFin& mf) const; // Tanimoto coefficient = sum of ai*bi / sum of (ai*ai + bi*bi - ai*bi)
	static double pvalue(double tc);

	double tc2(const MolFin& mf) const; // Tanimoto coefficient = sum of min(ai,bi)/sum of (ai + bi -min(ai,bi))
	//double tcw(const MolFin& mf) const; // weighted Tanimoto coefficient
	double tc_bin(const MolFin& mf) const; // Tanimoto coefficient using bit-string form of fingerprint

	static void get_triad(size_t i, int& c1, int& c2, int& c3, int& d12, int& d13, int& d23);

	void write(std::ostream& os) const;
	bool read(std::istream& is);

	void debug_on();
	void debug_off();
	int max_count() const;

private:
	std::string m_id;
	std::vector<int> m_fin;
	int m_natoms;
	bool m_debug;
	int m_csqr; // sum of square of components; for calculating Tanimoto coefficient

	/*
	static double m_weights[];
		// weight for each component:
		// wi = 1/((pi/p)+0.01), where
		// p = 1/size() (mean probability)
		// pi = observed probability of triad i
	static std::valarray<double> m_weights_sqr;
	*/

	static const int DISTANCE_RANGE_MIN;
	static const int DISTANCE_RANGE_MAX;
	static const double p_values[];
	static const double P_VALUES_TC_STEP;

	// Six atom classes (donor/acceptor/polar/hydrophobic_aromatic/hydrophobic_non_aromatic/none) are used:
	//
	// cation =>    donor
	// anion =>     acceptor
	// donor =>     donor
	// acceptor =>  acceptor
	// polar =>     polar
	// hydrophobic => hydrophobic_aromatic, hydrophobic_non_aromatic
	// none =>      none

	enum AtomTypes {
		UNDEFINED = 0,
		DONOR,
		ACCEPTOR,
		POLAR,
		HYDROPHOBIC_AROMATIC,
		HYDROPHOBIC_NON_AROMATIC,
		NONE,
		NO_OF_ATOM_TYPES
	};

	// Six distances are used:
	// graph distance 1 => 1
	// graph distance 2 => 2
	// graph distance 3 => 3
	// graph distance 4 => 4
	// graph distance 5-9 => 5
	// graph distance 10- => 6

	struct AtomTriadLess;
	struct AtomTriad {
		int c1; // class of atom1
		int c2; // class of atom2
		int c3; // class of atom3
		int d12; // distance atom1-atom2
		int d13; // distance atom1-atom3
		int d23; // distance atom2-atom3
			// c1 <= c2 <= c3
			// when all atom classes are identical;
			//  d12 <= d13 <= d23
			// when two of atom classes are identical;
			//  for c1 = c2 < c3, d13 <= d23
			//  for c1 < c2 = c3 , d12 <= d13

		int bit_pos() const; // returns the position of the bit corresponding to this atom triad
			// returns -1 if no relevant bit found

		void print(std::ostream& os) const;

		AtomTriad(int a1 = 0, int a2 = 0, int a3 = 0, int t12 = 0, int t13 = 0, int t23 = 0);
		void set(int a1, int a2, int a3, int t12, int t13, int t23);
		void set_as_is(int a1, int a2, int a3, int t12, int t13, int t23);
		void clear();
		static void set_bit_pos_tbl();
		static bool satisfy_triangle_inequality(int t12, int t13, int t23);
		static std::map<AtomTriad, int, AtomTriadLess> bit_pos_tbl; // positions of the bits corresponding to individual atom triads in a fingerprint
		static std::vector<AtomTriad> triad_at_pos; // atom triad at position i of fingerprint
	};

	struct AtomTriadLess {
		bool operator()(const AtomTriad& t1, const AtomTriad& t2) const;
	};

	AtomTriad m_atom_triad;

	double dsqr(const MolFin& mf) const;

	void calc_csqr();
	static int get_atom_class(const Atom* a);
	static int get_atom_dij(Molecule& mol, int i, int j);
};

inline const std::string& MolFin::id() const { return m_id; }
inline void MolFin::set_id(const std::string& id) { m_id = id; }
inline int MolFin::natoms() const { return m_natoms; }
inline int MolFin::size() const { return m_fin.size(); }
inline int MolFin::at(size_t i) const { return m_fin.at(i); }
inline void MolFin::debug_on() { m_debug = true; }
inline void MolFin::debug_off() { m_debug = false; }

}
#endif
