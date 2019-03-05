#ifndef CBA_BOND_MAKER_H
#define CBA_BOND_MAKER_H

#include <vector>
#include <map>
#include "molecule.h"
#include "protein.h"
#include "atom.h"

namespace Cba {

// BondMaker:
// defines covalent bonds between atoms of a small molecule

class BondMaker {
public:
	BondMaker();
	~BondMaker();

	int make_bonds(Molecule* mol);
	int make_bonds(Molecule& mol);

private:
	Molecule* m_mol;

	static const double MAX_BOND_LENGTH;
	static const double BOND_LENGTH_ALLOWANCE;

	static const double C_C_1;
	static const double C_C_2;
	static const double C_C_3;
	static const double C_C_5;
	static const double C_C_6;
	static const double C_N_1;
	static const double C_N_2;
	static const double C_N_3;
	static const double C_O_1;
	static const double C_O_2;
	static const double C_P_1;
	static const double C_S_1;
	static const double C_S_2;
	static const double N_N_1;
	static const double N_N_2;
	static const double N_O_1;
	static const double N_O_2;
	static const double N_P_1;
	static const double N_S_1;
	static const double O_O_1;
	static const double O_P_1;
	static const double O_P_2;
	static const double O_S_1;
	static const double O_S_2;
	static const double P_S_1;
	static const double S_S_1;

	double m_cc1_min;
	double m_cc2_min;
	double m_cn1_min;
	double m_cn2_min;
	double m_co1_min;
	double m_cs1_min;
	double m_nn1_min;
	double m_no1_min;
	double m_op1_min;
	double m_os1_min;
	std::vector<std::vector<Atom*> > m_ar_rings;
	std::map<Atom*, int> m_ar_nbonds_left;
	struct Bond {
		Atom* a1;
		Atom* a2;
		double d;
		int t; // 1, 2, or 3
	};
	std::vector<Bond> m_bonds;

	int make_bonds_for_protein(Protein* prot);
	int define_bond(Atom* a1, Atom* a2, double d);
	int define_bond_c_x(Atom* a1, Atom* a2, double d);
	int define_bond_c_c(double d);
	int define_bond_c_n(double d);
	int define_bond_c_o(double d);
	int define_bond_c_p(double d);
	int define_bond_c_s(double d);
	int define_bond_n_x(Atom* a1, Atom* a2, double d);
	int define_bond_n_n(double d);
	int define_bond_n_o(double d);
	int define_bond_n_p(double d);
	int define_bond_n_s(double d);
	int define_bond_o_x(Atom* a1, Atom* a2, double d);
	int define_bond_o_o(double d);
	int define_bond_o_p(double d);
	int define_bond_o_s(double d);
	int define_bond_p_x(Atom* a1, Atom* a2, double d);
	int define_bond_p_s(double d);
	int define_bond_s_x(Atom* a1, Atom* a2, double d);
	int define_bond_s_s(double d);
	Bond* get_bond(Atom* a1, Atom* a2);

	void make_bonds_consistent();
	void make_bonds_amide();
	void make_bonds_carboxyl();
	void make_bonds_phosphate();
	void make_bonds_nitro();

	void make_bonds_ar();
	void get_ar_rings();
	void get_ar_nbonds_left();
	void make_bonds_ar(std::vector<Atom*>& ring);
	void get_bonds_in_ar(std::vector<Atom*>& ring, std::vector<Bond>& bonds);
	int set_ar_bond0(std::vector<Bond>& bonds);
	void set_ar_bonds(int i0, std::vector<Bond>& bonds);
	void refine_ar_bonds(std::vector<Bond>& bonds);
	bool is_ar6(std::vector<Atom*>& ring);
	bool is_ar5(std::vector<Atom*>& ring);
	bool is_ar(std::vector<Atom*>& ring, int n, double dar);
};

inline bool BondMaker::is_ar6(std::vector<Atom*>& ring)
{ return is_ar(ring, 6, C_C_6); }

inline bool BondMaker::is_ar5(std::vector<Atom*>& ring)
{ return is_ar(ring, 5, C_C_5); }

}
#endif
