#include <vector>
#include <map>
#include <cmath>
#include "bond_maker.h"
#include "molecule.h"
#include "protein.h"
#include "atom.h"
#include "position.h"

using namespace Cba;
using std::vector;
using std::map;

const double BondMaker::MAX_BOND_LENGTH = 2.20;
const double BondMaker::BOND_LENGTH_ALLOWANCE = 0.25;

const double BondMaker::C_C_1 = 1.53;
const double BondMaker::C_C_2 = 1.31;
const double BondMaker::C_C_3 = 1.17;
const double BondMaker::C_C_5 = 1.40;
const double BondMaker::C_C_6 = 1.38;
const double BondMaker::C_N_1 = 1.47;
const double BondMaker::C_N_2 = 1.32;
const double BondMaker::C_N_3 = 1.14;
const double BondMaker::C_O_1 = 1.43;
const double BondMaker::C_O_2 = 1.21;
const double BondMaker::C_P_1 = 1.82;
const double BondMaker::C_S_1 = 1.81;
const double BondMaker::C_S_2 = 1.71;
const double BondMaker::N_N_1 = 1.40;
const double BondMaker::N_N_2 = 1.29;
const double BondMaker::N_O_1 = 1.42;
const double BondMaker::N_O_2 = 1.22;
const double BondMaker::N_P_1 = 1.68;
const double BondMaker::N_S_1 = 1.66;
const double BondMaker::O_O_1 = 1.46;
const double BondMaker::O_P_1 = 1.61;
const double BondMaker::O_P_2 = 1.48;
const double BondMaker::O_S_1 = 1.52;
const double BondMaker::O_S_2 = 1.42;
const double BondMaker::P_S_1 = 2.12;
const double BondMaker::S_S_1 = 2.05;

//
BondMaker::BondMaker() : m_mol(0)
{
	m_cc1_min = C_C_1 - 0.5*(C_C_1-C_C_2);
	m_cc2_min = C_C_2 - 0.5*(C_C_2-C_C_3);
	m_cn1_min = C_N_1 - 0.5*(C_N_1-C_N_2);
	m_cn2_min = C_N_2 - 0.5*(C_N_2-C_N_3);
	m_co1_min = C_O_1 - 0.5*(C_O_1-C_O_2);
	m_cs1_min = C_S_1 - 0.5*(C_S_1-C_S_2);
	m_nn1_min = N_N_1 - 0.5*(N_N_1-N_N_2);
	m_no1_min = N_O_1 - 0.5*(N_O_1-N_O_2);
	m_op1_min = O_P_1 - 0.5*(O_P_1-O_P_2);
	m_os1_min = O_S_1 - 0.5*(O_S_1-O_S_2);
}

//
BondMaker::~BondMaker() {}

//
int BondMaker::make_bonds(Molecule* mol)
{
	m_mol = mol;
	m_bonds.clear();
	if (m_mol->nbonds()) return m_mol->nbonds(); // bonds already set
	Protein* prot = dynamic_cast<Protein*>(m_mol);
	if (prot) return make_bonds_for_protein(prot);
	m_mol->remove_hydrogens();

	Atom* a1;
	Atom* a2;
	double d;
	int b;
	Bond bond;
	for (int i = 0; i < m_mol->natoms()-1; i++) {
		a1 = m_mol->atom(i);
		for (int j = i+1; j < m_mol->natoms(); j++) {
			a2 = m_mol->atom(j);
			d = a1->pos.distance_from_neighbor(a2->pos, MAX_BOND_LENGTH);
			if (d < 0.0) continue;
			if (b = define_bond(a1, a2, d)) {
				m_mol->add_bond(a1, a2, b);
				bond.a1 = a1;
				bond.a2 = a2;
				bond.d = d;
				bond.t = b;
				m_bonds.push_back(bond);
			}
		}
	}

	make_bonds_consistent();
	make_bonds_amide();
	make_bonds_carboxyl();
	make_bonds_phosphate();
	make_bonds_nitro();
	make_bonds_ar(); // for aromatic ring

	return m_mol->nbonds();
}

//
int BondMaker::make_bonds(Molecule& mol)
{   
	return make_bonds(&mol);
}

//
int BondMaker::make_bonds_for_protein(Protein* prot)
{
	return 0;
}

//
int BondMaker::define_bond(Atom* a1, Atom* a2, double d)
{
	if (a1->element == Atom::get_element_number("C"))
		return define_bond_c_x(a1, a2, d);
	if (a2->element == Atom::get_element_number("C"))
		return define_bond_c_x(a2, a1, d);
	if (a1->element == Atom::get_element_number("N"))
		return define_bond_n_x(a1, a2, d);
	if (a2->element == Atom::get_element_number("N"))
		return define_bond_n_x(a2, a1, d);
	if (a1->element == Atom::get_element_number("O"))
		return define_bond_o_x(a1, a2, d);
	if (a2->element == Atom::get_element_number("O"))
		return define_bond_o_x(a2, a1, d);
	if (a1->element == Atom::get_element_number("P"))
		return define_bond_p_x(a1, a2, d);
	if (a2->element == Atom::get_element_number("P"))
		return define_bond_p_x(a2, a1, d);
	if (a1->element == Atom::get_element_number("S"))
		return define_bond_s_x(a1, a2, d);
	if (a2->element == Atom::get_element_number("S"))
		return define_bond_s_x(a2, a1, d);
	if (d <= MAX_BOND_LENGTH) return 1;

	return 0;
}

//
int BondMaker::define_bond_c_x(Atom* a1, Atom* a2, double d)
{
	if (a2->element == Atom::get_element_number("C"))
		return define_bond_c_c(d);
	if (a2->element == Atom::get_element_number("N"))
		return define_bond_c_n(d);
	if (a2->element == Atom::get_element_number("O"))
		return define_bond_c_o(d);
	if (a2->element == Atom::get_element_number("P"))
		return define_bond_c_p(d);
	if (a2->element == Atom::get_element_number("S"))
		return define_bond_c_s(d);
	if (d <= MAX_BOND_LENGTH) return 1;
	return 0;
}

//
int BondMaker::define_bond_c_c(double d)
{
	if (d > C_C_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d >= m_cc1_min) return 1;
	if (d >= m_cc2_min) return 2;
	if (d >= C_C_3-BOND_LENGTH_ALLOWANCE) return 3;
	return 0;
}

//
int BondMaker::define_bond_c_n(double d)
{
	if (d > C_N_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d >= m_cn1_min) return 1;
	return 2;
	//if (d >= C_N_3) return 2;
	//if (d < C_N_3) return 3;
	//if (d >= m_cn2_min) return 2;
	//if (d >= C_N_3-BOND_LENGTH_ALLOWANCE) return 3;
	//return 0;
}

//
int BondMaker::define_bond_c_o(double d)
{
	if (d > C_O_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d >= m_co1_min) return 1;
	if (d >= C_O_2-BOND_LENGTH_ALLOWANCE) return 2;
	return 0;
}

//
int BondMaker::define_bond_c_p(double d)
{
	if (d > C_P_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d < C_P_1-BOND_LENGTH_ALLOWANCE) return 0;
	return 1;
}

//
int BondMaker::define_bond_c_s(double d)
{
	if (d > C_S_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d >= m_cs1_min) return 1;
	if (d >= C_S_2-BOND_LENGTH_ALLOWANCE) return 2;
	return 0;
}

//
int BondMaker::define_bond_n_x(Atom* a1, Atom* a2, double d)
{
	if (a2->element == Atom::get_element_number("N"))
		return define_bond_n_n(d);
	if (a2->element == Atom::get_element_number("O"))
		return define_bond_n_o(d);
	if (a2->element == Atom::get_element_number("P"))
		return define_bond_n_p(d);
	if (a2->element == Atom::get_element_number("S"))
		return define_bond_n_s(d);
	if (d <= MAX_BOND_LENGTH) return 1;
	return 0;
}

//
int BondMaker::define_bond_n_n(double d)
{
	if (d > N_N_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d >= m_nn1_min) return 1;
	if (d >= N_N_2-BOND_LENGTH_ALLOWANCE) return 2;
	return 0;
}

//
int BondMaker::define_bond_n_o(double d)
{
	if (d > N_O_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d >= m_no1_min) return 1;
	if (d >= N_O_2-BOND_LENGTH_ALLOWANCE) return 2;
	return 0;
}

//
int BondMaker::define_bond_n_p(double d)
{
	if (d > N_P_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d < N_P_1-BOND_LENGTH_ALLOWANCE) return 0;
	return 1;
}

//
int BondMaker::define_bond_n_s(double d)
{
	if (d > N_S_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d < N_S_1-BOND_LENGTH_ALLOWANCE) return 0;
	return 1;
}

//
int BondMaker::define_bond_o_x(Atom* a1, Atom* a2, double d)
{
	if (a2->element == Atom::get_element_number("O"))
		return define_bond_o_o(d);
	if (a2->element == Atom::get_element_number("P"))
		return define_bond_o_p(d);
	if (a2->element == Atom::get_element_number("S"))
		return define_bond_o_s(d);
	if (d <= MAX_BOND_LENGTH) return 1;
	return 0;
}

//
int BondMaker::define_bond_o_o(double d)
{
	if (d > O_O_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d < O_O_1-BOND_LENGTH_ALLOWANCE) return 0;
	return 1;
}

//
int BondMaker::define_bond_o_p(double d)
{
	if (d > O_P_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d >= m_op1_min) return 1;
	if (d >= O_P_2-BOND_LENGTH_ALLOWANCE) return 2;
	return 0;
}

//
int BondMaker::define_bond_o_s(double d)
{
	if (d > O_S_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d >= m_os1_min) return 1;
	if (d >= O_S_2-BOND_LENGTH_ALLOWANCE) return 2;
	return 0;
}

//
int BondMaker::define_bond_p_x(Atom* a1, Atom* a2, double d)
{
	if (a2->element == Atom::get_element_number("S"))
		return define_bond_p_s(d);
	if (d <= MAX_BOND_LENGTH) return 1;
	return 0;
}

//
int BondMaker::define_bond_p_s(double d)
{
	if (d > P_S_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d < P_S_1-BOND_LENGTH_ALLOWANCE) return 0;
	return 1;
}

//
int BondMaker::define_bond_s_x(Atom* a1, Atom* a2, double d)
{
	if (a2->element == Atom::get_element_number("S"))
		return define_bond_p_s(d);
	if (d <= MAX_BOND_LENGTH) return 1;
	return 0;
}

//
int BondMaker::define_bond_s_s(double d)
{
	if (d > S_S_1+BOND_LENGTH_ALLOWANCE) return 0;
	if (d < S_S_1-BOND_LENGTH_ALLOWANCE) return 0;
	return 1;
}

//
BondMaker::Bond* BondMaker::get_bond(Atom* a1, Atom* a2)
{
	for (int i = 0; i < m_bonds.size(); i++) {
		if ((a1 == m_bonds[i].a1 && a2 == m_bonds[i].a2) ||
			(a1 == m_bonds[i].a2 && a2 == m_bonds[i].a1)) {
			return &(m_bonds[i]);
		}
	}
	return 0;
}

//
void BondMaker::make_bonds_consistent()
{
	Atom* ac;
	vector<Atom*> aa;
	for (int i = 0; i < m_mol->natoms(); i++) {
		ac = m_mol->atom(i);
		if (ac->element != Atom::get_element_number("C")) continue;
		m_mol->get_atoms_bonded_to(ac, aa);
		if (aa.size() == 4) {
			for (int j = 0; j < aa.size(); j++) {
				if (m_mol->bond_between(ac, aa[j]) >= 2)
					m_mol->add_bond(ac, aa[j], 1);
			}
		} else if (aa.size() == 3) {
			for (int j = 0; j < aa.size(); j++) {
				if (m_mol->bond_between(ac, aa[j]) >= 3)
					m_mol->add_bond(ac, aa[j], 1);
			}
		}
	}
}

//
void BondMaker::make_bonds_amide()
{
	Atom* ac;
	Atom* ao;
	Atom* an;
	Atom* ax;
	vector<Atom*> aa;
	vector<Atom*> ab;
	for (int i = 0; i < m_mol->natoms(); i++) {
		ac = m_mol->atom(i);
		if (ac->element != Atom::get_element_number("C")) continue;

		ao = 0;
		an = 0;
		ax = 0;
		m_mol->get_atoms_bonded_to(ac, aa);
		if (aa.size() > 3) continue;
		for (int j = 0; j < aa.size(); j++) {
			//if (aa[j]->element == Atom::get_element_number("O") && m_mol->bond_between(ac, aa[j]) == 2) {
			if (aa[j]->element == Atom::get_element_number("O")) {
				m_mol->get_atoms_bonded_to(aa[j], ab);
				if (ab.size() == 1 && ao == 0) {
					ao = aa[j];
				} else {
					ax = aa[j];
				}
			} else if (aa[j]->element == Atom::get_element_number("N") &&
				m_mol->bond_between(ac, aa[j]) >= 2) {
				an = aa[j];
				m_mol->add_bond(ac, an, 1); // make any C-N bond single
			} else {
				ax = aa[j];
			}
		}
		if (an != 0 && ao != 0 && ax != 0) {
			m_mol->add_bond(ac, ao, 2); // X-C(-N,=O)
			m_mol->add_bond(ac, ax, 1); // X-C(-N,=O)
		}
	}
}

//
void BondMaker::make_bonds_carboxyl()
{
	Atom* ac;
	vector<Atom*> aa;
	int omin;
	double d, domin;
	for (int i = 0; i < m_mol->natoms(); i++) {
		ac = m_mol->atom(i);
		if (ac->element != Atom::get_element_number("C")) continue;

		m_mol->get_atoms_bonded_to(ac, aa);
		if (aa.size() != 3) continue;
		for (vector<Atom*>::iterator p = aa.begin(); p != aa.end(); ) {
			if ((*p)->element != Atom::get_element_number("O")) {
				p = aa.erase(p);
			} else {
				++p;
			}
		}
		if (aa.size() < 2) continue;
		omin = 0;
		domin = aa[0]->pos.distance_from(ac->pos);
		for (int j = 1; j < aa.size(); j++) {
			d = aa[j]->pos.distance_from(ac->pos);
			if (d < domin) {
				omin = j;
				domin = d;
			}
		}
		if (m_mol->bond_between(ac, aa[omin]) != 2) continue;
		for (int j = 0; j < aa.size(); j++) {
			if (j == omin) continue;
			m_mol->add_bond(ac, aa[j], 1); // C-O
		}
	}
}

//
void BondMaker::make_bonds_phosphate()
{
	Atom* ap;
	vector<Atom*> ao;
	vector<Atom*> ax;

	int b, n;
	int omin;
	double domin;
	Bond* bp;
	for (int i = 0; i < m_mol->natoms(); i++) {
		ap = m_mol->atom(i);
		if (ap->element != Atom::get_element_number("P")) continue;

		m_mol->get_atoms_bonded_to(ap, ao);
		if (ao.size() != 4) continue;
		for (n = 0; n < ao.size(); n++) {
			if (ao[n]->element != Atom::get_element_number("O")) break;
		}
		if (n != ao.size()) continue;

		// initially, set all P-O bonds single
		for (n = 0; n < ao.size(); n++) {
			b = m_mol->bond_between(ap, ao[n]);
			if (b != 1) m_mol->add_bond(ap, ao[n], 1); // P-O
		}

		omin = -1;
		domin = 999.9;
		for (n = 0; n < ao.size(); n++) {
			m_mol->get_atoms_bonded_to(ao[n], ax);
			if (ax.size() != 1) continue;
			bp = get_bond(ap, ao[n]);
			if (bp->d < domin) {
				omin = n;
				domin = bp->d;
			}
		}
		if (omin != -1) {
			m_mol->add_bond(ap, ao[omin], 2); // P=O
		}
	}
}

//
void BondMaker::make_bonds_nitro()
{
	Atom* an;
	Atom* ao1;
	Atom* ao2;
	Atom* ax;
	vector<Atom*> aa;
	for (int i = 0; i < m_mol->natoms(); i++) {
		an = m_mol->atom(i);
		if (an->element != Atom::get_element_number("N")) continue;
		m_mol->get_atoms_bonded_to(an, aa);
		if (aa.size() != 3) continue;

		ao1 = 0;
		ao2 = 0;
		ax = 0;
		for (int j = 0; j < aa.size(); j++) {
			if (aa[j]->element == Atom::get_element_number("O")) {
				if (ao1) ao2 = aa[j];
				else ao1 = aa[j];
			} else {
				ax = aa[j];
			}
		}
		if (ao1 == 0 || ao2 == 0 || ax == 0) continue;
		m_mol->get_atoms_bonded_to(ao1, aa);
		if (aa.size() != 1) continue;
		m_mol->get_atoms_bonded_to(ao2, aa);
		if (aa.size() != 1) continue;

		m_mol->add_bond(an, ao1, 2); // X-N=O,=O
		m_mol->add_bond(an, ao2, 2);
		m_mol->add_bond(an, ax, 1);
	}
}

//
void BondMaker::make_bonds_ar()
{
	get_ar_rings();
	get_ar_nbonds_left();
	for (int i = 0; i < m_ar_rings.size(); i++) {
		make_bonds_ar(m_ar_rings[i]);
	}
}

//
void BondMaker::get_ar_rings()
{
	m_ar_rings.clear();
	vector<vector<Atom*> > rings;
	m_mol->get_rings(rings);
	for (int i = 0; i < rings.size(); i++) {
		if (rings[i].size() == 6 && is_ar6(rings[i]))
			m_ar_rings.push_back(rings[i]);
		else if (rings[i].size() == 5 && is_ar5(rings[i]))
			m_ar_rings.push_back(rings[i]);
	}

	// initially, set all aromatic bonds single
	for (int i = 0; i < m_ar_rings.size(); i++) {
		vector<Atom*>& ring = m_ar_rings[i];
		for (int j = 0; j < ring.size()-1; j++) {
			m_mol->add_bond(ring[j], ring[j+1], 1);
		}
		m_mol->add_bond(ring.back(), ring[0], 1);
	}
}

//
void BondMaker::get_ar_nbonds_left()
{
	m_ar_nbonds_left.clear();

	map<Atom*, int>::iterator p;
	Atom* a;
	vector<Atom*> as;
	for (int i = 0; i < m_ar_rings.size(); i++) {
		const vector<Atom*>& ring = m_ar_rings[i];
		for (int j = 0; j < ring.size(); j++) {
			a = ring[j];
			p = m_ar_nbonds_left.find(a);
			if (p != m_ar_nbonds_left.end()) continue;

			if (a->element == Atom::get_element_number("C")) {
				m_ar_nbonds_left[a] = 4;
			} else if (a->element == Atom::get_element_number("N")) {
				m_ar_nbonds_left[a] = 3;
			} else {
				m_ar_nbonds_left[a] = 2;
			}

			m_mol->get_atoms_bonded_to(a, as);
			for (int k = 0; k < as.size(); k++) {
				m_ar_nbonds_left[a] -= m_mol->bond_between(a, as[k]);
			}
			if (m_ar_nbonds_left[a] < 0) m_ar_nbonds_left[a] = 0;
		}
	}
}

//
void BondMaker::make_bonds_ar(vector<Atom*>& ring)
{
	vector<Bond> bonds;
	get_bonds_in_ar(ring, bonds);
	int i0 = set_ar_bond0(bonds); // get the shortest bond
	if (i0 < 0) return;
	set_ar_bonds(i0, bonds);
	refine_ar_bonds(bonds);
	for (int i = 0; i < bonds.size(); i++) {
		m_mol->add_bond(bonds[i].a1, bonds[i].a2, bonds[i].t);
	}
}

//
void BondMaker::get_bonds_in_ar(vector<Atom*>& ring, vector<Bond>& bonds)
{
	Bond b;
	bonds.clear();
	for (int i = 0; i < ring.size()-1; i++) {
		b.a1 = ring[i];
		b.a2 = ring[i+1];
		b.t = m_mol->bond_between(b.a1, b.a2);
		bonds.push_back(b);
	}
	b.a1 = ring.back();
	b.a2 = ring[0];
	b.t = m_mol->bond_between(b.a1, b.a2);
	bonds.push_back(b);
	// 'bonds' have bonds 0->1, 1->2, ..., n-1->0 in this order

	Bond* bp = 0;
	for (int i = 0; i < bonds.size(); i++) {
		if (bp = get_bond(bonds[i].a1, bonds[i].a2)) {
			bonds[i].d = bp->d;
		}
	}
}

//
int BondMaker::set_ar_bond0(vector<Bond>& bonds)
{
	if (bonds.empty()) return -1;
	int imin = -1;
	double dmin = 999.9;
	int i1, i2;
	for (int i = 0; i < bonds.size(); i++) {
		if (m_ar_nbonds_left[bonds[i].a1] <= 0 ||
			m_ar_nbonds_left[bonds[i].a2] <= 0) continue;
		i1 = i - 1;
		i2 = i + 1;
		if (i1 < 0) i1 = bonds.size()-1;
		if (i2 >= bonds.size()) i2 = 0;
		if (bonds[i].t == 2 || bonds[i1].t == 2 || bonds[i2].t == 2) continue;

		if (bonds[i].d < dmin) {
			imin = i;
			dmin = bonds[i].d;
		}
	}
	if (imin >= 0) {
		bonds[imin].t = 2;
		m_ar_nbonds_left[bonds[imin].a1]--;
		m_ar_nbonds_left[bonds[imin].a2]--;
	}

	return imin;
}

//
void BondMaker::set_ar_bonds(int i0, vector<Bond>& bonds)
{
	if (i0 < 0) return;

	int i1, i2;
	for (int i = i0-2; i >= 0; ) {
		if (m_ar_nbonds_left[bonds[i].a1] <= 0 || m_ar_nbonds_left[bonds[i].a2] <= 0) {
			i--;
			continue;
		}

		i1 = i - 1;
		i2 = i + 1;
		if (i1 < 0) i1 = bonds.size()-1;
		if (i2 >= bonds.size()) i2 = 0;
		if (bonds[i].t == 2 || bonds[i1].t == 2 || bonds[i2].t == 2) {
			i--;
			continue;
		}

		bonds[i].t = 2;
		m_ar_nbonds_left[bonds[i].a1]--;
		m_ar_nbonds_left[bonds[i].a2]--;
		i -= 2;
	}

	for (int i = i0+2; i < bonds.size(); ) {
		if (m_ar_nbonds_left[bonds[i].a1] <= 0 || m_ar_nbonds_left[bonds[i].a2] <= 0) {
			i++;
			continue;
		}

		i1 = i - 1;
		i2 = i + 1;
		if (i1 < 0) i1 = bonds.size()-1;
		if (i2 >= bonds.size()) i2 = 0;
		if (bonds[i].t == 2 || bonds[i1].t == 2 || bonds[i2].t == 2) {
			i++;
			continue;
		}

		bonds[i].t = 2;
		m_ar_nbonds_left[bonds[i].a1]--;
		m_ar_nbonds_left[bonds[i].a2]--;
		i += 2;
	}
}

//
void BondMaker::refine_ar_bonds(vector<Bond>& bonds)
{
	// for NAD case
	Atom* ac;
	Atom* an;
	int i1, i2;
	for (int i = 0; i < bonds.size(); i++) {
		if (bonds[i].t != 1) continue;
		i1 = i - 1;
		i2 = i + 1;
		if (i1 < 0) i1 = bonds.size()-1;
		if (i2 >= bonds.size()) i2 = 0;
		if (bonds[i1].t != 1 || bonds[i2].t != 1) continue;

		if (bonds[i].a1->element == Atom::get_element_number("C") &&
			bonds[i].a2->element == Atom::get_element_number("N")) {
			ac = bonds[i].a1;
			an = bonds[i].a2;
		} else if (bonds[i].a1->element == Atom::get_element_number("N") &&
			bonds[i].a2->element == Atom::get_element_number("C")) {
			ac = bonds[i].a2;
			an = bonds[i].a1;
		} else {
			continue;
		}

		if (m_ar_nbonds_left[ac] != 2) continue;
		bonds[i].t = 2;
		m_ar_nbonds_left[ac]--;
		// m_ar_nbonds_left[an] is not changed, which causes N->N+
	}
}

//
bool BondMaker::is_ar(vector<Atom*>& ring, int n, double dar)
{
	Atom* a1;
	Atom* a2;
	double d, dd1, dd2, dd_ar;
	for (int i = 0; i < n-1; i++) {
		a1 = ring[i];
		if (a1->element != Atom::get_element_number("C")) continue;
		a2 = ring[i+1];
		if (a2->element != Atom::get_element_number("C")) continue;
		d = a1->pos.distance_from(a2->pos);
		if ((dd1 = C_C_1 - d) < 0) continue;
		if ((dd2 = d - C_C_2) < 0) continue;
		dd_ar = std::fabs(d - dar);
		if (dd_ar < dd1) return true;
	}
	a1 = ring[0];
	if (a1->element != Atom::get_element_number("C")) return false;
	a2 = ring[n-1];
	if (a2->element != Atom::get_element_number("C")) return false;
	d = a1->pos.distance_from(a2->pos);
	if ((dd1 = C_C_1 - d) < 0) return false;
	if ((dd2 = d - C_C_2) < 0) return false;
	dd_ar = std::fabs(d - dar);
	if (dd_ar < dd1) return true;
	return false;
}
