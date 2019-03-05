#include <vector>
#include <string>
#include <map>
//#include <list>
#include <set>
#include <iostream>
//#include <fstream>
#include <sstream>
//#include <algorithm>
#include <cctype>
#include <cmath>
//#include <valarray>
#include "matrix.h"
#include "atom.h"
#include "mol_fin.h"
#include "bit_string.h"

using namespace Cba;
using std::vector;
using std::string;
using std::map;
using std::set;
using std::ostream;
using std::istream;
using std::endl;

const int MolFin::DISTANCE_RANGE_MIN = 1;
const int MolFin::DISTANCE_RANGE_MAX = 6;

const double MolFin::P_VALUES_TC_STEP = 0.01;
const double MolFin::p_values[] = {
1.0, 0.896783, 0.845635, 0.801928, 0.76169, 0.723617, 0.686598, 0.65102, 0.61733, 0.583702,
0.55139, 0.520139, 0.489891, 0.460459, 0.43206, 0.405317, 0.379538, 0.354812, 0.331312, 0.30866,
0.286996, 0.267175, 0.24826, 0.23068, 0.2136, 0.197833, 0.182746, 0.168658, 0.155541, 0.142981,
0.13114, 0.120386, 0.110363, 0.100963, 0.0920558, 0.0838314, 0.0761871, 0.0691881, 0.0625889, 0.0566456,
0.051171, 0.0461083, 0.0415072, 0.0369433, 0.0330176, 0.0293661, 0.0260647, 0.0231046, 0.0204592, 0.0180808,
0.0159589, 0.0140385, 0.0123374, 0.0107778, 0.00943563, 0.00829862, 0.007128, 0.0061466, 0.00526952, 0.00456751,
0.00390086, 0.00334031, 0.00284165, 0.00243494, 0.00209013, 0.00176299, 0.00146592, 0.00124311, 0.00103799, 0.000838173,
0.000696709, 0.000595916, 0.000503965, 0.000413781, 0.000334208, 0.000267013, 0.000205122, 0.000173293, 0.000130854, 9.90246e-05,
8.13416e-05, 6.5427e-05, 5.12806e-05, 4.06708e-05, 3.18293e-05, 2.47562e-05, 1.7683e-05, 1.41464e-05, 1.23781e-05, 8.84148e-06,
5.30489e-06, 5.30489e-06, 1.7683e-06, 1.7683e-06, 1.7683e-06, 1.7683e-06, 0.0, 0.0, 0.0, 0.0
};

//#include "mol_fin_weights.h" // weight data
//std::valarray<double> MolFin::m_weights_sqr;

//
MolFin::MolFin() : m_natoms(0), m_csqr(0), m_debug(true)
{
	if (AtomTriad::bit_pos_tbl.size() == 0) {
		AtomTriad::set_bit_pos_tbl();
	}
	m_fin.resize(AtomTriad::bit_pos_tbl.size());
	for (int i = 0; i < AtomTriad::bit_pos_tbl.size(); i++) m_fin[i] = 0;
	/*
	if (m_weights_sqr.size() != AtomTriad::bit_pos_tbl.size()) {
		m_weights_sqr.resize(AtomTriad::bit_pos_tbl.size());
		for (int i = 0; i < AtomTriad::bit_pos_tbl.size(); i++)
			m_weights_sqr[i] = m_weights[i]*m_weights[i];
	}
	*/
}

//
MolFin::~MolFin() {}

//
void MolFin::get_triad(size_t i, int& c1, int& c2, int& c3, int& d12, int& d13, int& d23)
{
	MolFin mf; // to confirm initialization
	c1 = AtomTriad::triad_at_pos.at(i).c1;
	c2 = AtomTriad::triad_at_pos.at(i).c2;
	c3 = AtomTriad::triad_at_pos.at(i).c3;
	d12 = AtomTriad::triad_at_pos.at(i).d12;
	d13 = AtomTriad::triad_at_pos.at(i).d13;
	d23 = AtomTriad::triad_at_pos.at(i).d23;
}

//
double MolFin::dsqr(const MolFin& mf) const
{
	if (size() != mf.size()) return -1.0;
	double dsq = 0.0;
	for (int i = 0; i < size(); i++) {
		dsq += (m_fin[i]-mf.m_fin[i])*(m_fin[i]-mf.m_fin[i]);
	}
	return dsq;
}

//
double MolFin::rmsd(const MolFin& mf) const
{
	double dsq = dsqr(mf);
	if (dsq < 0.0) return dsq;
	return std::sqrt(dsq/size());
}

//
double MolFin::dist(const MolFin& mf) const
{
	double dsq = dsqr(mf);
	if (dsq < 0.0) return dsq;
	return std::sqrt(dsq);
}

//
double MolFin::tc(const MolFin& mf) const
{
	int c = 0;
	int n = size();
	for (int i = 0; i < n; i++) c += m_fin[i]*mf.m_fin[i];
	double s = (double)(m_csqr + mf.m_csqr - c);
	if (s <= 0.0) return 0.0;
	return c/s;
	//return std::sqrt(c/s);
}

//
double MolFin::pvalue(double tc)
{
	int n = (int)(tc/P_VALUES_TC_STEP);
	if (n >= 100) return 0.0;

	double p1, p2;
	p1 = p_values[n];
	if (n == 99) p2 = 0.0;
	else p2 = p_values[n+1];

	double d = tc - n*P_VALUES_TC_STEP;
	return (p1 + d*((p2-p1)/P_VALUES_TC_STEP));
}

//
double MolFin::tc2(const MolFin& mf) const
{
	int a = 0;
	int b = 0;
	int c = 0;
	int n = size();
	for (int i = 0; i < n; i++) {
		a += m_fin[i];
		b += mf.m_fin[i];
		c += std::min(m_fin[i], mf.m_fin[i]);
	}
	double s = (double)(a + b - c);
	if (s <= 0.0) return 0.0;
	return std::sqrt(c/s);
}

//
/*
double MolFin::tcw(const MolFin& mf) const
{
	double a = 0;
	double b = 0;
	double c = 0;
	int n = size();
	for (int i = 0; i < n; i++) {
		a += m_weights_sqr[i]*m_fin[i]*m_fin[i];
		b += m_weights_sqr[i]*mf.m_fin[i]*mf.m_fin[i];
		c += m_weights_sqr[i]*m_fin[i]*mf.m_fin[i];
	}
	double s = a + b - c;
	if (s <= 0.0) return 0.0;
	return std::sqrt(c/s);
}
*/

//
double MolFin::tc_bin(const MolFin& mf) const
{
	int c = 0;
	int n = size();
	int a = 0; // A and B
	int b = 0; // A or B
	for (int i = 0; i < n; i++) {
		if (m_fin[i] && mf.m_fin[i]) a++;
		if (m_fin[i] || mf.m_fin[i]) b++;
	}
	if (! b) return 0.0;
	return (double)a/(double)b;
}

//
void MolFin::make(Molecule& mol)
{
	mol.assign_atom_types();
	for (int i = 0; i < size(); i++) m_fin[i] = 0;
	m_natoms = mol.natoms();
	m_id = mol.id();
	int ci, cj, ck;
	int dij, dik, djk;

	// for all atom triads
	for (int i = 0; i < m_natoms-2; i++) {
		if ((ci = get_atom_class(mol.atom(i))) == 0) continue;
		for (int j = i+1; j < m_natoms-1; j++) {
			if ((cj = get_atom_class(mol.atom(j))) == 0) continue;
			if ((dij = get_atom_dij(mol,i,j)) == 0) continue;
			for (int k = j+1; k < m_natoms; k++) {
				if ((ck = get_atom_class(mol.atom(k))) == 0) continue;
				if ((dik = get_atom_dij(mol,i,k)) == 0) continue;
				if ((djk = get_atom_dij(mol,j,k)) == 0) continue;

				m_atom_triad.set(ci,cj,ck,dij,dik,djk);
				int p = m_atom_triad.bit_pos();
				if (0 <= p && p < size()) m_fin[p]++;
			}
		}
	}
	calc_csqr();
}

//
void MolFin::calc_csqr()
{
	m_csqr = 0;
	for (int i = 0; i < size(); i++) m_csqr += m_fin[i]*m_fin[i];
}

//
int MolFin::get_atom_class(const Atom* a)
{
	int c = a->pc_class;
	/*
	if (c == Atom::CATION) return Atom::DONOR;
	if (c == Atom::ANION) return Atom::ACCEPTOR;
	if (c == Atom::DONOR) return Atom::DONOR;
	if (c == Atom::ACCEPTOR) return Atom::ACCEPTOR;
	if (c == Atom::POLAR) return Atom::POLAR;
	if (c == Atom::HYDROPHOBIC) return Atom::HYDROPHOBIC;
	if (c == Atom::NONE) return Atom::HYDROPHOBIC;
	*/
	/*
	if (c == Atom::CATION) return Atom::DONOR;
	if (c == Atom::ANION) return Atom::ACCEPTOR;
	if (c == Atom::DONOR) return Atom::DONOR;
	if (c == Atom::ACCEPTOR) return Atom::ACCEPTOR;
	if (c == Atom::POLAR) return Atom::POLAR;
	if (c == Atom::HYDROPHOBIC) return Atom::HYDROPHOBIC;
	if (c == Atom::NONE) return 0;
	*/
	/*
	if (c == Atom::CATION) return Atom::CATION;
	if (c == Atom::ANION) return Atom::ANION;
	if (c == Atom::DONOR) return Atom::DONOR;
	if (c == Atom::ACCEPTOR) return Atom::ACCEPTOR;
	if (c == Atom::POLAR) return Atom::POLAR;
	if (c == Atom::HYDROPHOBIC) return Atom::HYDROPHOBIC;
	if (c == Atom::NONE) return Atom::HYDROPHOBIC;
	*/
	/*
	if (c == Atom::CATION) return Atom::CATION;
	if (c == Atom::ANION) return Atom::ANION;
	if (c == Atom::DONOR) return Atom::DONOR;
	if (c == Atom::ACCEPTOR) return Atom::ACCEPTOR;
	if (c == Atom::POLAR) return Atom::POLAR;
	if (c == Atom::HYDROPHOBIC) return Atom::HYDROPHOBIC;
	if (c == Atom::NONE) return 0;
	*/
	/*
	if (c == Atom::CATION) return Atom::CATION;
	if (c == Atom::ANION) return Atom::ANION;
	if (c == Atom::DONOR) return Atom::DONOR;
	if (c == Atom::ACCEPTOR) return Atom::ACCEPTOR;
	if (c == Atom::POLAR) return Atom::POLAR;
	if (c == Atom::HYDROPHOBIC) return Atom::HYDROPHOBIC;
	if (c == Atom::NONE) return Atom::NONE;
	*/
	/*
	if (c == Atom::CATION) return Atom::DONOR;
	if (c == Atom::ANION) return Atom::ACCEPTOR;
	if (c == Atom::DONOR) return Atom::DONOR;
	if (c == Atom::ACCEPTOR) return Atom::ACCEPTOR;
	if (c == Atom::POLAR) return Atom::POLAR;
	if (c == Atom::HYDROPHOBIC) return Atom::HYDROPHOBIC;
	if (c == Atom::NONE) return Atom::NONE;
	*/
	if (c == Atom::CATION) return DONOR;
	if (c == Atom::ANION) return ACCEPTOR;
	if (c == Atom::DONOR) return DONOR;
	if (c == Atom::ACCEPTOR) return ACCEPTOR;
	if (c == Atom::POLAR) return POLAR;
	if (c == Atom::HYDROPHOBIC) {
		if (a->properties.get_bit(Atom::ar) == 1) return HYDROPHOBIC_AROMATIC;
		return HYDROPHOBIC_NON_AROMATIC;
	}
	if (c == Atom::NONE) return NONE;
	return 0;
}

//
int MolFin::get_atom_dij(Molecule& mol, int i , int j)
{
	int d = mol.get_graph_distance(i, j);
	if (1 <= d && d <= 4) return d;
	if (5 <= d && d <= 9) return 5;
	if (d >= 10) return 6;
	return 0;
}

//
void MolFin::clear()
{
	m_id.erase();
	m_natoms = 0;
	m_csqr = 0;
	for (int i = 0; i < size(); i++) m_fin[i] = 0;
}

//
void MolFin::write(ostream& os) const
{
	map<int, vector<int> > counts; // count -> serial numbers of elements with that count
	for (int i = 0; i < size(); i++) {
		if (m_fin[i] > 0) counts[m_fin[i]].push_back(i);
	}

	os << m_id << endl;
	os << m_natoms << endl;
	os << m_csqr << endl;
	for (map<int, vector<int> >::const_iterator mci = counts.begin(); mci != counts.end(); ++mci) {
		os << mci->first;
		const vector<int>& elems = mci->second;
		for (vector<int>::const_iterator vci = elems.begin(); vci != elems.end(); ++vci) {
			os << ' ' << *vci;
		}
		os << endl;
	}
	os << '/' << endl;

	/*
	int nc = 50;
	for (int i = 0; i < size(); i++) {
		os << m_fin[i];
		if ((i+1)%nc == 0) os << endl;
		else os << ' ';
	}
	os << -1 << endl;
	*/
}

//
bool MolFin::read(istream& is)
{
	clear();

	if (! std::getline(is, m_id)) return false;

	string s;
	if (std::getline(is, s)) {
		m_natoms = std::atoi(s.c_str());
	} else {
		return false;
	}
	if (std::getline(is, s)) {
		m_csqr = std::atoi(s.c_str());
	} else {
		return false;
	}


	int nc, ne;
	while (std::getline(is, s)) {
		if (s.at(0) == '/') break;
		std::istringstream iss(s);
		if (! (iss >> nc)) return false;
		while (iss >> ne) {
			m_fin[ne] = nc;
		}
	}
	if (s.at(0) != '/') return false;
	//calc_csqr();

	/*
	int k, c;
	for (k = 0; is >> c; k++) {
		if (c < 0) break;
		m_fin[k] = c;
	}
	if (k != size()) {
		std::cerr << m_id << ": invalid fingerprint size: " << k << endl;
		return false;
	}
	*/
	return true;
}

//
int MolFin::max_count() const
{
	int mc = 0;
	for (int i = 0; i < size(); i++) {
		if (m_fin[i] > mc) mc = m_fin[i];
	}
	return mc;
}

///////////////////////////////////
// definitions of AtomTriad follow
///////////////////////////////////

//
map<MolFin::AtomTriad, int, MolFin::AtomTriadLess> MolFin::AtomTriad::bit_pos_tbl;
vector<MolFin::AtomTriad> MolFin::AtomTriad::triad_at_pos;

//
MolFin::AtomTriad::AtomTriad(int a1, int a2, int a3, int t12, int t13, int t23)
{
	set(a1, a2, a3, t12, t13, t23);
}

//
void MolFin::AtomTriad::set_as_is(int a1, int a2, int a3, int t12, int t13, int t23)
{
	c1 = a1;
	c2 = a2;
	c3 = a3;
	d12 = t12;
	d13 = t13;
	d23 = t23;
}

//
void MolFin::AtomTriad::set(int a1, int a2, int a3, int t12, int t13, int t23)
{
	clear();
	if (! satisfy_triangle_inequality(t12,t13,t23)) return;

	// c1 <= c2 <= c3
	// when all atom classes are identical;
	//  d12 <= d13 <= d23
	// when two of atom classes are identical;
	//  for c1 = c2 < c3, d13 <= d23
	//  for c1 < c2 = c3 , d12 <= d13

	// when two of atom classes are identical;
	//  for c1 = c2, d13 <= d23
	//  for c1 = c3, d12 <= d23
	//  for c2 = c3, d12 <= d13

	if (a1 < a2 && a2 < a3) { // when a1 < a2 < a3
		set_as_is(a1,a2,a3,t12,t13,t23);
	} else if (a1 < a3 && a3 < a2) { // when a1 < a3 < a2
		set_as_is(a1,a3,a2,t13,t12,t23);
	} else if (a2 < a1 && a1 < a3) { // when a2 < a1 < a3
		set_as_is(a2,a1,a3,t12,t23,t13);
	} else if (a2 < a3 && a3 < a1) { // when a2 < a3 < a1
		set_as_is(a2,a3,a1,t23,t12,t13);
	} else if (a3 < a1 && a1 < a2) { // when a3 < a1 < a2
		set_as_is(a3,a1,a2,t13,t23,t12);
	} else if (a3 < a2 && a2 < a1) { // when a3 < a2 < a1
		set_as_is(a3,a2,a1,t23,t13,t12);
	} else if (a1 == a2 && a2 < a3) { // when a1 == a2 < a3
		if (t13 <= t23) { set_as_is(a1,a2,a3,t12,t13,t23); }
		else { set_as_is(a2,a1,a3,t12,t23,t13); }
	} else if (a3 < a1 && a1 == a2) { // when a3 < a1 == a2
		if (t13 <= t23) { set_as_is(a3,a1,a2,t13,t23,t12); }
		else { set_as_is(a3,a2,a1,t23,t13,t12); }
	} else if (a1 == a3 && a3 < a2) { // when a1 == a3 < a2
		if (t12 <= t23) { set_as_is(a1,a3,a2,t13,t12,t23); }
		else { set_as_is(a3,a1,a2,t13,t23,t12); }
	} else if (a2 < a1 && a1 == a3) { // when a2 < a1 == a3
		if (t12 <= t23) { set_as_is(a2,a1,a3,t12,t23,t13); }
		else { set_as_is(a2,a3,a1,t23,t12,t13); }
	} else if (a1 < a2 && a2 == a3) { // when a1 < a2 == a3
		if (t12 <= t13) { set_as_is(a1,a2,a3,t12,t13,t23); }
		else { set_as_is(a1,a3,a2,t13,t12,t23); }
	} else if (a2 == a3 && a3 < a1) { // when a2 == a3 < a1
		if (t12 <= t13) { set_as_is(a2,a3,a1,t23,t12,t13); }
		else { set_as_is(a3,a2,a1,t23,t13,t12); }
	} else { // a1 == a2 == a3
		if (t12 < t13 && t13 < t23) set_as_is(a1,a2,a3,t12,t13,t23);
		else if (t12 < t23 && t23 < t13) set_as_is(a2,a1,a3,t12,t23,t13);
		else if (t13 < t12 && t12 < t23) set_as_is(a1,a3,a2,t13,t12,t23);
		else if (t13 < t23 && t23 < t12) set_as_is(a3,a1,a2,t13,t23,t12);
		else if (t23 < t12 && t12 < t13) set_as_is(a2,a3,a1,t23,t12,t13);
		else if (t23 < t13 && t13 < t12) set_as_is(a3,a2,a1,t23,t13,t12);
		else if (t12 == t13 && t13 < t23) set_as_is(a1,a2,a3,t12,t13,t23);
		else if (t23 < t12 && t12 == t13) set_as_is(a2,a3,a1,t23,t12,t13);
		else if (t12 == t23 && t23 < t13) set_as_is(a2,a1,a3,t12,t23,t13);
		else if (t13 < t12 && t12 == t23) set_as_is(a1,a3,a2,t13,t12,t23);
		else if (t12 < t13 && t13 == t23) set_as_is(a1,a2,a3,t12,t13,t23);
		else if (t13 == t23 && t23 < t12) set_as_is(a3,a1,a2,t13,t23,t12);
		else set_as_is(a1,a2,a3,t12,t13,t23);
	}
}

//
int MolFin::AtomTriad::bit_pos() const
{
	map<AtomTriad, int, AtomTriadLess>::const_iterator mci = AtomTriad::bit_pos_tbl.find(*this);
	if (mci == AtomTriad::bit_pos_tbl.end()) {
		std::cerr << "no bit position found for ";
		this->print(std::cerr);
		return -1;
	}
	return mci->second;
}

//
bool MolFin::AtomTriadLess::operator()(
	const AtomTriad& t1, const AtomTriad& t2) const
{
	if (t1.c1 < t2.c1) return true;
	if (t1.c1 == t2.c1) {
		if (t1.c2 < t2.c2) return true;
		if (t1.c2 == t2.c2) {
			if (t1.c3 < t2.c3) return true;
			if (t1.c3 == t2.c3) {
				if (t1.d12 < t2.d12) return true;
				if (t1.d12 == t2.d12) {
					if (t1.d13 < t2.d13) return true;
					if (t1.d13 == t2.d13) {
						if (t1.d23 < t2.d23) return true;
					}
				}
			}
		}
	}
	return false;
}

//
void MolFin::AtomTriad::set_bit_pos_tbl()
{
	if (! bit_pos_tbl.empty()) return;
	bit_pos_tbl.clear();
	triad_at_pos.clear();
	AtomTriad a;
	int n = 0;

	// c1 <= c2 <= c3
	// when all atom classes are identical;
	//  d12 <= d13 <= d23
	// when two of atom classes are identical;
	//  for c1 = c2 < c3, d13 <= d23
	//  for c1 < c2 = c3 , d12 <= d13
	for (int c1 = 1; c1 < NO_OF_ATOM_TYPES; c1++) {
		for (int c2 = c1; c2 < NO_OF_ATOM_TYPES; c2++) {
			for (int c3 = c2; c3 < NO_OF_ATOM_TYPES; c3++) {
				for (int d12 = DISTANCE_RANGE_MIN; d12 <= DISTANCE_RANGE_MAX; d12++) {
					for (int d13 = DISTANCE_RANGE_MIN; d13 <= DISTANCE_RANGE_MAX; d13++) {
						for (int d23 = DISTANCE_RANGE_MIN; d23 <= DISTANCE_RANGE_MAX; d23++) {
							if (c1 == c2 && c2 == c3) {
								if (d12 <= d13 && d13 <= d23 && satisfy_triangle_inequality(d12,d13,d23)) {
									a.set_as_is(c1,c2,c3,d12,d13,d23);
									bit_pos_tbl[a] = n++;
									triad_at_pos.push_back(a);
								}
							} else if (c1 == c2) {
								if (d13 <= d23 && satisfy_triangle_inequality(d12,d13,d23)) {
									a.set_as_is(c1,c2,c3,d12,d13,d23);
									bit_pos_tbl[a] = n++;
									triad_at_pos.push_back(a);
								}
							} else if (c2 == c3) {
								if (d12 <= d13 && satisfy_triangle_inequality(d12,d13,d23)) {
									a.set_as_is(c1,c2,c3,d12,d13,d23);
									bit_pos_tbl[a] = n++;
									triad_at_pos.push_back(a);
								}
							} else {
								if (satisfy_triangle_inequality(d12,d13,d23)) {
									a.set_as_is(c1,c2,c3,d12,d13,d23);
									bit_pos_tbl[a] = n++;
									triad_at_pos.push_back(a);
								}
							}
						}
					}
				}
			}
		}
	}
}

//
bool MolFin::AtomTriad::satisfy_triangle_inequality(int t12, int t13, int t23)
{
	int dmax, d1, d2;
	if (t12 >= t13 && t12 >= t23) {
		dmax = t12;
		d1 = t13;
		d2 = t23;
	} else if (t13 >= t12 && t13 >= t23) {
		dmax = t13;
		d1 = t12;
		d2 = t23;
	} else if (t23 >= t12 && t23 >= t13) {
		dmax = t23;
		d1 = t12;
		d2 = t13;
	}

	if (dmax == 6) {
		dmax = 10;
		if (d1 == 5) d1 = 9;
		else if (d1 == 6) d1 = 10;
		if (d2 == 5) d2 = 9;
		else if (d2 == 6) d2 = 10;
		if ((d1+d2) >= dmax) return true;
		else return false;
	}
	if ((d1+d2) >= dmax) return true;
	return false;
}

//
void MolFin::AtomTriad::print(ostream& os) const
{
	os << c1 << ' ' << c2 << ' ' << c3 << ' '
		<< d12 << ' ' << d13 << ' ' << d23 << std::endl;
}

//
void MolFin::AtomTriad::clear()
{
	c1 = c2 = c3 = 0;
	d12 = d13 = d23 = 0;
}
