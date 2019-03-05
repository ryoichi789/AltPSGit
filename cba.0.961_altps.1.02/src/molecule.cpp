#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include <utility>
#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "molecule.h"
#include "atom.h"
#include "matrix.h"
#include "asa_calculator.h"
#include "bit_string.h"

using namespace Cba;
using std::vector;
using std::list;
using std::map;
using std::set;
using std::string;
using std::queue;

//
Molecule::Molecule() : m_type(Molecule::UNDEFINED), m_atom_types_assigned(false)
{}

//
Molecule::~Molecule()
{
	clear();
}

// clear:
void Molecule::clear()
{
	m_id.erase();
	m_name.erase();
	m_type = Molecule::UNDEFINED;
	m_atoms.clear();
	m_bonds.clear();
	for (vector<Atom*>::iterator
		p = m_all_atoms.begin(); p != m_all_atoms.end(); ++p) {
		delete *p;
	}
	m_all_atoms.clear();
	for (vector<Bond*>::iterator
		p = m_all_bonds.begin(); p != m_all_bonds.end(); ++p) {
		delete *p;
	}
	m_all_bonds.clear();
	clear_bonds_to_atom();
	m_rings.clear();
	m_rotatable_bonds.clear();
	m_fragments.clear();
	m_graph_distances.clear();
	m_atom_types_assigned = false;
}

//
Molecule::Molecule(const Molecule& mol)
{
	mol.clone(this);
}

//
Molecule& Molecule::operator=(const Molecule& mol)
{
	if (this != &mol) {
		clear();
		mol.clone(this);
	}
	return *this;
}

// bond_between
int Molecule::bond_between(int a1, int a2) const
{
	return bond_between(atom(a1), atom(a2));
}

// bond_between
int Molecule::bond_between(const Atom* a1, const Atom* a2) const
{
	for (vector<Bond*>::const_iterator
		p = m_bonds.begin(); p != m_bonds.end(); ++p)
		if (((*p)->atom1 == a1 && (*p)->atom2 == a2) ||
			((*p)->atom1 == a2 && (*p)->atom2 == a1))
			return (*p)->v;
	return 0;
}

// get_atoms_bonded_to:
void Molecule::get_atoms_bonded_to(const Atom* a, vector<Atom*>& as)
{
	as.clear();
	for (vector<Bond*>::iterator
		p = m_bonds.begin(); p != m_bonds.end(); ++p)
		if ((*p)->atom1 == a)
			as.push_back((*p)->atom2);
		else if ((*p)->atom2 == a)
			as.push_back((*p)->atom1);
}

//
void Molecule::get_atoms_bonded_to(size_t i, vector<int>& as)
{
	const Atom* a = atom(i);
	vector<Atom*> atoms;
	get_atoms_bonded_to(a, atoms);
	as.clear();
	for (int i = 0; i < atoms.size(); i++) {
		as.push_back(atom_number_of(atoms[i]));
	}
	std::sort(as.begin(), as.end());
}

// clone:
Molecule* Molecule::clone(Molecule* mol) const
{
	if (mol == this) return 0;
	if (mol) mol->clear();
	else mol = new Molecule;

	mol->set_id(id());
	mol->set_name(name());
	mol->set_type(type());

	for (int i = 0; i < natoms(); i++)
		mol->add_atom(*atom(i));

	int a1, a2;
	for (int i = 0; i < nbonds(); i++) {
		a1 = atom_number_of(m_bonds[i]->atom1);
		a2 = atom_number_of(m_bonds[i]->atom2);
		mol->add_bond(mol->atom(a1), mol->atom(a2), m_bonds[i]->v);
	}

	mol->m_rings = m_rings;
	mol->m_rotatable_bonds = m_rotatable_bonds;
	mol->m_fragments = m_fragments;
	mol->m_atom_types_assigned = m_atom_types_assigned;
	mol->m_graph_distances = m_graph_distances;

	return mol;
}

//
bool Molecule::clone(Molecule& mol) const
{
	if (clone(&mol) == 0) return false;
	return true;
}

// shift_by:
void Molecule::shift_by(double x, double y, double z)
{
	for (int i = 0; i < m_all_atoms.size(); i++)
		m_all_atoms[i]->pos.shift_by(x, y, z);
}

// shift_by:
void Molecule::shift_by(const Position& p)
{
	for (int i = 0; i < m_all_atoms.size(); i++)
		m_all_atoms[i]->pos.shift_by(p);
}

// rotate_by:
void Molecule::rotate_by(const Matrix& rot)
{
	for (int i = 0; i < m_all_atoms.size(); i++)
		m_all_atoms[i]->pos.rotate_by(rot);
}

// add_atom:
void Molecule::add_atom(const Atom& a)
{
	Atom* atom = new Atom;
	*atom = a;
	m_all_atoms.push_back(atom);
	m_atoms.push_back(atom);
}

// add_bond:
void Molecule::add_bond(Atom* a1, Atom* a2, int v)
{
	if (a1 == a2) return;
	if (std::find(m_atoms.begin(), m_atoms.end(), a1) == m_atoms.end() ||
		std::find(m_atoms.begin(), m_atoms.end(), a2) == m_atoms.end())
		return;

	for (vector<Bond*>::iterator p = m_all_bonds.begin(); p != m_all_bonds.end(); ++p) {
		if (((*p)->atom1 == a1 && (*p)->atom2 == a2) ||
			((*p)->atom1 == a2 && (*p)->atom2 == a1)) {
			m_all_bonds.erase(p);
			break;
		}
	}
	for (vector<Bond*>::iterator p = m_bonds.begin(); p != m_bonds.end(); ++p) {
		if (((*p)->atom1 == a1 && (*p)->atom2 == a2) ||
			((*p)->atom1 == a2 && (*p)->atom2 == a1)) {
			delete *p;
			m_bonds.erase(p);
			break;
		}
	}

	Bond* b = new Bond;
	b->atom1 = a1;
	b->atom2 = a2;
	b->v = v;
	m_all_bonds.push_back(b);
	m_bonds.push_back(b);
}

// remove_hydrogens:
void Molecule::remove_hydrogens()
{
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ) {
		if ((*p)->element == Atom::get_element_number("H")) {
			remove_bonds_to(*p);
			p = m_atoms.erase(p);
		} else {
			++p;
		}
	}

	m_rings.clear();
	m_rotatable_bonds.clear();
	m_fragments.clear();
}

//
void Molecule::remove_atom(const Atom* a)
{
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		if (*p == a) {
			remove_bonds_to(*p);
			m_atoms.erase(p);
			break;
		}
	}

	m_rings.clear();
	m_rotatable_bonds.clear();
	m_fragments.clear();
}

// remove_bonds_to:
void Molecule::remove_bonds_to(const Atom* a)
{
	for (vector<Bond*>::iterator p = m_bonds.begin(); p != m_bonds.end(); ) {
		if ((*p)->atom1 == a || (*p)->atom2 == a) p = m_bonds.erase(p);
		else ++p;
	}
}

// atom_number_of:
int Molecule::atom_number_of(const Atom* a) const
{
	if (!a) return -1;
	for (int i = 0; i < natoms(); i++)
		if (a == atom(i)) return i;
	return -1;
}

//
int Molecule::bond_number_of(const Bond* b) const
{
	if (!b) return -1;
	for (int i = 0; i < nbonds(); i++)
		if (b == m_bonds[i]) return i;
	return -1;
}

// assign_atom_types:
void Molecule::assign_atom_types()
{
	if (m_atom_types_assigned) return;
	remove_hydrogens();
	define_bonds_to_atom();
	assign_atom_properties();
	assign_atom_pc_classes();
	clear_bonds_to_atom();
	m_atom_types_assigned = true;
}

// define_bonds_to_atom:
void Molecule::define_bonds_to_atom()
{
	clear_bonds_to_atom();
	Bond* b;
	for (vector<Bond*>::iterator
		p = m_bonds.begin(); p != m_bonds.end(); ++p) {
		b = *p;
		m_bonds_to_atom[b->atom1].push_back(b);
		m_bonds_to_atom[b->atom2].push_back(b);
	}
}

// clear_bonds_to_atom:
void Molecule::clear_bonds_to_atom()
{
	m_bonds_to_atom.clear();
}

// get_atoms_bonded_to_atom:
void Molecule::get_atoms_bonded_to_atom(Atom* a, vector<Atom*>& as)
{
	as.clear();
	map<Atom*,vector<Bond*> >::iterator p = m_bonds_to_atom.find(a);
	if (p != m_bonds_to_atom.end()) {
		vector<Bond*>& bs = p->second;
		for (vector<Bond*>::iterator q = bs.begin(); q != bs.end(); ++q) {
			if ((*q)->atom1 == a) as.push_back((*q)->atom2);
			else  as.push_back((*q)->atom1);
		}
	}
}

// assign_atom_properties:
void Molecule::assign_atom_properties()
{
	clear_atom_properties();
	assign_atom_sp_sp2_res();
	assign_atom_conj();
	assign_atom_sp3_res();
	assign_atom_nneighbors();
	assign_atom_aromaticity();
	assign_atom_electronegativity();
	assign_atom_amide();
	assign_atom_carboxylate();
}

//clear_atom_properties:
void Molecule::clear_atom_properties()
{
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p)
		(*p)->properties.clear();
}

// assign_atom_sp_sp2_res:
void Molecule::assign_atom_sp_sp2_res()
{
	const Bond* b;
	map<Atom*, size_t> n_double_bonds;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p)
		n_double_bonds[*p] = 0;
	for (vector<Bond*>::const_iterator
		p = m_bonds.begin(); p != m_bonds.end(); ++p) {
		b = *p;
		if (b->v == 2) {
			b->atom1->properties.set_bit(Atom::sp2, 1);
			b->atom2->properties.set_bit(Atom::sp2, 1);
			b->atom1->properties.set_bit(Atom::res, 1);
			b->atom2->properties.set_bit(Atom::res, 1);
			n_double_bonds[b->atom1]++;
			n_double_bonds[b->atom2]++;
		} else if (b->v == 3) {
			b->atom1->properties.set_bit(Atom::sp, 1);
			b->atom2->properties.set_bit(Atom::sp, 1);
			b->atom1->properties.set_bit(Atom::res, 1);
			b->atom2->properties.set_bit(Atom::res, 1);
		}
	}
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		if (n_double_bonds[*p] >= 2) {
			(*p)->properties.set_bit(Atom::sp, 1);
			(*p)->properties.set_bit(Atom::res, 1);
		}
	}
}

// assign_atom_conj:
void Molecule::assign_atom_conj()
{
	// conj any atom single-bonded to sp2 or sp
	const Bond* b;
	for (vector<Bond*>::const_iterator
		p = m_bonds.begin(); p != m_bonds.end(); ++p) {
		b = *p;
		if (b->v == 1) {
			if (b->atom1->properties.get_bit(Atom::sp) ||
				b->atom1->properties.get_bit(Atom::sp2))
				b->atom2->properties.set_bit(Atom::conj, 1);
			if (b->atom2->properties.get_bit(Atom::sp) ||
				b->atom2->properties.get_bit(Atom::sp2))
				b->atom1->properties.set_bit(Atom::conj, 1);
		}
	}
}

// assign_atom_sp3_res:
void Molecule::assign_atom_sp3_res()
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a = *p;
		if ((! a->properties.get_bit(Atom::sp)) &&
			(! a->properties.get_bit(Atom::sp2))) {
			a->properties.set_bit(Atom::sp3, 1);
			if (a->properties.get_bit(Atom::conj) &&
				a->element == Atom::get_element_number("N") ||
				a->element == Atom::get_element_number("O") ||
				a->element == Atom::get_element_number("S"))
				a->properties.set_bit(Atom::res, 1);
		}
	}
}

// assign_atom_nneighbors:
void Molecule::assign_atom_nneighbors()
{
	// number of atoms bonded to each atom (neighbors)
	Atom* a;
	int n;
	for (vector<Atom*>::iterator
		p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a = *p;
		n = m_bonds_to_atom[a].size();
		if (n == 0) a->properties.set_bit(Atom::x0, 1);
		else if (n == 1) a->properties.set_bit(Atom::x1, 1);
		else if (n == 2) a->properties.set_bit(Atom::x2, 1);
		else if (n == 3) a->properties.set_bit(Atom::x3, 1);
		else a->properties.set_bit(Atom::x4, 1);
	}
}

// assign_atom_aromaticity:
// prop_rule-19: res~res~res~res~res~res~@1 ? ar6 ar6 ar6 ar6 ar6 ar6
// prop_rule-20: res~res~res~res~res~@1 ? ar5 ar5 ar5 ar5 ar5
// prop_rule-21: ar5,ar6 ? ar
void Molecule::assign_atom_aromaticity()
{
	// find atoms in five- and/or six-membered rings
	Atom* a0;
	Atom* a5;
	Atom* a6;
	vector<Atom*> a5s, a6s; // atoms bonded to 4-bond/5-bond reachable atom
	vector<Atom*>::iterator p, q, r, s;
	vector<vector<Atom*> > paths;
	for (p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (! a0->properties.get_bit(Atom::res)) continue;
		if (a0->properties.get_bit(Atom::ar5) ||
			a0->properties.get_bit(Atom::ar6)) continue; // already assigned
		get_reachable_atoms_from(a0, 4, paths);
			// paths from a to 4-bond reachable atoms

		// for each path from atom a
		for (vector<vector<Atom*> >::iterator
			pp = paths.begin(); pp != paths.end(); ++pp) {
			vector<Atom*>& path = *pp;
			for (q = path.begin(); q != path.end(); ++q)
				if (! (*q)->properties.get_bit(Atom::res)) break;
			if (q != path.end()) continue; // go to next path

			get_atoms_bonded_to_atom(path[4], a5s);
			// for each putative 5-bond reachable atom
			for (r = a5s.begin(); r != a5s.end(); ++r) {
				a5 = *r;
				if (! a5->properties.get_bit(Atom::res)) continue;
				if ((! a0->properties.get_bit(Atom::ar5)) && a5 == a0) {
					// atom a0 is 5-bond reachable from itself,
					// i.e., it is part of a 5-memberd ring
					for (q = path.begin(); q != path.end(); ++q) {
						(*q)->properties.set_bit(Atom::ar, 1);
						(*q)->properties.set_bit(Atom::ar5, 1);
					}
					break; // go to next path
				}
				if (a0->properties.get_bit(Atom::ar6)) continue;
					// already assigned to six-membered ring
				if (a5 == path[3]) continue; // a5 is actually 3-bond reachable
				get_atoms_bonded_to_atom(a5, a6s);
				// for each putative 6-bond reachable atom
				for (s = a6s.begin(); s != a6s.end(); ++s) {
					a6 = *s;
					if (a6 == a0) {
						// atom a0 is 6-bond reachable from itself,
						// i.e., it is part of a 6-membered ring
						for (q = path.begin(); q != path.end(); ++q) {
							(*q)->properties.set_bit(Atom::ar, 1);
							(*q)->properties.set_bit(Atom::ar6, 1);
						}
						a5->properties.set_bit(Atom::ar, 1);
						a5->properties.set_bit(Atom::ar6, 1);
						break; // go to next path
					}
				}
				if (s != a6s.end()) break; // go to next path
			} // to next 5-bond reachable atom
		} // to next path
	} // to next atom
}

// get_reachable_atoms_from:
void Molecule::get_reachable_atoms_from(
	Atom* a, int nsteps, vector<vector<Atom*> >& paths)
{
	if (m_bonds_to_atom.empty()) define_bonds_to_atom();
	paths.clear();
	vector<Atom*> path;
	visit_atom(a, nsteps, 0, paths, path);
}

// visit_atom:
void Molecule::visit_atom(
	Atom* a, int nsteps, int msteps,
	vector<vector<Atom*> >& paths, vector<Atom*>& path)
{
	// argument 'msteps' has the number of steps already tried
	vector<Atom*> my_path(path);
	my_path.push_back(a);
	if (nsteps == msteps) { // this is the end point
		paths.push_back(my_path);
		return;
	}

	vector<Bond*>& bs = m_bonds_to_atom[a];
	vector<Atom*>::iterator q;
	Atom* a2;
	for (vector<Bond*>::iterator p = bs.begin(); p != bs.end(); ++p) {
		if ((*p)->atom1 == a) a2 = (*p)->atom2;
		else a2 = (*p)->atom1;
		vector<Atom*>::const_iterator q;
		for (q = path.begin(); q != path.end(); ++q)
			if (a2 == *q) break; // a2 already visited
		if (q == path.end()) // a2 not visited
			visit_atom(a2, nsteps, msteps+1, paths, my_path);
	}

	return;
}

// assign_atom_electronegativity:
void Molecule::assign_atom_electronegativity()
{
	Atom* a;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a = *p;
		if (a->element == Atom::get_element_number("N") ||
			a->element == Atom::get_element_number("O"))
			a->properties.set_bit(Atom::neg, 1);
	}
}

// assign_atom_amide:
void Molecule::assign_atom_amide()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	vector<vector<Atom*> > paths;
	unsigned char e1, e2;
	unsigned int b1, b2;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (a0->element != Atom::get_element_number("N")) continue;
		get_reachable_atoms_from(a0, 2, paths);

		// for each path from a0 to 2-bond reachable atoms
		for (vector<vector<Atom*> >::const_iterator
			pp = paths.begin(); pp != paths.end(); ++pp) {
			const vector<Atom*>& path = *pp;
			a1 = path[1];
			a2 = path[2];
			char e1 = a1->element;
			char e2 = a2->element;

			if (! (e1 == Atom::get_element_number("C") ||
				e1 == Atom::get_element_number("P") ||
				e1 == Atom::get_element_number("S")))
				continue;
			if (! (e2 == Atom::get_element_number("O") ||
				e2 == Atom::get_element_number("S")))
				continue;

			b1 = bond_between(a0, a1);
			if (b1 != 1) continue;
			b2 = bond_between(a1, a2);
			// N-(C,P,S)=(O,S) or N-((P&x4),(S&x3))-((O,S)&x1)
			if ((b2 == 2) || (b2 == 1 && a2->properties.get_bit(Atom::x1) &&
				(( e1 == Atom::get_element_number("P") &&
					a1->properties.get_bit(Atom::x4)) ||
				(e1 == Atom::get_element_number("S") &&
					 a1->properties.get_bit(Atom::x3))))) {
				a0->properties.set_bit(Atom::namide, 1);
				break; // go to next atom
			}
		}
	}
}

// assign_atom_carboxylate:
void Molecule::assign_atom_carboxylate()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	vector<Atom*> as;
	int nas, i, j;
	unsigned int b1, b2;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (a0->element != Atom::get_element_number("C")) continue;
		get_atoms_bonded_to_atom(a0, as);
		nas = as.size();
		for (i = 0; i < (nas-1); ++i) {
			a1 = as[i];
			if (a1->element != Atom::get_element_number("O") &&
				a1->element != Atom::get_element_number("S"))
				continue;

			for (j = (i+1); j < nas; ++j) {
				a2 = as[j];
				if (a2->element != Atom::get_element_number("O") &&
					a2->element != Atom::get_element_number("S"))
					continue;

				b1 = bond_between(a0, a1);
				b2 = bond_between(a0, a2);
				// C(=O,S)(-(O,S)&x1)
				if ((b1 == 2 && b2 == 1 && a2->properties.get_bit(Atom::x1)) ||
					(b2 == 2 && b1 == 1 && a1->properties.get_bit(Atom::x1))) {
					a0->properties.set_bit(Atom::cx, 1);
					break;
				}
			}
			if (j < nas) break;
		}
	}
}

// assign_atom_pc_classes:
void Molecule::assign_atom_pc_classes()
{
	assign_atom_pc_class_default();
	assign_atom_special_nitrogen_types();
	assign_atom_quaternary_nitrogen();
	assign_atom_sp3_amines();
	assign_atom_conjugated_ns();
	assign_atom_basic_conjugated_sp3_ns();
	assign_atom_trisubstituted_sp2_amines();
	assign_atom_disubstituted_sp_amines();
	assign_atom_trisubstituted_os();
	assign_atom_miscellaneous_anions();
	assign_atom_diaminopyrimidine();
	assign_atom_ns_in_imidazole();
	assign_atom_sulfonamides_phosphonamides();
	assign_atom_carbonyl_hydroxide();
	assign_atom_adjacent_keto_enols();
	assign_atom_misc_buried();
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_pc_class_default:
// default classes
// rule-01: * > 7 (default)
// rule-02: O > 4 (set O's to acceptors)
// rule-03: N > 3 (set N's to donors)
// rule-04: C,Si > 6
// rule-05: S,Se > 6
// rule-06: P,As > 6
// rule-07: F,Cl,Br,I&x1 > 6 (halogens to hydrophobic)
// rule-08: O&x1-* > 5 * (hydroxide (this may be overwritten later))
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_pc_class_default()
{
	Atom* a;
	unsigned char pc_class, e;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a = *p;
		pc_class = Atom::NONE;
		e = a->element;
		if (e == Atom::get_element_number("O")) {
			pc_class = Atom::ACCEPTOR;
		} else if (e == Atom::get_element_number("N")) {
			pc_class = Atom::DONOR;
		} else if (
			e == Atom::get_element_number("C") ||
			e == Atom::get_element_number("Si") ||
			e == Atom::get_element_number("S") ||
			e == Atom::get_element_number("Se") ||
			e == Atom::get_element_number("P") ||
			e == Atom::get_element_number("As")
			) {
			pc_class = Atom::HYDROPHOBIC;
		}

		if (a->properties.get_bit(Atom::x1)) { // x1
			vector<Bond*>& bs = m_bonds_to_atom[a];
			if (
				e == Atom::get_element_number("F") ||
				e == Atom::get_element_number("Cl") ||
				e == Atom::get_element_number("Br") ||
				e == Atom::get_element_number("I")
				) { // (F,Cl,Br,I)&x1
				pc_class = Atom::HYDROPHOBIC;
			} else if (e == Atom::get_element_number("O") && bs[0]->v == 1) {
				// O&x1-* hydroxide
				pc_class = Atom::POLAR;
			}
		}

		a->pc_class = pc_class;
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_special_nitrogen_types:
// special nitrogen properties
// rule-09: N&x1&sp > 4 (cyano nitrogens are acceptors)
// rule-10: N&x2&sp2 > 4 (-N= are acceptors)
// rule-11: N&x1&sp2 > 5 (HN=X are polar (overwritten later by guanidium etc.)
// rule-12: *=N&x2=N&x1 > * 7 4 (azide)
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_special_nitrogen_types()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	vector<Atom*> as;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (a0->element != Atom::get_element_number("N")) continue;

		// N&x1&sp (cyano nitrogens are acceptors)
		if (a0->properties.get_bit(Atom::x1) && a0->properties.get_bit(Atom::sp))
			a0->pc_class = Atom::ACCEPTOR;

		// N&x2&sp2 (-N= are acceptors)
		if (a0->properties.get_bit(Atom::x2) && a0->properties.get_bit(Atom::sp2))
			a0->pc_class = Atom::ACCEPTOR;

		// N&x1&sp2 (HN=X are polar)
		if (a0->properties.get_bit(Atom::x1) && a0->properties.get_bit(Atom::sp2))
			a0->pc_class = Atom::POLAR;

		// *=N&x2=N&x1 -> * 7 4 (azide)
		if (a0->properties.get_bit(Atom::x2) &&
			m_bonds_to_atom[a0][0]->v == 2 && m_bonds_to_atom[a0][1]->v == 2) {
			get_atoms_bonded_to_atom(a0, as);
			if (as[0]->element == Atom::get_element_number("N") &&
				as[0]->properties.get_bit(Atom::x1)) {
				a0->pc_class = Atom::NONE;
				as[0]->pc_class = Atom::ACCEPTOR;
			}
			if (as[1]->element == Atom::get_element_number("N") &&
				as[1]->properties.get_bit(Atom::x1)) {
				a0->pc_class = Atom::NONE;
				as[1]->pc_class = Atom::ACCEPTOR;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_quaternary_nitrogen:
// rule-13: N,P,As&x4&sp3 > 1 (quaternary)
// rule-14: N,P,As&x4&sp3-O&x1 > 7 4 (turns off quaternary N-oxides)
// rule-15: N,P,As&x4&sp3-S&x1 > 7 6 (turns off quaternary N-sulfides)
// rule-16: N,P,As&x3,x4=O > 7 4 (dative P-oxide)
// rule-17: N,P,As&x3,x4=S > 7 6 (dative P-sulfides)
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_quaternary_nitrogen()
{
	Atom* a0;
	vector<Atom*> as;

	vector<Atom*> npas; // N, P, or As
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (a0->element == Atom::get_element_number("N") ||
			a0->element == Atom::get_element_number("P") ||
			a0->element == Atom::get_element_number("As"))
			npas.push_back(a0);
	}

	// N,P,As&x4&sp3 > 1 (quaternary)
	vector<Atom*> npas_x4_sp3;
	vector<Atom*> npas_x3x4;
	for (vector<Atom*>::iterator p = npas.begin(); p != npas.end(); ++p) {
		a0 = *p;
		if (a0->properties.get_bit(Atom::x3) || a0->properties.get_bit(Atom::x4)) {
			npas_x3x4.push_back(a0);
			if (a0->properties.get_bit(Atom::x4) && a0->properties.get_bit(Atom::sp3)) {
				a0->pc_class = Atom::CATION;
				npas_x4_sp3.push_back(a0);
			}
		}
	}

	// N,P,As&x4&sp3-O&x1 > 7 4 (turns off quaternary N-oxides)
	// N,P,As&x4&sp3-S&x1 > 7 6 (turns off quaternary N-sulfides)
	Atom* a1;
	unsigned int bt;
	unsigned char e1;
	for (vector<Atom*>::iterator
		p = npas_x4_sp3.begin(); p != npas_x4_sp3.end(); ++p) {
		a0 = *p;
		get_atoms_bonded_to_atom(a0, as);
		for (int i = 0; i < as.size(); ++i) {
			a1 = as[i];
			if (! a1->properties.get_bit(Atom::x1)) continue;
			bt = m_bonds_to_atom[a0][i]->v;
			if (bt != 1) continue;
			e1 = a1->element;
			if (e1 == Atom::get_element_number("O")) {
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::ACCEPTOR;
				break;
			} else if (e1 == Atom::get_element_number("S")) {
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::HYDROPHOBIC;
			}
		}
	}

	// N,P,As&(x3,x4)=O > 7 4 (dative P-oxide)
	// N,P,As&(x3,x4)=S > 7 6 (dative P-sulfides)
	for (vector<Atom*>::iterator
		p = npas_x3x4.begin(); p != npas_x3x4.end(); ++p) {
		a0 = *p;
		get_atoms_bonded_to_atom(a0, as);
		for (int i = 0; i < as.size(); ++i) {
			bt = m_bonds_to_atom[a0][i]->v;
			if (bt != 2) continue;
			a1 = as[i];
			e1 = a1->element;
			if (e1 == Atom::get_element_number("O")) {
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::ACCEPTOR;
			} else if (e1 == Atom::get_element_number("S")) {
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::HYDROPHOBIC;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_sp3_amines:
// sp3 amines and exceptions
// rule-18: N&x1,x2,x3&sp3!conj > 1 (any nonconjugated sp3 N is a basic amine)
// rule-19: N&x1,x2&sp3!conj-N&conj > 3 * (except next to a conjugated N)
// rule-20: N&x3&sp3!conj-N&conj > 7 * (except next to a conjugated N)
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_sp3_amines()
{
	Atom* a0;
	Atom* a1;
	vector<Atom*> as;

	// any nonconjugated sp3 N is a basic amine
	vector<Atom*> n_x1x2_sp3_nonconj;
	vector<Atom*> n_x3_sp3_nonconj;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (a0->element != Atom::get_element_number("N")) continue;
		if ((! a0->properties.get_bit(Atom::sp3)) ||
			a0->properties.get_bit(Atom::conj)) continue;
		if (a0->properties.get_bit(Atom::x1) ||
			a0->properties.get_bit(Atom::x2)) {
			a0->pc_class = Atom::CATION;
			n_x1x2_sp3_nonconj.push_back(a0);
		} else if (a0->properties.get_bit(Atom::x3)) {
			a0->pc_class = Atom::CATION;
			n_x3_sp3_nonconj.push_back(a0);
		}
	}

	// except next to a conjugated N
	for (vector<Atom*>::iterator
		p = n_x1x2_sp3_nonconj.begin(); p != n_x1x2_sp3_nonconj.end(); ++p) {
		a0 = *p;
		get_atoms_bonded_to_atom(a0, as);
		for (int i = 0; i < as.size(); ++i) {
			if (m_bonds_to_atom[a0][i]->v != 1) continue;
			a1 = as[i];
			if (a1->element == Atom::get_element_number("N") &&
				a1->properties.get_bit(Atom::conj)) {
				a0->pc_class = Atom::DONOR;
			}
		}
	}
	for (vector<Atom*>::iterator
		p = n_x3_sp3_nonconj.begin(); p != n_x3_sp3_nonconj.end(); ++p) {
		a0 = *p;
		get_atoms_bonded_to_atom(a0, as);
		for (int i = 0; i < as.size(); ++i) {
			if (m_bonds_to_atom[a0][i]->v != 1) continue;
			a1 = as[i];
			if (a1->element == Atom::get_element_number("N") &&
				a1->properties.get_bit(Atom::conj)) {
				a0->pc_class = Atom::NONE;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_conjugated_ns:
// rule-21: N&x1&sp3&conj > 3 (NH2-X= are donors)
// rule-22: N&x2&sp3&conj > 3 (-NH-= are donors)
// rule-23: N&x3&sp3&conj > 7 (-N(-)-X= are planar)
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_conjugated_ns()
{
	Atom* a0;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (a0->element != Atom::get_element_number("N")) continue;
		if (! (a0->properties.get_bit(Atom::sp3) &&
			a0->properties.get_bit(Atom::conj))) continue;
		if (a0->properties.get_bit(Atom::x1))
			a0->pc_class = Atom::DONOR; // x1; NH2-X= are donors
		if (a0->properties.get_bit(Atom::x2))
			a0->pc_class = Atom::DONOR; // x2; -NH-= are donors
		if (a0->properties.get_bit(Atom::x3))
			a0->pc_class = Atom::NONE; // x3; -N(-)-X= are planar
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_basic_conjugated_sp3_ns:
// there are special cases where conjugated sp3 N's are basic
// rule-24: C(=N&x1,x2&!ar,conj)-N&sp3&!ar,namide > 7 1 1 (amidine)
// rule-25: N&sp3&!ar,namide-C(=N&x1,x2&!ar,conj)-N&sp3&!ar,namide > 1 7 1 1 (guanidinium)
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_basic_conjugated_sp3_ns()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	vector<Atom*> as;
	vector<Atom*> n1s; // for =N&(x1,x2)&!(ar,conj)
	vector<Atom*> n2s; // for -N&sp3&!(ar,namide)
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (a0->element != Atom::get_element_number("C")) continue;
		//my $bonded = $bonds_to_atom->[$a0];
		get_atoms_bonded_to_atom(a0, as);
		n1s.clear();
		n2s.clear();
		for (int i = 0; i < as.size(); ++i) {
			a1 = as[i];
			if (a1->element != Atom::get_element_number("N")) continue;
			if (a1->properties.get_bit(Atom::ar)) continue;
			if (m_bonds_to_atom[a0][i]->v == 2 &&
				(! a1->properties.get_bit(Atom::conj)) &&
				(a1->properties.get_bit(Atom::x1) || a1->properties.get_bit(Atom::x2))) {
				// =N&(x1,x2)&!(ar,conj)
				n1s.push_back(a1);
			} else if (m_bonds_to_atom[a0][i]->v == 1 &&
				a1->properties.get_bit(Atom::sp3) &&
				(! a1->properties.get_bit(Atom::namide))) {
				// -N&sp3&!(ar,namide)
				n2s.push_back(a1);
			}
		}

		// C(=N&(x1,x2)&!(ar,conj))-N&sp3&!(ar,namide) > 7 1 1 (amidine)
		// N&sp3&!(ar,namide)-C(=N&(x1,x2)&!(ar,conj))-N&sp3&!(ar,namide) > 1 7 1 1 (guanidinium)
		for (vector<Atom*>::iterator p = n1s.begin(); p != n1s.end(); ++p) {
			a1 = *p;
			for (vector<Atom*>::iterator q = n2s.begin(); q != n2s.end(); ++q) {
				a2 = *q;
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::CATION;
				a2->pc_class = Atom::CATION;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_trisubstituted_sp2_amines:
// trisubstituted sp2 amines are cations with exceptions
// rule-26: N,P,As&x3&sp2 > 1 (trisubstituted aromatic amine)
// rule-27: N,P,As&x3&sp2-O&x1 > 7 4 (turns off tertiary in nitro and N-oxides)
// rule-28: N,P,As&x3&sp2-S&x1 > 7 6 (turns off tertiary in N-sulfides)
// rule-29: N&x3&sp(=O)=O > 7 4 4 (turns off tertiary N in alternative representation in dative nitro)
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_trisubstituted_sp2_amines()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	vector<Atom*> as;
	unsigned char e0, e1;

	// N,P,As&x3&sp2 > 1 (trisubstituted aromatic amine)
	vector<Atom*> npas_x3_sp2;
	vector<Atom*> n_x3_sp;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (! a0->properties.get_bit(Atom::x3)) continue;
		e0 = a0->element;
		if ((e0 == Atom::get_element_number("N") ||
			e0 == Atom::get_element_number("P") ||
			e0 == Atom::get_element_number("As")) &&
			a0->properties.get_bit(Atom::sp2)) {
			a0->pc_class = Atom::CATION;
			npas_x3_sp2.push_back(a0);
		} else if (e0 == Atom::get_element_number("N") &&
					a0->properties.get_bit(Atom::sp)) {
			n_x3_sp.push_back(a0);
		}
	}

	// N,P,As&x3&sp2-O&x1 > 7 4 (turns off tertiary in nitro and N-oxides)
	// N,P,As&x3&sp2-S&x1 > 7 6 (turns off tertiary in N-sulfides)
	for (vector<Atom*>::iterator
		p = npas_x3_sp2.begin(); p != npas_x3_sp2.end(); ++p) {
		a0 = *p;
		get_atoms_bonded_to_atom(a0, as);
		for (int i = 0; i < as.size(); ++i) {
			if (m_bonds_to_atom[a0][i]->v != 1) continue; // single-bonded
			a1 = as[i];
			if (! a1->properties.get_bit(Atom::x1)) continue;
			e1 = a1->element;
			if (e1 == Atom::get_element_number("O")) {
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::ACCEPTOR;
			} else if (e1 == Atom::get_element_number("S")) {
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::HYDROPHOBIC;
				break; // -S&x1 has higher priority
			}
		}
	}

	// N&x3&sp(=O)=O > 7 4 4 (turns off tertiary N in alternative representation in dative nitro)
	vector<Atom*> db_os; // double-bonded O's
	for (vector<Atom*>::iterator
		p = n_x3_sp.begin(); p != n_x3_sp.end(); ++p) {
		a0 = *p;
		get_atoms_bonded_to_atom(a0, as);
		db_os.clear();
		for (int i = 0; i < as.size(); ++i) {
			a1 = as[i];
			if (m_bonds_to_atom[a0][i]->v == 2 &&
				a1->element == Atom::get_element_number("O"))
				db_os.push_back(a1);
		}
		for (int i = 0; i < (db_os.size()-1); ++i) {
			a1 = db_os[i];
			for (int j = (i+1); j < db_os.size(); ++j) {
				a2 = db_os[j];
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::ACCEPTOR;
				a2->pc_class = Atom::ACCEPTOR;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_disubstituted_sp_amines:
// disubstituted sp amines and exceptions
// rule-30: N&x2&sp > 1 (diazonium)
// rule-31: N&x2&sp-O&x1 > 7 4 (nitrile oxide)
// rule-32: N&x2&sp#C&x1 > 7 4 (isonitrile)
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_disubstituted_sp_amines()
{
	Atom* a0;
	Atom* a1;
	vector<Atom*> as;

	// N&x2&sp > 1 (diazonium)
	vector<Atom*> n_x2_sp;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (! (a0->properties.get_bit(Atom::x2) &&
				a0->properties.get_bit(Atom::sp))) continue;
		if (a0->element == Atom::get_element_number("N")) {
			a0->pc_class = Atom::CATION;
			n_x2_sp.push_back(a0);
		}
	}

	for (vector<Atom*>::iterator p = n_x2_sp.begin(); p != n_x2_sp.end(); ++p) {
		a0 = *p;
		get_atoms_bonded_to_atom(a0, as);

		// N&x2&sp-O&x1 > 7 4 (nitrile oxide)
		for (int i = 0; i < as.size(); ++i) {
			if (m_bonds_to_atom[a0][i]->v != 1) continue; // single-bonded
			a1 = as[i];
			if (! a1->properties.get_bit(Atom::x1)) continue;
			if (a1->element == Atom::get_element_number("O")) {
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::ACCEPTOR;
			}
		}

		// N&x2&sp#C&x1 > 7 4 (isonitrile) (higher priority)
		for (int i = 0; i < as.size(); ++i) {
			if (m_bonds_to_atom[a0][i]->v != 3) continue; // triple-bonded
			a1 = as[i];
			if (! a1->properties.get_bit(Atom::x1)) continue;
			if (a1->element == Atom::get_element_number("C")) {
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::ACCEPTOR;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_trisubstituted_os:
// trisubstituted oxygens, sulfurs and exceptions
// rule-33: O,S&x3&sp3 > 1 (trisubstituted)
// rule-34: O,S&x2&sp2 > 1 (trivalent)
// rule-35: S&x3,x4&sp3-O&x1 > 7 4 (turns off sulfoxides and sulfones)
// rule-36: S&x3,x4&sp3-S&x1 > 7 6
// rule-37: S&x3,x4&sp2sp=O > 7 4 (turns off dative sulfoxides and sulfones)
// rule-38: S&x3,x4&sp2sp=S > 7 6
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_trisubstituted_os()
{
	Atom* a0;
	Atom* a1;
	unsigned char e0, e1;
	vector<Atom*> as;

	// O,S&x3&sp3 > 1 (trisubstituted)
	// O,S&x2&sp2 > 1 (trivalent)
	vector<Atom*> s_x3x4;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		e0 = a0->element;
		if (! (e0 == Atom::get_element_number("O") ||
				e0 == Atom::get_element_number("S"))) continue;
		if (a0->properties.get_bit(Atom::x3) && a0->properties.get_bit(Atom::sp3))
			a0->pc_class = Atom::CATION;
		else if (a0->properties.get_bit(Atom::x2) && a0->properties.get_bit(Atom::sp2))
			a0->pc_class = Atom::CATION;

		if (e0 == Atom::get_element_number("S") &&
			(a0->properties.get_bit(Atom::x3) || a0->properties.get_bit(Atom::x4)))
			s_x3x4.push_back(a0); // S&x3,x4
	}

	// S&x3,x4&sp3-O&x1 > 7 4 (turns off sulfoxides and sulfones)
	// S&x3,x4&sp3-S&x1 > 7 6
	// S&x3,x4&sp2sp=O > 7 4 (turns off dative sulfoxides and sulfones)
	// S&x3,x4&sp2sp=S > 7 6
	unsigned int bond_type;
	for (vector<Atom*>::iterator p = s_x3x4.begin(); p != s_x3x4.end(); ++p) {
		a0 = *p;
		get_atoms_bonded_to_atom(a0, as);
		bool sp3 = a0->properties.get_bit(Atom::sp3);
		bool sp2sp = (a0->properties.get_bit(Atom::sp) ||
						a0->properties.get_bit(Atom::sp2));
		for (int i = 0; i < as.size(); ++i) {
			a1 = as[i];
			e1 = a1->element;
			bond_type = m_bonds_to_atom[a0][i]->v;
			if (sp3 && bond_type == 1 &&
				e1 == Atom::get_element_number("O") &&
				a1->properties.get_bit(Atom::x1)) { // -O&x1
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::ACCEPTOR;
			} else if (sp3 && bond_type == 1 &&
				e1 == Atom::get_element_number("S") &&
				a1->properties.get_bit(Atom::x1)) { // -S&x1
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::HYDROPHOBIC;
			} else if (sp2sp && bond_type == 2 &&
				e1 == Atom::get_element_number("O")) { // =O
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::ACCEPTOR;
			} else if (sp2sp && bond_type == 2 &&
				e1 == Atom::get_element_number("S")) { // =S
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::HYDROPHOBIC;
			}
		}
	}
}

// assign_atom_miscellaneous_anions:
void Molecule::assign_atom_miscellaneous_anions()
{
	assign_atom_ma_conjugated_sulhydryls();
	assign_atom_ma_carboxylates_sequiv();
	assign_atom_ma_phosphate_arsenate();
	assign_atom_ma_muscimol();
	assign_atom_ma_tetrazole();
}

// rule-39: S&x1-ar6 > 2 * (conjugated sulhydryls)
void Molecule::assign_atom_ma_conjugated_sulhydryls()
{
	// S&x1-ar6 > 2 * (conjugated sulhydryls)
	Bond* b;
	Atom* a0;
	Atom* a1;
	for (vector<Bond*>::iterator p = m_bonds.begin(); p != m_bonds.end(); ++p) {
		b = *p;
		a0 = b->atom1;
		a1 = b->atom2;
		if (b->v != 1) continue; // single bond
		if (a0->element == Atom::get_element_number("S") &&
			a0->properties.get_bit(Atom::x1) && a1->properties.get_bit(Atom::ar6)) {
			a0->pc_class = Atom::ANION;
		}
		if (a1->element == Atom::get_element_number("S") &&
			a1->properties.get_bit(Atom::x1) && a0->properties.get_bit(Atom::ar6)) {
			a1->pc_class = Atom::ANION;
		}
	}
}

// assign_atom_ma_carboxylates_sequiv:
// rule-40: C&x3(=O,S)-O,S&x1 > 7 2 2 (carboxylates and S equivalent)
void Molecule::assign_atom_ma_carboxylates_sequiv()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	vector<Atom*> as;
	vector<Atom*> a1s; // =O,S
	vector<Atom*> a2s; // -(O,S)&x1
	unsigned char e1;

	// C&x3(=O,S)-(O,S)&x1 > 7 2 2 (carboxylates and S equivalent)
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (! a0->properties.get_bit(Atom::x3)) continue;
		if (a0->element != Atom::get_element_number("C")) continue;
		get_atoms_bonded_to_atom(a0, as);
		a1s.clear();
		a2s.clear();
		for (int i = 0; i < as.size(); ++i) {
			a1 = as[i];
			e1 = a1->element;
			if (! (e1 == Atom::get_element_number("O") ||
					e1 == Atom::get_element_number("S"))) continue;
			if (m_bonds_to_atom[a0][i]->v == 2) {
				a1s.push_back(a1);
			} else if (m_bonds_to_atom[a0][i]->v == 1 &&
				a1->properties.get_bit(Atom::x1)) {
				a2s.push_back(a1);
			}
		}

		for (vector<Atom*>::iterator p = a1s.begin(); p != a1s.end(); ++p) {
			a1 = *p;
			for (vector<Atom*>::iterator q = a2s.begin(); q != a2s.end(); ++q) {
				a2 = *q;
				a0->pc_class = Atom::NONE;
				a1->pc_class = Atom::ANION;
				a2->pc_class = Atom::ANION;
			}
		}
	}
}

// assign_atom_ma_phosphate_arsenate:
// rule-41: P,As&x3,x4(~O,S&x1)~O,S&x1 > 7 2 2 (phosphate/arsenate monoanion and S equiv)
// rule-42: P,As&x4(~O,S&x1)(~O,S&x1)~O,S&x1 > 7 2 2 2 (phosphate/arsenate dianion and S equiv)
// rule-43: S,Se&x3(~O,S&x1)~O,S&x1 > 7 2 2 (sulfite/selenite anion and S equivalent)
// rule-44: S,Se&x4(~O,S&x1)(~O,S&x1)~O,S&x1 > 7 2 2 2
void Molecule::assign_atom_ma_phosphate_arsenate()
{
	Atom* a0;
	Atom* a1;
	vector<Atom*> as;
	vector<Atom*> osx1; // ~(O,S)&x1
	unsigned char e0, e1;

	// P,As&x3,x4(~O,S&x1)~O,S&x1 > 7 2 2 (phosphate/arsenate monoanion and S equiv)
	// P,As&x4(~O,S&x1)(~O,S&x1)~O,S&x1 > 7 2 2 2 (phosphate/arsenate dianion and S equiv)
	// S,Se&x3(~O,S&x1)~O,S&x1 > 7 2 2 (sulfite/selenite anion and S equivalent)
	// S,Se&x4(~O,S&x1)(~O,S&x1)~O,S&x1 > 7 2 2 2
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (! (a0->properties.get_bit(Atom::x3) ||
				a0->properties.get_bit(Atom::x4))) continue;
		e0 = a0->element;
		if (! (e0 == Atom::get_element_number("P") ||
				e0 == Atom::get_element_number("As") ||
				e0 == Atom::get_element_number("S") ||
				e0 == Atom::get_element_number("Se"))) continue;
		get_atoms_bonded_to_atom(a0, as);
		osx1.clear();
		for (int i = 0; i < as.size(); ++i) {
			a1 = as[i];
			e1 = a1->element;
			if ((e1 == Atom::get_element_number("O") ||
					e1 == Atom::get_element_number("S")) &&
					a1->properties.get_bit(Atom::x1)) {
				osx1.push_back(a1);
			}
		}
		if (osx1.size() < 2) continue;
		if ((e0 == Atom::get_element_number("P") ||
				e0 == Atom::get_element_number("As")) && osx1.size() == 2) {
			a0->pc_class = Atom::NONE;
			osx1[0]->pc_class = Atom::ANION;
			osx1[1]->pc_class = Atom::ANION;
		} else if ((e0 == Atom::get_element_number("S") ||
				e0 == Atom::get_element_number("Se")) &&
				a0->properties.get_bit(Atom::x3) && osx1.size() == 2) {
			a0->pc_class = Atom::NONE;
			osx1[0]->pc_class = Atom::ANION;
			osx1[1]->pc_class = Atom::ANION;
		} else if (osx1.size() == 3 && a0->properties.get_bit(Atom::x4)) {
			a0->pc_class = Atom::NONE;
			osx1[0]->pc_class = Atom::ANION;
			osx1[1]->pc_class = Atom::ANION;
			osx1[2]->pc_class = Atom::ANION;
		}
	}
}

// assign_atom_ma_muscimol:
// rule-45: C&x3&sp2(~O&x1)~N&x2~O~*~*~@1 > 6 2 2 4 * * (anion as in muscimol)
void Molecule::assign_atom_ma_muscimol()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	Atom* a3;
	Atom* a4;
	Atom* a5;
	Atom* a6;
	vector<Atom*> as;
	vector<Atom*> a5s;
	vector<vector<Atom*> > paths;
	bool is_ring;

	// C&x3&sp2(~O&x1)~N&x2~O~*~*~@1 > 6 2 2 4 * * (anion as in muscimol)
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (! (a0->properties.get_bit(Atom::x3) &&
			a0->properties.get_bit(Atom::sp2) &&
			a0->element == Atom::get_element_number("C"))) continue;
		get_reachable_atoms_from(a0, 4, paths);
			// paths to 4-bond reachable atoms
		for (vector<vector<Atom*> >::iterator
			p = paths.begin(); p != paths.end(); ++p) {
			vector<Atom*>& path = *p;
			a1 = path[1]; // 1-bond reachable atom in this path
			a2 = path[2]; // 2-bond reachable atom in this path
			a3 = path[3]; // 3-bond reachable atom in this path
			a4 = path[4]; // 4-bond reachable atom in this path
			get_atoms_bonded_to_atom(a4, a5s);

			is_ring = false;
			for (vector<Atom*>::iterator
				q = a5s.begin(); q != a5s.end(); ++q) {
				// for each putative 5-bond reachable atom
				if (*q == a0) {
					is_ring = true;
					break;
				}
			}
			if (! is_ring) continue;

			a6 = 0;
			get_atoms_bonded_to_atom(a0, as);
			for (vector<Atom*>::iterator q = as.begin(); q != as.end(); ++q) {
				if ((*q)->element == Atom::get_element_number("O") &&
					(*q)->properties.get_bit(Atom::x1)) {
					a6 = *q;
					break;
				}
			}
			if (!a6) continue;

			if (a1->element == Atom::get_element_number("N") &&
				a1->properties.get_bit(Atom::x2) &&
				a2->element == Atom::get_element_number("O")) {
				a0->pc_class = Atom::HYDROPHOBIC;
				a6->pc_class = Atom::ANION;
				a1->pc_class = Atom::ANION;
				a2->pc_class = Atom::ACCEPTOR;
			} else if (a4->element == Atom::get_element_number("N") &&
				a4->properties.get_bit(Atom::x2) &&
				a3->element == Atom::get_element_number("O")) {
				a0->pc_class = Atom::HYDROPHOBIC;
				a6->pc_class = Atom::ANION;
				a4->pc_class = Atom::ANION;
				a3->pc_class = Atom::ACCEPTOR;
			}
		}
	}
}

// assign_atom_ma_tetrazole:
// rule-46: C=N&x2-N&x2=N&x2-N&x2-@1 > 6 2 2 2 2 (tetrazole)
// rule-47: C=N&x2-N&x2-N&x2=N&x2-@1 > 6 2 2 2 2 (tetrazole tautomer)
void Molecule::assign_atom_ma_tetrazole()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	Atom* a3;
	Atom* a4;
	vector<Atom*> a5s;
	vector<vector<Atom*> > paths;
	bool is_ring;
	unsigned int bt1, bt2, bt3, bt4, bt5;

	// C=N&x2-N&x2=N&x2-N&x2-@1 > 6 2 2 2 2 (tetrazole)
	// C=N&x2-N&x2-N&x2=N&x2-@1 > 6 2 2 2 2 (tetrazole tautomer)
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (a0->element != Atom::get_element_number("C")) continue;
		get_reachable_atoms_from(a0, 4, paths);
			// paths to 4-bond reachable atoms
		for (vector<vector<Atom*> >::iterator
			p = paths.begin(); p != paths.end(); ++p) {
			vector<Atom*>& path = *p;

			a1 = path[1]; // 1-bond reachable atom in this path
			if (! (a1->element == Atom::get_element_number("N") &&
					a1->properties.get_bit(Atom::x2))) continue;
			a2 = path[2]; // 2-bond reachable atom in this path
			if (! (a2->element == Atom::get_element_number("N") &&
					a2->properties.get_bit(Atom::x2))) continue;
			a3 = path[3]; // 3-bond reachable atom in this path
			if (! (a3->element == Atom::get_element_number("N") &&
					a3->properties.get_bit(Atom::x2))) continue;
			a4 = path[4]; // 4-bond reachable atom in this path
			if (! (a4->element == Atom::get_element_number("N") &&
					a4->properties.get_bit(Atom::x2))) continue;

			get_atoms_bonded_to_atom(a4, a5s);
				// atoms bonded to the 4-bond reachable atom
			is_ring = false;
			for (vector<Atom*>::iterator
				q = a5s.begin(); q != a5s.end(); ++q) {
				// for each putative 5-bond reachable atom
				if (*q == a0) {
					is_ring = true;
					break;
				}
			}
			if (! is_ring) continue;

			bt1 = bond_between(a0, a1);
			bt2 = bond_between(a1, a2);
			bt3 = bond_between(a2, a3);
			bt4 = bond_between(a3, a4);
			bt5 = bond_between(a4, a0);
			if ((bt1 == 2 && bt2 == 1 && bt3 == 2 && bt4 == 1 && bt5 == 1) ||
				(bt1 == 1 && bt2 == 1 && bt3 == 2 && bt4 == 1 && bt5 == 2) ||
				(bt1 == 2 && bt2 == 1 && bt3 == 1 && bt4 == 2 && bt5 == 1) ||
				(bt1 == 1 && bt2 == 2 && bt3 == 1 && bt4 == 1 && bt5 == 2)) {
				a0->pc_class = Atom::HYDROPHOBIC;
				a1->pc_class = Atom::ANION;
				a2->pc_class = Atom::ANION;
				a3->pc_class = Atom::ANION;
				a4->pc_class = Atom::ANION;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// diaminopyrimidine
// rule-48: N&x2&sp2~C&x3&sp2(-N&x1)~N&x2&sp2~C&x3&sp2(-N&x1)~*~*~@1 > 5 * 3 5 * 3 * *
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_diaminopyrimidine()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	Atom* a3;
	Atom* a4;
	Atom* a5;
	vector<Atom*> a6s;
	vector<Atom*> a3nx1s;
	vector<Atom*> a1nx1s;
	vector<Atom*> a5nx1s;
	vector<vector<Atom*> > paths;
	bool is_ring;

	// N&x2&sp2~C&x3&sp2(-N&x1)~N&x2&sp2~C&x3&sp2(-N&x1)~*~*~@1 > 5 * 3 5 * 3 * *
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (! (a0->element == Atom::get_element_number("N") &&
				a0->properties.get_bit(Atom::x2) &&
				a0->properties.get_bit(Atom::sp2))) continue;

		get_reachable_atoms_from(a0, 5, paths);
			// paths to 5-bond reachable atoms
		for (vector<vector<Atom*> >::iterator
			p = paths.begin(); p != paths.end(); ++p) {
			vector<Atom*>& path = *p;
			a5 = path[5];
			get_atoms_bonded_to_atom(a5, a6s);
				// atoms bonded to the 5-bond reachable atom
			is_ring = false;
			for (vector<Atom*>::iterator
				q = a6s.begin(); q != a6s.end(); ++q) {
				// for each putative 6-bond reachable atom
				if (*q == a0) {
					is_ring = true;
					break;
				}
			}
			if (! is_ring) continue;
			a1 = path[1]; // 1-bond reachable atom in this path
			a2 = path[2]; // 2-bond reachable atom in this path
			a3 = path[3]; // 3-bond reachable atom in this path
			a4 = path[4]; // 4-bond reachable atom in this path
			// needs to check both 0->1->2->3->4->5 and 0->5->4->3->2->1
			// position 3 is C&x3&sp2 with (-N&x1) in any case
			if (! (a3->element == Atom::get_element_number("C") &&
					a3->properties.get_bit(Atom::x3) &&
					a3->properties.get_bit(Atom::sp2)))
					continue;

			get_sb_nx1s(a3, a3nx1s);
			if (a3nx1s.empty()) continue;
			Atom* a3nx1 = a3nx1s[0];
			// for 0->1->2->3->4->5
			if (a1->element == Atom::get_element_number("C") &&
				a1->properties.get_bit(Atom::x3) &&
				a1->properties.get_bit(Atom::sp2) &&
				a2->element == Atom::get_element_number("N") &&
				a2->properties.get_bit(Atom::x2) &&
				a2->properties.get_bit(Atom::sp2)) {
				get_sb_nx1s(a1, a1nx1s);
				if (! a1nx1s.empty()) {
					a0->pc_class = Atom::POLAR;
					a1nx1s[0]->pc_class = Atom::DONOR;
					a2->pc_class = Atom::POLAR;
					a3nx1->pc_class = Atom::DONOR;
				}
			}
			// for 0->5->4->3->2->1
			if (a5->element == Atom::get_element_number("C") &&
				a5->properties.get_bit(Atom::x3) &&
				a5->properties.get_bit(Atom::sp2) &&
				a4->element == Atom::get_element_number("N") &&
				a4->properties.get_bit(Atom::x2) &&
				a4->properties.get_bit(Atom::sp2)) {
				get_sb_nx1s(a5, a5nx1s);
				if (! a5nx1s.empty()) {
					a0->pc_class = Atom::POLAR;
					a5nx1s[0]->pc_class = Atom::DONOR;
					a4->pc_class = Atom::POLAR;
					a3nx1->pc_class = Atom::DONOR;
				}
			}
		}
	}
}

// get_sb_nx1s:
// returns an array of serial numbers of N&x1's single-bonded to atom $a0
void Molecule::get_sb_nx1s(Atom* a0, vector<Atom*>& nx1s)
{
	nx1s.clear();
	Atom* a1;
	if (a0) return;
	vector<Atom*> as;
	get_atoms_bonded_to_atom(a0, as);
	for (int i = 0; i < as.size(); ++i) {
		if (m_bonds_to_atom[a0][i]->v != 1) continue; // must be single-bonded
		a1 = as[i];
		if (a1->element == Atom::get_element_number("N") && 
			a1->properties.get_bit(Atom::x1)) {
			nx1s.push_back(a1);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_ns_in_imidazole:
// make N's in imidazole equivalent
// rule-49: N&x2&sp3-C=N&x2-C=C-@1 > 5 * 5 * *
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_ns_in_imidazole()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	Atom* a3;
	Atom* a4;
	vector<Atom*> a5s;
	vector<vector<Atom*> > paths;
	bool is_ring;
	unsigned int bt1, bt2, bt3, bt4, bt5;

	// N&x2&sp3-C=N&x2-C=C-@1 > 5 * 5 * *
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (! (a0->element == Atom::get_element_number("N") &&
				a0->properties.get_bit(Atom::x2) &&
				a0->properties.get_bit(Atom::sp3))) continue;

		get_reachable_atoms_from(a0, 4, paths);
			// paths to 4-bond reachable atoms
		for (vector<vector<Atom*> >::iterator
			p = paths.begin(); p != paths.end(); ++p) {
			vector<Atom*>& path = *p;
			a4 = path[4]; // 4-bond reachable atom
			get_atoms_bonded_to_atom(a4, a5s);
				// atoms bonded to the 4-bond reachable atom
			is_ring = false;
			for (vector<Atom*>::iterator
				q = a5s.begin(); q != a5s.end(); ++q) {
				// for each putative 5-bond reachable atom
				if (*q == a0) {
					is_ring = true;
					break;
				}
			}
			if (! is_ring) continue; // try next path

			// needs to check both 0->1->2->3->4 and 0->4->3->2->1
			a1 = path[1]; // 1-bond reachable atom
			a2 = path[2]; // 2-bond reachable atom
			a3 = path[3]; // 3-bond reachable atom
			bt1 = bond_between(a0, a1);
			bt5 = bond_between(a0, a4);
			if (! (bt1 == 1 && bt5 == 1)) continue;
			bt2 = bond_between(a1, a2);
			bt4 = bond_between(a4, a3);
			if (! (bt2 == 2 && bt4 == 2)) continue;
			bt3 = bond_between(a2, a3);
			if (bt3 != 1) continue;
			// for 0->1->2->3->4
			if (a1->element == Atom::get_element_number("C") &&
				a2->element == Atom::get_element_number("N") &&
				a2->properties.get_bit(Atom::x2) &&
				a3->element == Atom::get_element_number("C") &&
				a4->element == Atom::get_element_number("C")) {
				a0->pc_class = Atom::POLAR;
				a2->pc_class = Atom::POLAR;
			}
			// for 0->4->3->2->1
			if (a4->element == Atom::get_element_number("C") &&
				a3->element == Atom::get_element_number("N") &&
				a3->properties.get_bit(Atom::x2) &&
				a2->element == Atom::get_element_number("C") &&
				a1->element == Atom::get_element_number("C")) {
				a0->pc_class = Atom::POLAR;
				a3->pc_class = Atom::POLAR;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_sulfonamides_phosphonamides:
// stabilized sulfonamides and phosphonamides
// rule-50: ar-S,P&x4(~O,S&x1)(~O,S&x1)~N&x1 > * 7 4 4 5
// rule-51: *-S,P&x4(~O,S&x1)(~O,S&x1)-N&x2-sp2,sp > * 7 4 4 5 *
// rule-52: ar-S,P&x4(~O,S&x1)(~O,S&x1)-N&x2-sp2,sp > * 7 2 2 2 *
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_sulfonamides_phosphonamides()
{
	Atom* a0;
	Atom* a1;
	Atom* a4;
	Atom* atmp;
	vector<Atom*> as;
	vector<Atom*> osx1s;
	vector<Atom*> ars_sb;
	vector<Atom*> sbs;
	vector<Atom*> nx1s;
	vector<Atom*> nx2_sp2sps;
	vector<Atom*> bonds_nx2;

	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p; // S,P&x4 = a0
		if (! ((a0->element == Atom::get_element_number("S") ||
				a0->element == Atom::get_element_number("P")) &&
				a0->properties.get_bit(Atom::x4))) continue;
		get_atoms_bonded_to_atom(a0, as);
		get_bonded_osx1s(a0, osx1s);
		if (osx1s.size() < 2) continue;

		ars_sb.clear();
		sbs.clear();
		nx1s.clear();
		nx2_sp2sps.clear();
		for (int i = 0; i < as.size(); ++i) {
			atmp = as[i];
			if (atmp->element == Atom::get_element_number("N") &&
				atmp->properties.get_bit(Atom::x1))
				nx1s.push_back(atmp);
			if (m_bonds_to_atom[a0][i]->v == 1) {
				sbs.push_back(atmp); // *-
				if (atmp->properties.get_bit(Atom::ar))
					ars_sb.push_back(atmp); // ar-
				if (atmp->element == Atom::get_element_number("N") &&
					atmp->properties.get_bit(Atom::x2)) { // -N&x2
					get_atoms_bonded_to_atom(atmp, bonds_nx2);
					if (bonds_nx2[0] == a0 &&
						(bonds_nx2[1]->properties.get_bit(Atom::sp) ||
							bonds_nx2[1]->properties.get_bit(Atom::sp2))) {
						nx2_sp2sps.push_back(atmp);
					} else if (bonds_nx2[1] == a0 &&
						(bonds_nx2[0]->properties.get_bit(Atom::sp) ||
							bonds_nx2[0]->properties.get_bit(Atom::sp2))) {
						nx2_sp2sps.push_back(atmp);
					}
				}
			}
		}

		// ar-S,P&x4(~O,S&x1)(~O,S&x1)~N&x1 > * 7 4 4 5
		for (vector<Atom*>::iterator p = ars_sb.begin(); p != ars_sb.end(); ++p) {
			a1 = *p;
			for (vector<Atom*>::iterator q = nx1s.begin(); q != nx1s.end(); ++q) {
				a4 = *q;
				if (a4 == a1) continue;
				for (vector<Atom*>::iterator
					r = osx1s.begin(); r != osx1s.end(); ++r) {
					(*r)->pc_class = Atom::ACCEPTOR;
				}
				a0->pc_class = Atom::NONE;
				a4->pc_class = Atom::POLAR;
			}
		}

		// *-S,P&x4(~O,S&x1)(~O,S&x1)-N&x2-sp2,sp > * 7 4 4 5 *
		for (vector<Atom*>::iterator p = sbs.begin(); p != sbs.end(); ++p) {
			a1 = *p;
			for (vector<Atom*>::iterator
				q = nx2_sp2sps.begin(); q != nx2_sp2sps.end(); ++q) {
				a4 = *q;
				if (a4 == a1) continue;
				for (vector<Atom*>::iterator
					r = osx1s.begin(); r != osx1s.end(); ++r) {
					(*r)->pc_class = Atom::ACCEPTOR;
				}
				a0->pc_class = Atom::NONE;
				a4->pc_class = Atom::POLAR;
			}
		}

		// ar-S,P&x4(~O,S&x1)(~O,S&x1)-N&x2-sp2,sp > * 7 2 2 2 *
		for (vector<Atom*>::iterator p = ars_sb.begin(); p != ars_sb.end(); ++p) {
			a1 = *p;
			for (vector<Atom*>::iterator
				q = nx2_sp2sps.begin(); q != nx2_sp2sps.end(); ++q) {
				a4 = *q;
				if (a4 == a1) continue;
				for (vector<Atom*>::iterator
					r = osx1s.begin(); r != osx1s.end(); ++r) {
					(*r)->pc_class = Atom::ANION;
				}
				a0->pc_class = Atom::NONE;
				a4->pc_class = Atom::ANION;
			}
		}
	}
}

// get_bonded_osx1s:
// returns an array of serial numbers of O,S&x1's bonded to a0
void Molecule::get_bonded_osx1s(Atom* a0, vector<Atom*>& osx1s)
{
	Atom* a1;
	vector<Atom*> as;
	get_atoms_bonded_to_atom(a0, as);
	osx1s.clear();
	for (int i = 0; i < as.size(); ++i) {
		a1 = as[1];
		if ((a1->element == Atom::get_element_number("O") ||
				a1->element == Atom::get_element_number("S")) &&
				a1->properties.get_bit(Atom::x1)) {
			osx1s.push_back(a1);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_carbonyl_hydroxide:
// carbonyl-hydroxide tautomers
// rule-53: O&x1-C!cx=sp2-C!cx=O > 5 * * * 5
// rule-54: O&x1-C&ar6~ar6~C&ar6=O&x1 > 5 * * * 5
// rule-55: O&x1-C&ar6~ar6~ar6~ar6~C&ar6=O&x1 > 5 * * * * * 5
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_carbonyl_hydroxide()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	Atom* a3;
	Atom* a4;
	Atom* a5;
	Atom* a6;
	vector<vector<Atom*> > paths;
	unsigned int bt1, bt2, bt3, bt4, bt5, bt6;

	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (! (a0->element == Atom::get_element_number("O") &&
				a0->properties.get_bit(Atom::x1))) continue;
		get_reachable_atoms_from(a0, 4, paths);
			// paths to 4-bond reachable atoms
		for (vector<vector<Atom*> >::iterator
			p = paths.begin(); p != paths.end(); ++p) {
			vector<Atom*>& path = *p;

			a4 = path[4]; // 4-bond reachable atom
			if (a4->element != Atom::get_element_number("O")) continue;

			a3 = path[3]; // 3-bond reachable atom
			bt4 = bond_between(a3, a4);
			if (bt4 != 2) continue; // must be =O
			if (a3->element != Atom::get_element_number("C")) continue;

			a1 = path[1]; // 1-bond reachable atom
			if (a1->element != Atom::get_element_number("C")) continue;
			bt1 = bond_between(a0, a1);
			if (bt1 != 1) continue; // must be -C

			a2 = path[2]; // 2-bond reachable atom
			bt2 = bond_between(a1, a2);
			bt3 = bond_between(a2, a3);

			// now, O&x1-C~*~C=O satisfied
			// O&x1-C!cx=sp2-C!cx=O > 5 * * * 5
			if ((! a1->properties.get_bit(Atom::cx)) && bt2 == 2 &&
				a2->properties.get_bit(Atom::sp2) && bt3 == 1 &&
				(! a3->properties.get_bit(Atom::cx))) {
				a0->pc_class = Atom::POLAR;
				a4->pc_class = Atom::POLAR;
			}
			// O&x1-C&ar6~ar6~C&ar6=O&x1 > 5 * * * 5
			if (a1->properties.get_bit(Atom::ar6) &&
				a2->properties.get_bit(Atom::ar6) &&
				a3->properties.get_bit(Atom::ar6) &&
				a4->properties.get_bit(Atom::x1)) {
				a0->pc_class = Atom::POLAR;
				a4->pc_class = Atom::POLAR;
			}
		}

		get_reachable_atoms_from(a0, 6, paths); // paths to 6-bond reachable atoms
		for (vector<vector<Atom*> >::iterator
			p = paths.begin(); p != paths.end(); ++p) {
			vector<Atom*>& path = *p;
			// O&x1-C&ar6~ar6~ar6~ar6~C&ar6=O&x1 > 5 * * * * * 5

			a6 = path[6]; // 6-bond reachable atom
			if (! (a6->element == Atom::get_element_number("O") &&
				a6->properties.get_bit(Atom::x1))) continue;
			a5 = path[5]; // 5-bond reachable atom
			bt6 = bond_between(a5, a6);
			if (bt6 != 2) continue; // =O&x1

			a1 = path[1]; // 1-bond reachable atom
			if (! (a1->element == Atom::get_element_number("C") &&
				a1->properties.get_bit(Atom::ar6))) continue;
			bt1 = bond_between(a0, a1);
			if (bt1 != 1) continue; // -C&ar6

			a2 = path[2]; // 2-bond reachable atom
			if (! a2->properties.get_bit(Atom::ar6)) continue;

			a3 = path[3]; // 3-bond reachable atom
			if (! a3->properties.get_bit(Atom::ar6)) continue;

			a4 = path[4]; // 4-bond reachable atom
			if (! a4->properties.get_bit(Atom::ar6)) continue;

			if (! (a5->element == Atom::get_element_number("C") &&
				a5->properties.get_bit(Atom::ar6))) continue;

			a0->pc_class = Atom::POLAR;
			a6->pc_class = Atom::POLAR;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_adjacent_keto_enols:
// special ionizations of adjacent keto/enols
// rule-56: O&x1-C!cx=C&x2,x3-C!cx(=O)-O,N-* > 2 7 6 7 4 * *
// rule-57: O=C!cx-C&x2,x3-C!cx(=O)-O,N-* > 2 7 6 7 4 * *
// rule-58: O&x1-C&sp2~C&sp2-C(=O)-C&sp2~C&sp2-O&x1 > 2 7 * 7 2 * 7 2
// rule-59: O&x1-C&sp2~C&sp2~C(-O&x1)~C&sp2-C=O > 2 7 * 7 2 * 7 2
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_adjacent_keto_enols()
{
	assign_atom_ake1();
	assign_atom_ake2();
}

// assign_atom_ake1:
// rule-56: O&x1-C!cx=C&x2,x3-C!cx(=O)-O,N-* > 2 7 6 7 4 * *
// rule-57: O   =C!cx-C&x2,x3-C!cx(=O)-O,N-* > 2 7 6 7 4 * *
void Molecule::assign_atom_ake1()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	Atom* a3;
	Atom* a4;
	Atom* a5;
	vector<Atom*> bonds_a3;
	vector<vector<Atom*> > paths;
	unsigned int bt1, bt2, bt3, bt4, bt5, bt6;
	bool has_dbo; // if a3 has double-bonded O or not

	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (a0->element != Atom::get_element_number("O")) continue;
		get_reachable_atoms_from(a0, 5, paths);
			// paths to 5-bond reachable atoms

		for (vector<vector<Atom*> >::iterator
			p = paths.begin(); p != paths.end(); ++p) {
			vector<Atom*>& path = *p;

			a1 = path[1];
			bt1 = bond_between(a0, a1);
			if (! (a1->element == Atom::get_element_number("C") &&
				(!a1->properties.get_bit(Atom::cx)))) continue;

			a2 = path[2];
			bt2 = bond_between(a1, a2);
			if (! (a2->element == Atom::get_element_number("C") &&
				(a2->properties.get_bit(Atom::x2) ||
				 a2->properties.get_bit(Atom::x3)))) continue;

			a3 = path[3];
			bt3 = bond_between(a2, a3);
			if (! (a3->element == Atom::get_element_number("C") &&
				(!a1->properties.get_bit(Atom::cx)))) continue;
			if (bt3 != 1) continue;

			a4 = path[4];
			bt4 = bond_between(a3, a4);
			if (! (a4->element == Atom::get_element_number("O") ||
				a4->element == Atom::get_element_number("N"))) continue;
			if (bt4 != 1) continue;

			a5 = path[5];
			bt5 = bond_between(a4, a5);
			if (bt5 != 1) continue;

			get_atoms_bonded_to_atom(a3, bonds_a3);
			has_dbo = false; // if a3 has double-bonded O or not
			for (int i = 0; i < bonds_a3.size(); ++i) {
				if (bonds_a3[i]->element == Atom::get_element_number("O") &&
					m_bonds_to_atom[a3][i]->v == 2) {
					has_dbo = true;
					break;
				}
			}
			if (! has_dbo) continue;
			if ((a0->properties.get_bit(Atom::x1) &&
				bt1 == 1 && bt2 == 2) || (bt1 == 2 && bt2 == 1)) {
				a0->pc_class = Atom::ANION;
				a1->pc_class = Atom::NONE;
				a2->pc_class = Atom::HYDROPHOBIC;
				a3->pc_class = Atom::NONE;
				a4->pc_class = Atom::ACCEPTOR;
			}
		}
	}
}

// assign_atom_ake2:
// rule-58: O&x1 - C&sp2 ~ C&sp2 - C(=O)    -C&sp2 ~ C&sp2 -O&x1 > 2 7 * 7 2 * 7 2
// rule-59: O&x1 - C&sp2 ~ C&sp2 ~ C(-O&x1) ~C&sp2 - C     =O    > 2 7 * 7 2 * 7 2
void Molecule::assign_atom_ake2()
{
	Atom* a0;
	Atom* a1;
	Atom* a2;
	Atom* a3;
	Atom* a4;
	Atom* a5;
	Atom* a6;
	Atom* a7;
	vector<Atom*> bonds_a3;
	vector<Atom*> a3_dbo; // O double-bonded to a3
	vector<Atom*> a3_sbox1; // O&x1 single-bonded to a3
	vector<vector<Atom*> > paths;
	unsigned int bt1, bt2, bt3, bt4, bt5, bt6;

	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (! (a0->element != Atom::get_element_number("O") &&
				a0->properties.get_bit(Atom::x1))) continue;

		get_reachable_atoms_from(a0, 6, paths);
			// paths to 6-bond reachable atoms
		for (vector<vector<Atom*> >::iterator
			p = paths.begin(); p != paths.end(); ++p) {
			vector<Atom*>& path = *p;

			a1 = path[1];
			bt1 = bond_between(a0, a1);
			if (! (a1->element == Atom::get_element_number("C") &&
					a1->properties.get_bit(Atom::sp2))) continue;
			if (bt1 != 1) continue;

			a2 = path[2];
			bt2 = bond_between(a1, a2);
			if (! (a2->element == Atom::get_element_number("C") &&
					a2->properties.get_bit(Atom::sp2))) continue;

			a3 = path[3];
			bt3 = bond_between(a2, a3);
			if (a3->element != Atom::get_element_number("C")) continue;

			a4 = path[4];
			bt4 = bond_between(a3, a4);
			if (! (a4->element == Atom::get_element_number("C") &&
					a4->properties.get_bit(Atom::sp2))) continue;

			a5 = path[5];
			bt5 = bond_between(a4, a5);
			if (a5->element != Atom::get_element_number("C")) continue;

			a6 = path[6];
			bt6 = bond_between(a5, a6);
			if (a6->element != Atom::get_element_number("O")) continue;

			get_atoms_bonded_to_atom(a3, bonds_a3);
			//const vector<unsigned int>& bonds_type_a3 = m_types_of_bonds_to_atom[a3];
			a3_dbo.clear();
			a3_sbox1.clear();
			for (int i = 0; i < bonds_a3.size(); ++i) {
				a7 = bonds_a3[i];
				if (a7->element != Atom::get_element_number("O")) continue;
				if (m_bonds_to_atom[a3][i]->v == 2) a3_dbo.push_back(a7);
				if (m_bonds_to_atom[a3][i]->v == 1 &&
					a7->properties.get_bit(Atom::x1))
					a3_sbox1.push_back(a7);
			}

			if (bt3 == 1 && (! a3_dbo.empty()) && bt4 == 1 &&
				a5->properties.get_bit(Atom::sp2) && bt6 == 1 &&
				a6->properties.get_bit(Atom::x1)) {
				a0->pc_class = Atom::ANION;
				a1->pc_class = Atom::NONE;
				a3->pc_class = Atom::NONE;
				for (vector<Atom*>::iterator
					q = a3_dbo.begin(); q != a3_dbo.end(); ++q) {
					(*q)->pc_class = Atom::ANION;
				}
				a5->pc_class = Atom::NONE;
				a6->pc_class = Atom::ANION;
			} else if ((! a3_sbox1.empty()) && bt5 == 1 && bt6 == 2) {
				a0->pc_class = Atom::ANION;
				a1->pc_class = Atom::NONE;
				a3->pc_class = Atom::NONE;
				for (vector<Atom*>::iterator
					q = a3_sbox1.begin(); q != a3_sbox1.end(); ++q) {
					(*q)->pc_class = Atom::ANION;
				}
				a5->pc_class = Atom::NONE;
				a6->pc_class = Atom::ANION;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// assign_atom_misc_buried:
// miscellaneous fixups for "buried" atoms
// rule-60: C&x3&sp2~neg > 7 * (a flat carbon surrounded by > 1 neg atoms)
// rule-61: C,S&x2&sp~neg > 7 * (an sp carbon or sulfur surrounded by > 1 neg atoms)
// rule-62: N&x2&sp > 7 (an sp nitrogen surrounded by neighbors)
//////////////////////////////////////////////////////////////////////////////
void Molecule::assign_atom_misc_buried()
{
	Atom* a0;
	vector<Atom*> as;
	bool has_neg; // if it has an 'neg' atom bonded or not

	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		has_neg = false;
		get_atoms_bonded_to_atom(a0, as);
		for (vector<Atom*>::iterator p = as.begin(); p != as.end(); ++p) {
			if ((*p)->properties.get_bit(Atom::neg)) {
				has_neg = true;
				break;
			}
		}
		if (! has_neg) continue;

		if ((a0->element == Atom::get_element_number("C") &&
				a0->properties.get_bit(Atom::x3) &&
				a0->properties.get_bit(Atom::sp2)) ||
			((a0->element == Atom::get_element_number("C") ||
			  a0->element == Atom::get_element_number("S")) &&
				a0->properties.get_bit(Atom::x2) &&
				a0->properties.get_bit(Atom::sp))) {
			a0->pc_class = Atom::NONE;
		}
	}

	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p;
		if (a0->element == Atom::get_element_number("N") &&
			a0->properties.get_bit(Atom::x2) &&
			a0->properties.get_bit(Atom::sp)) {
			a0->pc_class = Atom::NONE;
		}
	}
}

// calc_atom_asas:
void Molecule::calc_atom_asas()
{
	calc_atom_asas(
		AsaCalculator::SOLVENT_RADIUS, AsaCalculator::ATOM_SLICING_STEP);
}

// calc_atom_asas:
void Molecule::calc_atom_asas(double rsolv, double zstep)
{
	if (natoms() == 0) return;
	vector<Atom*> atoms;
	for (int i = 0; i < natoms(); i++)
		atoms.push_back(atom(i));

	AsaCalculator ac;
	ac.calc(atoms, rsolv, zstep);
}

/*
// superpose:
void Molecule::superpose(const Superposer& s)
{
	for (int i = 0; i < m_atoms.size(); i++)
		s.transform(m_atoms[i].pos);
}
*/

//
void Molecule::print(std::ostream& os) const
{
	print_id(os);
	print_name(os);
	print_type(os);
	print_atoms(os);
	print_bonds(os);
}

//
void Molecule::print_entry_terminator(std::ostream& os) const
{
	os << "/" << std::endl;
}

//
void Molecule::print_id(std::ostream& os) const
{
	os << "I\t" << id() << std::endl;
}

//
void Molecule::print_name(std::ostream& os) const
{
	os << "N\t" << name() << std::endl;
}

//
void Molecule::print_type(std::ostream& os) const
{
	os << "T\t";
	switch (type()) {
	case PROTEIN:
		os << "protein"; break;
	case PEPTIDE:
		os << "peptide"; break;
	case DNA:
		os << "dna"; break;
	case RNA:
		os << "rna"; break;
	case SMALL:
		os << "small"; break;
	case WATER:
		os << "water"; break;
	case ATOM:
		os << "atom"; break;
	default:
		os << "undefined"; break;
	}
	os << std::endl;
}

//
void Molecule::print_atoms(std::ostream& os) const
{
	for (int i = 0; i < natoms(); i++)
		print_atom(os, i);
}

//
void Molecule::print_atom(std::ostream& os, size_t i) const
{
	using std::setw;
	using std::fixed;
	using std::left;
	using std::right;
	using std::setprecision;
	os << "A";

	const Atom* a = atom(i);
	os << "\t" << setw(4) << left << i;
	os << "\t" << setw(4) << left << a->number;
	os << "\t" << a->name;
	os << "\t" << Atom::symbols[a->element];
	os << "\t" << fixed << right << setprecision(4) << setw(10) << a->pos.x();
	os << "\t" << fixed << right << setprecision(4) << setw(10) << a->pos.y();
	os << "\t" << fixed << right << setprecision(4) << setw(10) << a->pos.z();
	os << "\t" << fixed << right << setprecision(4) << setw(10) << a->asa;
	os << "\t" << fixed << right << setprecision(2) << setw(4) << a->vdwr();
	os << "\t" << (int)a->pc_class;
	os << "\t" << (int)a->type();
	os << "\t" << a->properties.str();
	os << std::endl;
}

//
void Molecule::print_bonds(std::ostream& os) const
{
	for (int i = 0; i < nbonds(); i++)
		print_bond(os, i);
}

//
void Molecule::print_bond(std::ostream& os, size_t i) const
{
	const Bond* b = m_bonds[i];

	os << "B";
	os << "\t" << std::setw(4) << std::left << atom_number_of(b->atom1);
	os << "\t" << std::setw(4) << std::left << atom_number_of(b->atom2);
	os << "\t" << b->v;
	os << std::endl;
}

//
void Molecule::print_in_sdf(std::ostream& os, bool annot)
{
	os << id() << std::endl;
	print_in_sdf_header_line2(os);
	os << name() << std::endl;

	if (annot) assign_atom_types();
	print_in_sdf_count_line(os);

	for (int i = 0; i < natoms(); i++) {
		print_in_sdf_atom(os, atom(i));
	}

	for (int i = 0; i < m_bonds.size(); i++) {
		print_in_sdf_bond(os,
			atom_number_of(m_bonds[i]->atom1)+1,
			atom_number_of(m_bonds[i]->atom2)+1,
			m_bonds[i]->v);
	}

	os << "M  END" << std::endl;
	if (annot) print_in_sdf_annot(os);
	os << "$$$$" << std::endl;
}

//
void Molecule::print_in_sdf_header_line2(std::ostream& os) const
{
	os << "YMcba_mol 04180500003D" << std::endl;
}

//
void Molecule::print_in_sdf_count_line(std::ostream& os) const
{
	os << std::setw(3) << std::right << natoms();
	os << std::setw(3) << std::right << nbonds() << std::endl;
}

//
void Molecule::print_in_sdf_atom(std::ostream& os, const Atom* atom) const
{
	os << std::fixed << std::right << std::setprecision(4) << std::setw(10) << atom->pos.x();
	os << std::fixed << std::right << std::setprecision(4) << std::setw(10) << atom->pos.y();
	os << std::fixed << std::right << std::setprecision(4) << std::setw(10) << atom->pos.z();
	os << ' ';
	os << std::fixed << std::left << std::setw(3) << Atom::symbols[atom->element];
	os << " 0";
	os << std::setw(3) << std::right << atom->charge;
	for (int i = 0; i < 10; i++) os << "  0";
	os << std::endl;
}

//
void Molecule::print_in_sdf_bond(std::ostream& os, int i, int j, int t) const
{
	os << std::fixed << std::right << std::setw(3) << i;
	os << std::fixed << std::right << std::setw(3) << j;
	os << std::fixed << std::right << std::setw(3) << t;
	for (int i = 0; i < 4; i++) os << "  0";
	os << std::endl;
}

//
void Molecule::print_in_sdf_annot(std::ostream& os)
{
	assign_atom_types();
	os << "> <atom.names.classes>" << std::endl;
	for (int i = 0; i < natoms(); i++) {
		os << std::setw(5) << std::left << atom(i)->name << ' ';
		os << (int)(atom(i)->pc_class) << ' ';
		os << atom(i)->properties.str() << std::endl;
	}
	os << std::endl;

	find_rings();
	os << "> <rings>" << std::endl;
	os << m_rings.size() << std::endl;
	for (int i = 0; i < m_rings.size(); i++) {
		const vector<int>& ring = m_rings[i];
		for (int j = 0; j < ring.size(); j++) {
			os << ring[j]+1;
			if (j == ring.size()-1) os << std::endl;
			else os << ' ';
		}
	}
	os << std::endl;

	find_rotatable_bonds();
	os << "> <rotatable.bonds>" << std::endl;
	os << m_rotatable_bonds.size() << std::endl;
	for (int i = 0; i < m_rotatable_bonds.size(); i++) {
		os << std::setw(3) << std::left << m_rotatable_bonds[i]+1 << ' ';
		os << std::setw(3) << std::right << atom_number_of(m_bonds[m_rotatable_bonds[i]]->atom1)+1 << ' ';
		os << std::setw(3) << std::right << atom_number_of(m_bonds[m_rotatable_bonds[i]]->atom2)+1 << std::endl;
	}
	os << std::endl;

	os << "> <fragments>" << std::endl;
	int n = find_fragments();
	os << n << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m_fragments[i].size(); j++) {
			os << m_fragments[i][j]+1;
			if (j == m_fragments[i].size()-1) os << std::endl;
			else os << ' ';
		}
	}
	os << std::endl;
}

//
int Molecule::find_rings()
{
	m_rings.clear();
	vector<vector<Atom*> > rings_work;
	vector<int> ring;
	for (int i = 3; i <= 8; i++) {
		get_n_membered_rings(rings_work, i);
		for (int j = 0; j < rings_work.size(); j++) {
			ring.clear();
			for (int k = 0; k < rings_work[j].size(); k++) {
				ring.push_back(atom_number_of(rings_work[j][k]));
			}
			m_rings.push_back(ring);
		}
	}
	return m_rings.size();
}

//
int Molecule::find_rotatable_bonds()
{
	if (! m_rotatable_bonds.empty()) return m_rotatable_bonds.size();
	if (! m_atom_types_assigned) assign_atom_types();
	if (m_rings.empty()) find_rings();
	m_rotatable_bonds.clear();

	// only single bonds can be rotatable
	for (int i = 0; i < m_bonds.size(); i++) {
		if (m_bonds[i]->v == 1) m_rotatable_bonds.push_back(i);
	}

	// 'terminal' bonds cannot be rotatable
	Bond* b;
	for (vector<int>::iterator p = m_rotatable_bonds.begin(); p != m_rotatable_bonds.end(); ) {
		b = m_bonds[*p];
		if (b->atom1->properties.get_bit(Atom::x1) == 1 ||
			b->atom2->properties.get_bit(Atom::x1) == 1) {
			p = m_rotatable_bonds.erase(p);
		} else {
			++p;
		}
	}

	// ring bonds cannot be rotatable
	for (vector<int>::iterator p = m_rotatable_bonds.begin(); p != m_rotatable_bonds.end(); ) {
		if (is_ring_bond(*p)) {
			p = m_rotatable_bonds.erase(p);
		} else {
			++p;
		}
	}

	// amide bonds cannot be rotatable
	Atom* c;
	Atom* n;
	Atom* o;
	for (vector<int>::iterator p = m_rotatable_bonds.begin(); p != m_rotatable_bonds.end(); ) {
		b = m_bonds[*p];
		c = 0; n = 0; o = 0;

		if (b->atom1->element == Atom::get_element_number("C") &&
			b->atom2->element == Atom::get_element_number("N")) {
			c = b->atom1;
			n = b->atom2;
		} else if (b->atom1->element == Atom::get_element_number("N") &&
			b->atom2->element == Atom::get_element_number("C")) {
			c = b->atom2;
			n = b->atom1;
		}

		if (c && n && bond_between(c, n) == 1) {
			vector<Atom*> as;
			get_atoms_bonded_to(c, as);
			for (vector<Atom*>::iterator q = as.begin(); q != as.end(); ++q) {
				if ((*q)->element == Atom::get_element_number("O") && bond_between(c, *q) == 2) {
					o = *q;
					break;
				}
			}
		}

		if (o) {
			p = m_rotatable_bonds.erase(p);
		} else {
			++p;
		}
	}

	return m_rotatable_bonds.size();
}

//
bool Molecule::is_ring_bond(int b) const
{
	int nmem;
	for (vector<vector<int> >::const_iterator p = m_rings.begin(); p != m_rings.end(); ++p) {
		const vector<int>& ring = *p;
		nmem = ring.size();
		for (int i = 0; i < nmem; i++) {
			int j = i + 1;
			if (j == nmem) j = 0;
			if ((atom(ring[i]) == m_bonds[b]->atom1 && atom(ring[j]) == m_bonds[b]->atom2) ||
				(atom(ring[i]) == m_bonds[b]->atom2 && atom(ring[j]) == m_bonds[b]->atom1))
				return true;
		}
	}
	return false;
}

//
bool Molecule::is_ring_bond(const Bond* b) const
{
	int nmem;
	for (vector<vector<int> >::const_iterator p = m_rings.begin(); p != m_rings.end(); ++p) {
		const vector<int>& ring = *p;
		nmem = ring.size();
		for (int i = 0; i < nmem; i++) {
			int j = i + 1;
			if (j == nmem) j = 0;
			if ((atom(ring[i]) == b->atom1 && atom(ring[j]) == b->atom2) ||
				(atom(ring[i]) == b->atom2 && atom(ring[j]) == b->atom1))
				return true;
		}
	}
	return false;
}

//
int Molecule::get_rings(vector<vector<int> >& rings)
{
	if (m_rings.empty()) find_rings();
	rings = m_rings;
	return rings.size();
}

//
int Molecule::get_rings(vector<vector<Atom*> >& rings)
{
	rings.clear();
	if (m_rings.empty()) find_rings();
	vector<Atom*> ring;
	for (int i = 0; i < m_rings.size(); i++) {
		ring.clear();
		for (int j = 0; j < m_rings[i].size(); j++) {
			ring.push_back(atom(m_rings[i][j]));
		}
		rings.push_back(ring);
	}
	return rings.size();
}

//
int Molecule::get_n_membered_rings(vector<vector<Atom*> >& rings, int n)
{
	rings.clear();
	if (n < 3) return 0;
	map<Atom*, bool> in_ring;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p)
		in_ring[*p] = false;

	Atom* a0;
	vector<vector<Atom*> > paths;
	vector<Atom*> atoms_bonded;
	for (vector<Atom*>::iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a0 = *p; // starting atom
		if (in_ring[a0]) continue;
		paths.clear();
		get_reachable_atoms_from(a0, n-1, paths);
		for (int i = 0; i < paths.size(); i++) {
			vector<Atom*>& path = paths[i];
			get_atoms_bonded_to(path.at(n-1), atoms_bonded);
			for (int j = 0; j < atoms_bonded.size(); j++) {
				if (atoms_bonded[j] == a0) { // n-membered ring found
					rings.push_back(path);
					for (int k = 0; k < path.size(); k++)
						in_ring[path[k]] = true;
					break;
				}
			}
			if (in_ring[a0]) break;
		}
	}

	return rings.size();
}

//
int Molecule::get_aromatic_rings(vector<vector<int> >& rings)
{
	rings.clear();
	if (m_rings.empty()) find_rings();
	assign_atom_types();

	//bool is_ar_ring;
	int nars;
	for (int i = 0; i < m_rings.size(); i++) {
		const vector<int>& ring = m_rings[i];
		//is_ar_ring = true;
		nars = 0;
		for (int j = 0; j < ring.size(); j++) {
			if (atom(ring[j])->properties.get_bit(Atom::ar) == 1) {
				nars++;
				if (nars > 2) break;
			}
			//if (atom(ring[j])->properties.get_bit(Atom::ar) != 1) {
			//	is_ar_ring = false;
			//	break;
			//}
		}
		//if (is_ar_ring) {
		if (nars > 2) {
			rings.push_back(ring);
		}
	}

	return rings.size();
}

//
int Molecule::get_aromatic_rings(vector<vector<Atom*> >& rings)
{
	rings.clear();
	vector<vector<int> > rings_work;
	get_aromatic_rings(rings_work);

	vector<Atom*> ring;
	for (int i = 0; i < rings_work.size(); i++) {
		ring.clear();
		for (int j = 0; j < rings_work[i].size(); j++) {
			ring.push_back(atom(rings_work[i][j]));
		}
		rings.push_back(ring);
	}

	return rings.size();
}

//
int Molecule::find_fragments()
{
	if (! m_fragments.empty()) return m_fragments.size();
	m_fragments.clear();
	vector<int> fragment;
	find_rotatable_bonds();

	vector<int> atoms;
	for (int i = 0; i < m_atoms.size(); i++) atoms.push_back(i);
	vector<int> bonds;
	if (m_rotatable_bonds.empty()) {
		for (int j = 0; j < m_bonds.size(); j++) bonds.push_back(j);
	} else {
		for (int j = 0; j < m_rotatable_bonds[0]; j++) bonds.push_back(j);
		for (int i = 1; i < m_rotatable_bonds.size(); i++) {
			for (int j = m_rotatable_bonds[i-1]+1; j < m_rotatable_bonds[i]; j++) bonds.push_back(j);
		}
		for (int j = m_rotatable_bonds.back()+1; j < m_bonds.size(); j++) bonds.push_back(j);
	}

	queue<int> atom_queue;
	vector<int>::iterator p;
	while (! atoms.empty()) {
		fragment.clear();
		p = atoms.begin();
		atom_queue.push(*p);
		atoms.erase(p);
		while (! atom_queue.empty()) {
			for (p = atoms.begin(); p != atoms.end(); ) {
				int i;
				Bond* b;
				for (i = 0; i < bonds.size(); i++) {
					b = m_bonds[bonds[i]];
					if ((b->atom1 == m_atoms[atom_queue.front()] && b->atom2 == m_atoms[*p]) ||
						(b->atom2 == m_atoms[atom_queue.front()] && b->atom1 == m_atoms[*p])) {
						atom_queue.push(*p);
						break;
					}
				}
				if (i < bonds.size()) {
					p = atoms.erase(p);
				} else {
					++p;
				}
			}
			fragment.push_back(atom_queue.front());
			atom_queue.pop();
		}
		m_fragments.push_back(fragment);
	}

	return m_fragments.size();
}

//
int Molecule::get_fragments(vector<vector<int> >& fragments)
{
	find_fragments();
	fragments = m_fragments;
	return fragments.size();
}

//
int Molecule::get_fragments(vector<vector<Atom*> >& fragments)
{
	fragments.clear();
	find_fragments();
	vector<Atom*> fragment;
	for (int i = 0; i < m_fragments.size(); i++) {
		fragment.clear();
		for (int j = 0; j < m_fragments[i].size(); j++) {
			fragment.push_back(atom(m_fragments[i][j]));
		}
		fragments.push_back(fragment);
	}
	return fragments.size();
}

//
void Molecule::get_bond(size_t i, BondBetweenAtoms& bond) const
{
	if (i >= m_bonds.size()) {
		bond.atom1 = 0;
		bond.atom2 = 0;
		bond.v = 0;
	} else {
		bond.atom1 = atom_number_of(m_bonds[i]->atom1);
		bond.atom2 = atom_number_of(m_bonds[i]->atom2);
		bond.v = m_bonds[i]->v;
	}
}

//
int Molecule::get_bonds(vector<BondBetweenAtoms>& bonds) const
{
	bonds.clear();
	BondBetweenAtoms b;
	for (int i = 0; i < m_bonds.size(); i++) {
		get_bond(i, b);
		bonds.push_back(b);
	}
	return bonds.size();
}

//
int Molecule::get_rotatable_bonds(std::vector<int>& rotatable_bonds)
{
	rotatable_bonds.clear();
	find_rotatable_bonds();
	rotatable_bonds = m_rotatable_bonds;
	return rotatable_bonds.size();
}

//
void Molecule::get_fragment(size_t i, vector<int>& fragment)
{
	fragment.clear();
	if (i < find_fragments()) {
		fragment = m_fragments[i];
	}
}

//
void Molecule::get_fragment(size_t i, vector<Atom*>& fragment)
{
	fragment.clear();
	if (i < find_fragments()) {
		for (int k = 0; k < m_fragments[i].size(); k++)
			fragment.push_back(atom(m_fragments[i][k]));
	}
}

//
void Molecule::get_fragment(size_t i, Molecule& fragment)
{
	fragment.clear();
	if (i >= find_fragments()) return;
	set_fragment_id(i, fragment);
	set_fragment_name(i, fragment);
	fragment.set_type(Molecule::SMALL);
	set_fragment_atoms(m_fragments[i], fragment);
	set_fragment_bonds(m_fragments[i], fragment);
	fragment.m_atom_types_assigned = true;
}

//
void Molecule::set_fragment_id(int i, Molecule& fragment)
{
	std::ostringstream oss;
	oss << id();
	oss << '_' << i;
	fragment.set_id(oss.str());
}

//
void Molecule::set_fragment_name(int i, Molecule& fragment)
{
	std::ostringstream oss;
	oss << "fragment " << i << " of " << name();
	fragment.set_name(oss.str());
}

//
void Molecule::set_fragment_atoms(const vector<int>& atoms, Molecule& fragment)
{
	Atom* a;
	for (int i = 0; i < atoms.size(); i++) {
		a = atom(atoms[i]);
		fragment.add_atom(*a);
	}
}

//
void Molecule::set_fragment_bonds(const vector<int>& atoms, Molecule& fragment)
{
	int v;
	for (int i = 0; i < atoms.size()-1; i++) {
		for (int j = i+1; j < atoms.size(); j++) {
			if (v = bond_between(atoms[i], atoms[j])) {
				fragment.add_bond(fragment.atom(i), fragment.atom(j), v);
			}
		}
	}
}

//
void Molecule::calc_graph_distances()
{
	if (natoms() == m_graph_distances.size()) return;
	initialize_graph_distances();

	vector<vector<int> > linked_atoms; // atoms linked to i-th atom
	get_linked_atoms(linked_atoms);

	for (int i = 0; i < natoms(); i++) {
		calc_graph_distances_from(i, linked_atoms);
	}
}

//
void Molecule::initialize_graph_distances()
{
	m_graph_distances.clear();
	m_graph_distances.resize(natoms());
	for (int i = 0; i < natoms(); i++) {
		m_graph_distances[i].resize(natoms());
		for (int j = 0; j < natoms(); j++) {
			m_graph_distances[i][j] = -1;
		}
	}
}

//
void Molecule::get_linked_atoms(vector<vector<int> >& linked_atoms)
{
	map<const Atom*, int> atom_numbers;
	for (int i = 0; i < natoms(); i++) {
		atom_numbers[atom(i)] = atom_number_of(atom(i));
	}

	linked_atoms.clear();
	linked_atoms.resize(natoms());
	vector<Atom*> atoms;
	for (int i = 0; i < natoms(); i++) {
		atoms.clear();
		get_atoms_bonded_to(atom(i), atoms);
		for (int j = 0; j < atoms.size(); j++) {
			linked_atoms[i].push_back(atom_numbers[atoms[j]]);
		}
	}
}

//
void Molecule::calc_graph_distances_from(int atomi, const vector<vector<int> >& linked_atoms)
{
	if (atomi < 0 || natoms() <= atomi) return;

	/////// breadth-first saerch //////
	queue<int> q;
	q.push(atomi);
	vector<int>& dist = m_graph_distances[atomi];
	dist[atomi] = 0;
	vector<int> p;
	for (int i = 0; i < natoms(); i++) p.push_back(0);
	p[atomi] = 1;

	while (! q.empty()) {
		int a = q.front();
		const vector<int>& as = linked_atoms[a];
		for (int i = 0; i < as.size(); i++) {
			int b = as[i];
			if (p[b] == 0) { // if not yet reached in the search path
				q.push(b);
				p[b] = 1;
				dist[b] = dist[a] + 1;
			}
		}
		q.pop();
	}
}

//
int Molecule::get_graph_distance(int a1, int a2)
{
	if (a1 < 0 || natoms() <= a1) return -1;
	if (a2 < 0 || natoms() <= a2) return -1;
	calc_graph_distances();
	return m_graph_distances[a1][a2];
}

//
int Molecule::get_graph_distance(const Atom* a1, const Atom* a2)
{
	return get_graph_distance(atom_number_of(a1), atom_number_of(a2));
}

//
bool Molecule::is_a_complex()
{
	for (int i = 1; i < natoms(); i++) {
		if (get_graph_distance(0,i) < 0) return true;
	}
	return false;
}

//
int Molecule::get_component_molecules(vector<Molecule>& mols)
{
	mols.clear();
	list<int> atoms;
	for (int i = 0; i < natoms(); i++) atoms.push_back(i);

	Molecule mol;
	vector<int> atoms_connected;
	while (! atoms.empty()) {
		mol.clear();
		mol.set_id(id());
		mol.set_name(name());
		mol.set_type(type());
		atoms_connected.clear();

		// add atoms
		int a0 = atoms.front();
		mol.add_atom(*atom(a0));
		atoms.pop_front();
		atoms_connected.push_back(a0);
		for (list<int>::iterator li = atoms.begin(); li != atoms.end(); ) {
			if (get_graph_distance(a0,*li) > 0) {
				mol.add_atom(*atom(*li));
				atoms_connected.push_back(*li);
				li = atoms.erase(li);
			} else {
				++li;
			}
		}
		mol.m_atom_types_assigned = m_atom_types_assigned;

		// add bonds
		for (int i = 0; i < (atoms_connected.size()-1); i++) {
			for (int j = i+1; j < atoms_connected.size(); j++) {
				int k = bond_between(atoms_connected[i],atoms_connected[j]);
				if (k > 0) {
					mol.add_bond(mol.atom(i), mol.atom(j), k);
				}
			}
		}

		mols.push_back(mol);
	}

	return mols.size();
}

//
bool Molecule::is_identical_to(Molecule& mol)
{
	vector<vector<int> > a;
	return is_identical_to(mol, a);
}

//
bool Molecule::is_identical_to(Molecule& mol, vector<vector<int> >& alignments)
{
	alignments.clear();
	if (natoms() != mol.natoms()) return false;
	if (nbonds() != mol.nbonds()) return false;
	if (! atom_counts_equal(mol)) return false;
	if (! bond_counts_equal(mol)) return false;

	//assign_atom_types();
	//mol.assign_atom_types();
	//calc_graph_distances();
	//mol.calc_graph_distances();

	vector<vector<int> > equivalent_atom_candidates;
		// [i][j]: for i-th atom of this molecule, j-th candidate from mol
	get_equivalent_atom_candidates(mol, equivalent_atom_candidates);
	for (int i = 0; i < natoms(); i++) {
		if (equivalent_atom_candidates[i].empty()) return false;
	}

	vector<int> alignment;
	alignment.resize(natoms());
	search_equivalent_atom_alignments(0, mol, equivalent_atom_candidates, alignment, alignments);
	if (alignments.empty()) return false;
	return true;
}

//
bool Molecule::atom_counts_equal(const Molecule& mol) const
{
	int nelements = 256;
	vector<int> c1(nelements, 0);
	vector<int> c2(nelements, 0);
	for (int i = 0; i < m_atoms.size(); i++) {
		if (m_atoms[i]->element < nelements) {
			c1[m_atoms[i]->element]++;
		}
	}
	for (int i = 0; i < mol.m_atoms.size(); i++) {
		if (mol.m_atoms[i]->element < nelements) {
			c2[mol.m_atoms[i]->element]++;
		}
	}
	for (int i = 0; i < nelements; i++) {
		if (c1[i] != c2[i]) return false;
	}
	return true;
}

//
bool Molecule::bond_counts_equal(const Molecule& mol) const
{
	int nbondtypes = 4;
	vector<int> c1(nbondtypes, 0);
	vector<int> c2(nbondtypes, 0);
	for (int i = 0; i < m_bonds.size(); i++) {
		if (m_bonds[i]->v < nbondtypes) {
			c1[m_bonds[i]->v]++;
		}
	}
	for (int i = 0; i < mol.m_bonds.size(); i++) {
		if (mol.m_bonds[i]->v < nbondtypes) {
			c2[mol.m_bonds[i]->v]++;
		}
	}
	if (c1 != c2) return false;
	return true;
}

//
void Molecule::get_equivalent_atom_candidates(Molecule& mol, vector<vector<int> >& equiv_candidates)
{
	equiv_candidates.clear();
	equiv_candidates.resize(natoms());
	int na = m_atoms.size();
	if (na != mol.m_atoms.size()) return;
	int nb = m_bonds.size();
	if (nb != mol.m_bonds.size()) return;

	vector<vector<int> > bondtypes1(na);
	vector<vector<int> > bondtypes2(na);
	int nbondtypes = 4;
	for (int i = 0; i < na; i++) {
		bondtypes1[i].resize(nbondtypes);
		bondtypes2[i].resize(nbondtypes);
		for (int j = 0; j < nbondtypes; j++) bondtypes1[i][j] = bondtypes2[i][j] = 0;
	}
	map<const Atom*, int> an1;
	map<const Atom*, int> an2;
	for (int i = 0; i < na; i++) {
		an1[m_atoms[i]] = i;
		an2[mol.m_atoms[i]] = i;
	}
	map<const Atom*, int>::const_iterator mci;
	for (int i = 0; i < nb; i++) {
		if (0 < m_bonds[i]->v && m_bonds[i]->v < nbondtypes) {
			if ((mci = an1.find(m_bonds[i]->atom1)) != an1.end()) (bondtypes1[mci->second][m_bonds[i]->v])++;
			if ((mci = an1.find(m_bonds[i]->atom2)) != an1.end()) (bondtypes1[mci->second][m_bonds[i]->v])++;
		}
		if (0 < mol.m_bonds[i]->v && mol.m_bonds[i]->v < nbondtypes) {
			if ((mci = an2.find(mol.m_bonds[i]->atom1)) != an2.end()) (bondtypes2[mci->second][mol.m_bonds[i]->v])++;
			if ((mci = an2.find(mol.m_bonds[i]->atom2)) != an2.end()) (bondtypes2[mci->second][mol.m_bonds[i]->v])++;
		}
	}

	for (int i = 0; i < na; i++) {
		for (int j = 0; j < na; j++) {
			if (m_atoms[i]->element != mol.m_atoms[j]->element) continue;
			if (bondtypes1[i] != bondtypes2[j]) continue;
			equiv_candidates[i].push_back(j);
		}
	}

	/*
	vector<list<int> > gd1;
	for (int i = 0; i < natoms(); i++) {
		for (int k = 0; k < natoms(); k++) {
			if (k == i) continue;
			gd1[i].push_back(m_graph_distances[i][k]);
		}
		gd1[i].sort();
	}
	vector<list<int> > gd2;
	for (int i = 0; i < mol.natoms(); i++) {
		for (int k = 0; k < mol.natoms(); k++) {
			if (k == i) continue;
			gd2[i].push_back(mol.m_graph_distances[i][k]);
		}
		gd2[i].sort();
	}

	for (int i = 0; i < natoms(); i++) {
		for (int j = 0; j < mol.natoms(); j++) {
			if (atom(i)->element != mol.atom(j)->element) continue;
			if (atom(i)->properties != mol.atom(j)->properties) continue;
			if (atom(i)->pc_class != mol.atom(j)->pc_class) continue;
			if (gd1[i] != gd2[j]) continue;
			equiv_candidates[i].push_back(j);
		}
	}
	*/
}

//
void Molecule::search_equivalent_atom_alignments(int iatom,
	const Molecule& mol, const vector<vector<int> >& equivalent_candidates,
	vector<int>& alignment, vector<vector<int> >& alignments)
{
	if (iatom == natoms()) {
		alignments.push_back(alignment);
		return;
	}

	int na = natoms();
	if (na != mol.natoms()) return;
	vector<vector<int> > bonds1(na);
	vector<vector<int> > bonds2(na);
	for (int i = 0; i < na; i++) {
		bonds1[i].resize(na);
		bonds2[i].resize(na);
	}
	for (int i = 0; i < na; i++) {
		bonds1[i][i] = 0;
		bonds2[i][i] = 0;
		for (int j = (i+1); j < na; j++) {
			bonds1[i][j] = bonds1[j][i] = bond_between(m_atoms[i], m_atoms[j]);
			bonds2[i][j] = bonds2[j][i] = mol.bond_between(mol.m_atoms[i], mol.m_atoms[j]);
		}
	}

	int jatom, katom;
	for (int equiv_candidate = 0; equiv_candidate < equivalent_candidates[iatom].size(); equiv_candidate++) {
		jatom = equivalent_candidates[iatom][equiv_candidate];
		for (katom = 0; katom < iatom; katom++) {
			if ((alignment[katom] == jatom) ||
				(bonds1[katom][iatom] != bonds2[alignment[katom]][jatom]))
				break;
		}
		if (katom == iatom) {
			alignment[iatom] = jatom;
			search_equivalent_atom_alignments(iatom+1, mol, equivalent_candidates, alignment, alignments);
		}
	}
	/*
	int jatom, katom;
	for (int equiv_candidate = 0; equiv_candidate < equivalent_candidates[iatom].size(); equiv_candidate++) {
		jatom = equivalent_candidates[iatom][equiv_candidate];
		for (katom = 0; katom < iatom; katom++) {
			if ((alignment[katom] == jatom) ||
				(m_graph_distances[katom][iatom] != mol.m_graph_distances[alignment[katom]][jatom]))
				break;
		}
		if (katom == iatom) {
			alignment[iatom] = jatom;
			search_equivalent_atom_alignments(iatom+1, mol, equivalent_candidates, alignment, alignments);
		}
	}
	*/
}

//
void Molecule::get_morgan_numbering(vector<int>& ns)
{
	vector<vector<int> > hist;
	get_morgan_numbering(ns, hist);
}

//
void Molecule::get_morgan_numbering(vector<int>& ns, vector<vector<int> >& hist)
{
	int nas = natoms();
	vector<vector<int> > atom_neighbors(nas);
	for (int i = 0; i < nas; i++) {
		get_atoms_bonded_to(i, atom_neighbors[i]);
	}

	int nsteps = get_morgan_hist(hist, atom_neighbors);
	do_morgan_numbering(ns, nsteps, hist, atom_neighbors);
}

//
int Molecule::get_morgan_hist(vector<vector<int> >& hist, const vector<vector<int> >& atom_neighbors)
{
	int nas = natoms();
	hist.clear();
	hist.resize(2*nas);
	for (int i = 0; i < hist.size(); i++) hist[i].resize(nas+1);
	set<int> s;
	for (int i = 0; i < nas; i++) {
		hist[0][i] = atom_neighbors[i].size();
		s.insert(hist[0][i]);
	}
	hist[0][nas] = s.size();

	int nsteps;
	for (nsteps = 1; ; nsteps++) {
		s.clear();
		for (int i = 0; i < nas; i++) {
			hist[nsteps][i] = 0;
			for (int j = 0; j < atom_neighbors[i].size(); j++) {
				hist[nsteps][i] += hist[nsteps-1][atom_neighbors[i][j]];
			}
			s.insert(hist[nsteps][i]);
		}
		hist[nsteps][nas] = s.size();
		if (hist[nsteps][nas] <= hist[nsteps-1][nas]) break;
	}
	hist.resize(nsteps+1);
	nsteps--;
	return nsteps;
}

//
void Molecule::do_morgan_numbering(vector<int>& ns,
	int nsteps, vector<vector<int> >& hist, vector<vector<int> >& atom_neighbors)
{
	ns.clear();
	int nas = natoms();
	vector<int> atom_done(nas, 0);
	vector<int> atmp;
	for (int i = 0; i < nas; i++) atmp.push_back(i);
	sort_atoms_in_morgan(-1, atmp, hist[nsteps]);
	queue<int> q;
	q.push(atmp[0]);
	atom_done[q.front()] = 1;

	while (! q.empty()) {
		atmp.clear();
		const vector<int>& neighbors = atom_neighbors[q.front()];
		for (int i = 0; i < neighbors.size(); i++) {
			if (atom_done[neighbors[i]] == 0) {
				atmp.push_back(neighbors[i]);
			}
		}
		sort_atoms_in_morgan(q.front(), atmp, hist[nsteps]);
		for (int i = 0; i < atmp.size(); i++) {
			q.push(atmp[i]);
			atom_done[atmp[i]] = 1;
		}
		ns.push_back(q.front());
		q.pop();

		// in the case this contains two or more molecules
		if (q.empty()) {
			atmp.clear();
			for (int i = 0; i < nas; i++) {
				if (atom_done[i] == 0) atmp.push_back(i);
			}
			if (! atmp.empty()) {
				sort_atoms_in_morgan(-1, atmp, hist[nsteps]);
				q.push(atmp[0]);
				atom_done[q.front()] = 1;
			}
		}
	}
}

//
void Molecule::sort_atoms_in_morgan(int a, std::vector<int>& as, const std::vector<int>& ns)
{
	int nas = as.size();
	if (nas <= 1) return;
	int maxa, atmp, b1, b2, amax, aj;
	for (int i = 0; i < nas-1; i++) {
		maxa = i;
		for (int j = i+1; j < nas; j++) {
			aj = as[j];
			amax = as[maxa];
			if (ns[aj] > ns[amax]) {
				maxa = j;
			} else if (ns[aj] == ns[amax]) {
				if (a >= 0) {
					b1 = bond_between(a, aj);
					b2 = bond_between(a, amax);
				} else {
					b1 = b2 = 0;
				}
				if (b1 > b2) {
					maxa = j;
				} else if (b1 == b2) {
					const Atom* a1 = atom(aj);
					const Atom* a2 = atom(amax);
					if (a1->element < a2->element) {
						maxa = j;
					}
				}
			}
		}
		if (maxa != i) {
			atmp = as[i];
			as[i] = as[maxa];
			as[maxa] = atmp;
		}
	}
}
