#include <vector>
#include <map>
#include <string>
//#include <iostream>
//#include <fstream>
//#include <algorithm>
#include "pdb_small_mol_equiv.h"
#include "molecule.h"
#include "atom.h"

using namespace Cba;
using std::vector;
using std::map;
using std::string;

//
PdbSmallMolEquiv::PdbSmallMolEquiv() : m_mol1(0), m_mol2(0) {}

//
PdbSmallMolEquiv::~PdbSmallMolEquiv() {}

//
bool PdbSmallMolEquiv::set_equiv_mols(Molecule* mol1, Molecule* mol2)
{
	if (mol1->natoms() != mol2->natoms()) {
		m_mol1 = 0;
		m_mol2 = 0;
		return false;
	}

	m_mol1 = mol1;
	m_mol2 = mol2;
	if (! gather_atom_names(m_mol1, m_atom_names1)) return false;
	if (! gather_atom_names(m_mol2, m_atom_names2)) return false;
	if (! see_atom_name_equiv()) return false;
	return true;
}

//
bool PdbSmallMolEquiv::gather_atom_names(Molecule* mol, map<string, Atom*>& atom_names)
{
	atom_names.clear();
	map<string, Atom*>::iterator p;
	for (int i = 0; i < mol->natoms(); i++) {
		if ((p = atom_names.find(mol->atom(i)->name)) != atom_names.end())
			return false; // each atom must have a unique name
		atom_names[mol->atom(i)->name] = mol->atom(i);
	}
	return true;
}

//
bool PdbSmallMolEquiv::see_atom_name_equiv()
{
	m_atom1_to_atom2.clear();
	m_atom2_to_atom1.clear();
	map<string, Atom*>::iterator p = m_atom_names1.begin();
	map<string, Atom*>::iterator q = m_atom_names2.begin();
	for (; p != m_atom_names1.end() && q != m_atom_names2.end(); ++p, ++q) {
		if (p->first != q->first) return false;
		m_atom1_to_atom2[p->second] = q->second;
		m_atom2_to_atom1[q->second] = p->second;
	}
	if (p != m_atom_names1.end() || q != m_atom_names2.end()) return false;
	return true;
}

//
Atom* PdbSmallMolEquiv::get_equiv_atom(int i, Atom* atom)
{
	if (m_mol1 == 0 || m_mol2 == 0) return 0;
	map<Atom*, Atom*>::iterator p;
	if (i == 1) {
		if ((p = m_atom1_to_atom2.find(atom)) == m_atom1_to_atom2.end()) return 0;
	} else if (i == 2) {
		if ((p = m_atom2_to_atom1.find(atom)) == m_atom2_to_atom1.end()) return 0;
	} else {
		return 0;
	}
	return p->second;
}

//
Atom* PdbSmallMolEquiv::get_atom_with_name(int i, const std::string& name)
{
	if (m_mol1 == 0 || m_mol2 == 0) return 0;
	map<string, Atom*>::iterator p;
	if (i == 1) {
		if ((p = m_atom_names1.find(name)) == m_atom_names1.end()) return 0;
	} else if (i == 2) {
		if ((p = m_atom_names2.find(name)) == m_atom_names2.end()) return 0;
	} else {
		return 0;
	}
	return p->second;
}
