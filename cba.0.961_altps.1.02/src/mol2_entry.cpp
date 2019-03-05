#include <vector>
#include <string>
#include "mol2_entry.h"
#include "molecule.h"
#include "atom.h"

using namespace Cba;
using std::vector;
using std::string;

//
Mol2Entry::Mol2Entry() {}

//
Mol2Entry::~Mol2Entry() {}

//
void Mol2Entry::clear()
{
	m_mol.clear();
	m_atoms.clear();
	m_bonds.clear();
}

//
void Mol2Entry::MoleculeRecord::clear()
{
	name.erase();
	natoms = nbonds = 0;
	type.erase();
	charge_type.erase();
	status_bit.erase();
	comment.erase();
}

//
bool Mol2Entry::create_molecule(Molecule& mol) const
{
	mol.clear();
	mol.set_id(name());
	mol.set_name(name());
	mol.set_type(Molecule::SMALL);

	Atom a;
	for (vector<AtomRecord>::const_iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a.name = p->name;
		set_atom_element(a, p->type);
		a.pos.move_to(p->x, p->y, p->z);
		mol.add_atom(a);
	}

	int a1, a2;
	for (vector<BondRecord>::const_iterator p = m_bonds.begin(); p != m_bonds.end(); ++p) {
		if ((a1 = get_atom_serial_number(p->atom1)) < 0) continue;
		if ((a2 = get_atom_serial_number(p->atom2)) < 0) continue;
		mol.add_bond(mol.atom(a1), mol.atom(a2), get_bond_type(p->type));
	}
}

//
void Mol2Entry::set_atom_element(Atom& a, const string& atom_type) const
{
	// according to "Mol2 File Format" document (Sybyl 6.8.1, 5/2002, Tripos)
	if (atom_type == "C.3") a.element = Atom::get_element_number("C");
	else if (atom_type == "C.2") a.element = Atom::get_element_number("C");
	else if (atom_type == "C.1") a.element = Atom::get_element_number("C");
	else if (atom_type == "C.ar") a.element = Atom::get_element_number("C");
	else if (atom_type == "C.cat") a.element = Atom::get_element_number("C");
	else if (atom_type == "N.3") a.element = Atom::get_element_number("N");
	else if (atom_type == "N.2") a.element = Atom::get_element_number("N");
	else if (atom_type == "N.1") a.element = Atom::get_element_number("N");
	else if (atom_type == "N.ar") a.element = Atom::get_element_number("N");
	else if (atom_type == "N.am") a.element = Atom::get_element_number("N");
	else if (atom_type == "N.pl3") a.element = Atom::get_element_number("N");
	else if (atom_type == "N.4") a.element = Atom::get_element_number("N");
	else if (atom_type == "O.3") a.element = Atom::get_element_number("O");
	else if (atom_type == "O.2") a.element = Atom::get_element_number("O");
	else if (atom_type == "O.co2") a.element = Atom::get_element_number("O");
	else if (atom_type == "O.spc") a.element = Atom::get_element_number("O");
	else if (atom_type == "O.t3p") a.element = Atom::get_element_number("O");
	else if (atom_type == "S.3") a.element = Atom::get_element_number("S");
	else if (atom_type == "S.2") a.element = Atom::get_element_number("S");
	else if (atom_type == "S.O") a.element = Atom::get_element_number("S");
	else if (atom_type == "S.O2") a.element = Atom::get_element_number("S");
	else if (atom_type == "P.3") a.element = Atom::get_element_number("P");
	else if (atom_type == "H") a.element = Atom::get_element_number("H");
	else if (atom_type == "H.spc") a.element = Atom::get_element_number("H");
	else if (atom_type == "H.t3p") a.element = Atom::get_element_number("H");
	else if (atom_type == "F") a.element = Atom::get_element_number("F");
	else if (atom_type == "Cl") a.element = Atom::get_element_number("Cl");
	else if (atom_type == "Br") a.element = Atom::get_element_number("Br");
	else if (atom_type == "I") a.element = Atom::get_element_number("I");
	else if (atom_type == "Si") a.element = Atom::get_element_number("Si");
	else if (atom_type == "LP") a.element = 0;
	else if (atom_type == "Du") a.element = 0;
	else if (atom_type == "Na") a.element = Atom::get_element_number("Na");
	else if (atom_type == "K") a.element = Atom::get_element_number("K");
	else if (atom_type == "Ca") a.element = Atom::get_element_number("Ca");
	else if (atom_type == "Li") a.element = Atom::get_element_number("Li");
	else if (atom_type == "Al") a.element = Atom::get_element_number("Al");
	else if (atom_type == "Any") a.element = 0;
	else if (atom_type == "Hal") a.element = 0;
	else if (atom_type == "Het") a.element = 0;
	else if (atom_type == "Hev") a.element = 0;
	else a.element = 0;
}

//
int Mol2Entry::get_atom_serial_number(int atom1) const
{
	int n = 0;
	for (vector<AtomRecord>::const_iterator p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		if (p->id == atom1) return n;
		n++;
	}
	return -1;
}

//
int Mol2Entry::get_bond_type(const string& bond_type) const
{
	if (bond_type == "1") return 1;
	if (bond_type == "2") return 2;
	if (bond_type == "3") return 3;
	if (bond_type == "am") return 1;
	if (bond_type == "ar") return 4; // aromatic = 4 as defined in SDF
	if (bond_type == "du") return 1;
	if (bond_type == "un") return 1;
	if (bond_type == "nc") return 0;
	return 0;
}
