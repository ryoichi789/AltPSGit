#include <vector>
#include <string>
#include "sdf_entry.h"
#include "exception.h"

using namespace Cba;
using std::vector;
using std::string;

//
SdfEntry::SdfEntry() {}

// name_of_annot:
const string& SdfEntry::name_of_annot(size_t i) const
{
	if (i >= nannots()) throw Exception();
	return m_annots[i].name;
}

// annot:
const vector<string>& SdfEntry::annot(size_t i) const
{
	if (i >= nannots()) throw Exception();
	return m_annots[i].annot;
}

// create_molecule:
bool SdfEntry::create_molecule(Molecule& mol) const
{
	mol.clear();
	mol.set_id(id());
	mol.set_name(name());
	mol.set_type(Molecule::SMALL);

	Atom a;
	for (vector<AtomRecord>::const_iterator
		p = m_atoms.begin(); p != m_atoms.end(); ++p) {
		a.name = p->symbol;
		a.element = Atom::get_element_number(a.name.c_str());
		a.pos.move_to(p->x, p->y, p->z);
		mol.add_atom(a);
	}

	for (vector<BondRecord>::const_iterator
		p = m_bonds.begin(); p != m_bonds.end(); ++p) {
		mol.add_bond(mol.atom(p->atom1), mol.atom(p->atom2), p->bond_type);
	}

	return true;
}

// clear:
void SdfEntry::clear()
{
	m_header.clear();
	m_id.erase();
	m_name.erase();
	m_atoms.clear();
	m_bonds.clear();
	m_annots.clear();
}

// Header::clear:
void SdfEntry::Header::clear()
{
	molecule_name.erase();
	user_initials.erase();
	program_name.erase();
	date.erase();
	dimension = 0;
	comments.erase();
}
