#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cstring>
#include "mol2_reader.h"
#include "mol2_entry.h"

using namespace Cba;
using std::string;

const char* Mol2Reader::record_names[] = {
	"@<TRIPOS>ALT_TYPE",
	"@<TRIPOS>ANCHOR_ATOM",
	"@<TRIPOS>ASSOCIATED_ANNOTATION",
	"@<TRIPOS>ATOM",
	"@<TRIPOS>BOND",
	"@<TRIPOS>CENTER_OF_MASS",
	"@<TRIPOS>CENTROID",
	"@<TRIPOS>COMMENT",
	"@<TRIPOS>CRYSIN",
	"@<TRIPOS>DICT",
	"@<TRIPOS>DATA_FILE",
	"@<TRIPOS>EXTENSION_POINT",
	"@<TRIPOS>FF_PBC",
	"@<TRIPOS>FFCON_ANGLE",
	"@<TRIPOS>FFCON_DIST",
	"@<TRIPOS>FFCON_MULTI",
	"@<TRIPOS>FFCON_RANGE",
	"@<TRIPOS>FFCON_TORSION",
	"@<TRIPOS>LINE",
	"@<TRIPOS>LSPLANE",
	"@<TRIPOS>MOLECULE",
	"@<TRIPOS>NORMAL",
	"@<TRIPOS>QSAR_ALIGN_RULE",
	"@<TRIPOS>RING_CLOSURE",
	"@<TRIPOS>ROTATABLE_BOND",
	"@<TRIPOS>SEARCH_DIST",
	"@<TRIPOS>SEARCH_OPTIONS",
	"@<TRIPOS>SET",
	"@<TRIPOS>SUBSTRUCTURE",
	"@<TRIPOS>U_FEAT",
	"@<TRIPOS>UNITY_ATOM_ATTR",
	"@<TRIPOS>UNITY_BOND_ATTR",
	0
};

//
bool Mol2Reader::read(Mol2Entry& ent, std::istream& from)
{
	ent.clear();
	int record_type = get_first_record_type(from);
	while (record_type < NO_OF_RECORD_TYPES) {
		if (record_type == MOLECULE) {
			record_type = read_molecule_records(ent, from);
		} else if (record_type == ATOM) {
			record_type = read_atom_records(ent, from);
		} else if (record_type == BOND) {
			record_type = read_bond_records(ent, from);
		} else {
			record_type = skip_records(from);
		}
	}
	if (ent.natoms()) return true;
	return false;
}

//
int Mol2Reader::get_first_record_type(std::istream& from)
{
	string line;
	while (std::getline(from, line)) {
		if (line.length() && line.at(0) == '@') return get_record_type(line);
	}
	return NO_OF_RECORD_TYPES;
}

//
int Mol2Reader::get_record_type(const string& line)
{
	int i;
	for (i = 0; record_names[i]; i++) {
		if (line.length() != std::strlen(record_names[i])) continue;
		if (line == record_names[i]) break;
	}
	return i;
}

//
int Mol2Reader::read_molecule_records(Mol2Entry& ent, std::istream& from)
{
	string line;

	if (! std::getline(from, line)) return NO_OF_RECORD_TYPES;
	ent.m_mol.name = line;

	if (! std::getline(from, line)) return NO_OF_RECORD_TYPES;
	std::istringstream iss(line);
	iss >> ent.m_mol.natoms >> ent.m_mol.nbonds;

	if (! std::getline(from, line)) return NO_OF_RECORD_TYPES;
	ent.m_mol.type = line;

	if (! std::getline(from, line)) return NO_OF_RECORD_TYPES;
	ent.m_mol.charge_type = line;

	if (! std::getline(from, line)) return NO_OF_RECORD_TYPES;
	ent.m_mol.status_bit = line;

	if (! std::getline(from, line)) return NO_OF_RECORD_TYPES;
	ent.m_mol.comment = line;

	return skip_records(from);
}

//
int Mol2Reader::read_atom_records(Mol2Entry& ent, std::istream& from)
{
	string line;
	Mol2Entry::AtomRecord atom_record;
	for (int i = 0; i < ent.m_mol.natoms; ) {
		if (! std::getline(from, line)) return NO_OF_RECORD_TYPES;
		if (line.length() && line.at(0) == '#') continue;
		read_atom_record(atom_record, line);
		ent.m_atoms.push_back(atom_record);
		i++;
	}
	return skip_records(from);
}

//
void Mol2Reader::read_atom_record(Mol2Entry::AtomRecord& atom_record, const string& from)
{
	std::istringstream iss(from);
	iss >> atom_record.id >> atom_record.name
		>> atom_record.x >> atom_record.y >> atom_record.z
		>> atom_record.type
		>> atom_record.subst_id >> atom_record.subst_name
		>> atom_record.charge;
}

//
int Mol2Reader::read_bond_records(Mol2Entry& ent, std::istream& from)
{
	string line;
	Mol2Entry::BondRecord bond_record;
	for (int i = 0; i < ent.m_mol.nbonds; ) {
		if (! std::getline(from, line)) return NO_OF_RECORD_TYPES;
		if (line.length() && line.at(0) == '#') continue;
		read_bond_record(bond_record, line);
		ent.m_bonds.push_back(bond_record);
		i++;
	}
	return skip_records(from);
}

//
void Mol2Reader::read_bond_record(Mol2Entry::BondRecord& bond_record, const string& from)
{
	std::istringstream iss(from);
	iss >> bond_record.id
		>> bond_record.atom1 >> bond_record.atom2
		>> bond_record.type;
}

//
int Mol2Reader::skip_records(std::istream& from)
{
	string line;
	while (std::getline(from, line)) {
		if (line.length() && line.at(0) == '@') return get_record_type(line);
	}
	return NO_OF_RECORD_TYPES;
}
