#include <iostream>
#include <string>
#include <cstring>
#include <cctype>
#include "cba_entry.h"
#include "atom.h"
#include "molecule.h"
#include "protein.h"

using namespace Cba;
using std::vector;
using std::string;
using std::istream;

//
CbaEntry::CbaEntry() {}

//
CbaEntry::~CbaEntry() {}

//
void CbaEntry::clear()
{
	m_id.erase();
	m_name.erase();
	m_type.erase();
	m_atoms.clear();
	m_bonds.clear();
	m_residues.clear();
	m_comments.clear();
}

//
bool CbaEntry::read(istream& is)
{
	clear();

	string line;

	// I record
	while (std::getline(is, line))
		if (line.at(0) == 'I') break; // skip to the first record
	if (line.empty() || line.at(0) != 'I') return false;
	read_i(line);

	// other records
	while (std::getline(is, line)) {
		if (line.empty() || line.at(0) == '/') break; // end of entry
		switch (line.at(0)) {
		case 'N': read_n(line); break;
		case 'T': read_t(line); break;
		case 'A': read_a(line); break;
		case 'B': read_b(line); break;
		case 'R': read_r(line); break;
		case 'C': read_c(line); break;
		}
	}

	return true;
}

//
void CbaEntry::read_i(const string& line)
{
	if (line.length() >= 3 && line.at(1) == '\t')
		m_id = line.substr(2);
}

//
void CbaEntry::read_n(const string& line)
{
	if (line.length() >= 3 && line.at(1) == '\t')
		m_name = line.substr(2);
}

//
void CbaEntry::read_t(const string& line)
{
	if (line.length() >= 3 && line.at(1) == '\t')
		m_type = line.substr(2);
}

//
void CbaEntry::read_a(const string& line)
{
	if (line.length() < 3 && line.at(1) != '\t') return;
	Atom a;
	double x, y, z;

	int p0 = 2; // starting position of a field in the line
	int p1 = p0; // ending position of a field
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	a.number = std::atoi(line.substr(p0, p1-p0).c_str());
	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	a.name = line.substr(p0, p1-p0);
	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	a.element = Atom::get_element_number(line.substr(p0, p1-p0).c_str());
	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	x = std::atof(line.substr(p0, p1-p0).c_str());
	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	y = std::atof(line.substr(p0, p1-p0).c_str());
	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	z = std::atof(line.substr(p0, p1-p0).c_str());
	a.pos.move_to(x, y, z);
	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	a.asa = std::atof(line.substr(p0, p1-p0).c_str());
	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	//a.vdwr = std::atof(line.substr(p0, p1-p0).c_str());
	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	a.pc_class = std::atoi(line.substr(p0, p1-p0).c_str());
	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	//a.type = std::atoi(line.substr(p0, p1-p0).c_str());
	p0 = ++p1;
	if (line.at(p0) != '0' && line.at(p0) != '1') return;
	a.properties = line.substr(p0, Atom::number_of_properties).c_str();

	m_atoms.push_back(a);
}

//
void CbaEntry::read_b(const string& line)
{
	if (line.length() < 3 && line.at(1) != '\t') return;
	BRecord b;

	int p0 = 2; // starting position of a field in the line
	int p1 = p0; // ending position of a field
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	b.a1 = std::atoi(line.substr(p0, p1-p0).c_str());

	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	b.a2 = std::atoi(line.substr(p0, p1-p0).c_str());

	p0 = ++p1;
	if (! std::isdigit(line.at(p0))) return;
	b.v = std::atoi(line.substr(p0, 1).c_str());

	m_bonds.push_back(b);
}

//
void CbaEntry::read_r(const string& line)
{
	if (line.length() < 3 && line.at(1) != '\t') return;
	RRecord r;

	int p0 = 2; // starting position of a field in the line
	int p1 = p0; // ending position of a field
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	r.serial_num = std::atoi(line.substr(p0, p1-p0).c_str());

	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	r.serial_num_orig = line.substr(p0, p1-p0);

	p0 = ++p1;
	if (((p1 = line.find('\t', p0)) == string::npos) || (p0 == p1)) return;
	r.name = line.substr(p0, p1-p0);

	p0 = ++p1;
	r.end_atom = std::atoi(line.substr(p0, line.length()-p0).c_str());

	m_residues.push_back(r);
}

//
bool CbaEntry::create_molecule(Molecule& m) const
{
	m.clear();
	m.set_id(m_id);
	m.set_name(m_name);
	if (m_type == "protein") m.set_type(Molecule::PROTEIN);
	else if (m_type == "peptide") m.set_type(Molecule::PEPTIDE);
	else if (m_type == "dna") m.set_type(Molecule::DNA);
	else if (m_type == "rna") m.set_type(Molecule::RNA);
	else if (m_type == "small") m.set_type(Molecule::SMALL);
	else if (m_type == "water") m.set_type(Molecule::WATER);
	else if (m_type == "atom") m.set_type(Molecule::ATOM);
	else m.set_type(Molecule::UNDEFINED);

	for (vector<Atom>::const_iterator
		vci = m_atoms.begin(); vci != m_atoms.end(); ++vci)
		m.add_atom(*vci);

	for (vector<BRecord>::const_iterator
		vci = m_bonds.begin(); vci != m_bonds.end(); ++vci) {
		const BRecord& brec = *vci;
		m.add_bond(m.atom(brec.a1), m.atom(brec.a2), brec.v);
	}

	return true;
}

//
bool CbaEntry::create_protein(Protein& p) const
{
	p.clear();
	create_molecule(p);

	// define residues
	Protein::Residue residue;
	Atom* atom;
	int a = 0; // serial number of atom
	for (int i = 0; i < m_residues.size(); i++) {
		residue.clear();
		for (int j = a; j <= m_residues[i].end_atom; j++) {
			atom = p.atom(j);
			if (atom->name.compare(" N  ") == 0) residue.n = atom;
			else if (atom->name.compare(" H  ") == 0) residue.h = atom;
			else if (atom->name.compare(" CA ") == 0) residue.ca = atom;
			else if (atom->name.compare(" C  ") == 0) residue.c = atom;
			else if (atom->name.compare(" O  ") == 0) residue.o = atom;
			else if (atom->name.compare(" OXT") == 0) residue.oxt = atom;
			else if (atom->name.compare(" CB ") == 0) residue.cb = atom;
			else residue.side_chain_atoms.push_back(atom);
		}
		residue.name = m_residues[i].name;
		residue.number = m_residues[i].serial_num_orig;
		p.add_residue(residue);
	}

	return true;
}

//
void CbaEntry::read_c(const string& line)
{
	if (line.length() < 3 && line.at(1) != '\t') return;
	m_comments.push_back(line.substr(2));
}
