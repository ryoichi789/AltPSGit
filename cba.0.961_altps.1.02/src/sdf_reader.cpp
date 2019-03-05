#include <string>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "sdf_reader.h"
#include "sdf_entry.h"

using namespace Cba;

//
SdfReader::SdfReader(std::istream& is) : m_is(is), m_serial_number(0) {}

// read:
bool SdfReader::read(SdfEntry& ent)
{
	ent.clear();
	if (read_header(ent) && read_counts_line(ent) &&
		read_atoms(ent) && read_bonds(ent) && skip_properties_block()) {
		read_annots(ent);
		set_id(ent);
		set_name(ent);
		return true;
	}
	return false;
}

// read_header:
bool SdfReader::read_header(SdfEntry& ent)
{
	using std::getline;

	// line 1
	if (! getline(m_is, ent.m_header.molecule_name)) return false;

	// line 2
	std::string line;
	if (! getline(m_is, line)) return false;
	ent.m_header.user_initials = line.substr(0, 2);
	ent.m_header.program_name = line.substr(2, 8);
	ent.m_header.date = line.substr(10, 10);
	if (line.substr(20, 2).compare("2D") == 0)
		ent.m_header.dimension = 2;
	else if (line.substr(20, 2).compare("3D") == 0)
		ent.m_header.dimension = 3;
	else
		ent.m_header.dimension = 0;

	// line 3
	if (! getline(m_is, ent.m_header.comments)) return false;

	return true;
}

// read_counts_line:
bool SdfReader::read_counts_line(SdfEntry& ent)
{
	std::string line;
	if (! getline(m_is, line)) return false;
	ent.m_natoms = std::atoi(line.substr(0, 3).c_str());
	ent.m_nbonds = std::atoi(line.substr(3, 3).c_str());
	return true;
}

// read_atoms:
bool SdfReader::read_atoms(SdfEntry& ent)
{
	SdfEntry::AtomRecord a;
	std::string line;

	for (int i = 0; i < ent.m_natoms; i++) {
		if (! std::getline(m_is, line)) return false;
		a.x = std::atof(line.substr(0, 10).c_str());
		a.y = std::atof(line.substr(10, 10).c_str());
		a.z = std::atof(line.substr(20, 10).c_str());
		a.symbol = line.substr(31, 3);
		a.symbol.erase(a.symbol.find(' ')); // remove trailing blanks
		a.charge = 0;
		int c = std::atoi(line.substr(36, 3).c_str());
		if (c == 1) a.charge = 3;
		else if (c == 2) a.charge = 2;
		else if (c == 3) a.charge = 1;
		else if (c == 5) a.charge = -1;
		else if (c == 6) a.charge = -2;
		else if (c == 7) a.charge = -3;
		ent.m_atoms.push_back(a);
	}

	return true;
}

// read_bonds:
bool SdfReader::read_bonds(SdfEntry& ent)
{
	SdfEntry::BondRecord b;
	std::string line;

	for (int i = 0; i < ent.m_nbonds; i++) {
		if (! getline(m_is, line)) return false;
		b.atom1 = std::atoi(line.substr(0, 3).c_str()) - 1;
		b.atom2 = std::atoi(line.substr(3, 3).c_str()) - 1;
		b.bond_type = std::atoi(line.substr(6, 3).c_str());
		ent.m_bonds.push_back(b);
	}

	return true;
}

// skip_properties_block:
bool SdfReader::skip_properties_block()
{
	std::string line;
	while (std::getline(m_is, line))
		if (line.substr(0, 6).compare("M  END") == 0) break;
	if (line.substr(0, 6).compare("M  END")) return false;
	return true;
}

// read_annots:
void SdfReader::read_annots(SdfEntry& ent)
{
	SdfEntry::Annot annot;
	std::string line;
	int p1, p2;
	while (std::getline(m_is, line)) {
		if (line.at(0) == '>') { // beginning of a data item
			p1 = line.find_first_of('<', 1);
			if (p1++ == std::string::npos) continue;
			p2 = line.find_first_of('>', p1);
			if (p2-- == std::string::npos) continue;
			annot.name = line.substr(p1, p2-p1+1);

			annot.registry_number = 0;
			p1 = line.find_first_of('(', p2+2);
			if (p1++ != std::string::npos) {
				p2 = line.find_first_of(')', p1);
				if (p2-- != std::string::npos) {
					int i;
					for (i = p1; i <= p2; i++)
						if (! std::isdigit(line.at(i))) break;
					if (i > p2) annot.registry_number =
							std::atoi(line.substr(p1,p2-p1+1).c_str());
				}
			}

			annot.annot.clear();
			while (std::getline(m_is, line)) {
				if (line.length() == 0) break;
				if (line.length() == 1 && std::isspace(line.at(0))) break;
				annot.annot.push_back(line);
			}

			ent.m_annots.push_back(annot);
		} else if (line.substr(0, 4).compare("$$$$") == 0) {
			break; // end of an entry
		}
	}
}

// set_id:
void SdfReader::set_id(SdfEntry& ent)
{
	if (ent.m_header.molecule_name.length()) {
		ent.m_id = ent.m_header.molecule_name;
		std::replace(ent.m_id.begin(), ent.m_id.end(), ' ', '_');
	} else {
		std::ostringstream oss;
		oss << ++m_serial_number;
		ent.m_id = oss.str();
	}
}

// set_name:
void SdfReader::set_name(SdfEntry& ent)
{
	if (ent.m_header.molecule_name.length())
		ent.m_name = ent.m_header.molecule_name;
	else
		ent.m_name = "NA";
}
