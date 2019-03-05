#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "show_atom_types.h"
#include "pdb_reader.h"
#include "pdb_entry.h"
#include "sdf_reader.h"
#include "sdf_entry.h"
#include "molecule.h"
#include "protein.h"
#include "atom.h"

using namespace Cba;

ShowAtomTypes the_app;

// ShowAtomTypes:
ShowAtomTypes::ShowAtomTypes() {}

// set_param_default:
void ShowAtomTypes::set_param_default()
{
	m_prog = "show_atom_types";
	m_outfile = 0;
	m_pdb = false;
}

// set_param:
void ShowAtomTypes::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (m_files.empty()) {
		usage();
		std::exit(2);
	}
}

// run:
int ShowAtomTypes::run()
{
	std::ostream* os;
	if (m_outfile) os = new std::ofstream(m_outfile);
	else os = &std::cout;

	for (std::vector<const char*>::const_iterator
		p = m_files.begin(); p != m_files.end(); ++p)
		show_atom_types(*p, *os);

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void ShowAtomTypes::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tf: SDF or PDB file [string] (multiple)" << endl;
	cerr << "\to: file for output [string] (unique; stdout by default)" << endl;
	cerr << "\tpdb: set true if file is in PDB format [t/f] (unique; f by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void ShowAtomTypes::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "f") == 0) {
		m_obj->m_files.push_back(param.value);
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_outfile = param.value;
	} else if (strcmp(param.key, "pdb") == 0) {
		if (param.value[0] == 't') m_obj->m_pdb = true;
		else m_obj->m_pdb = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

// show_atom_types:
void ShowAtomTypes::show_atom_types(const char* file, std::ostream& os) const
{
	std::ifstream ifs(file);
	if (m_pdb) {
		PdbReader reader;
		PdbEntry ent;
		Protein prot;
		reader.read(ent, ifs);
		for (int i = 0; i < ent.nmolecules(); i++) {
			if (ent.type_of_molecule(i) == PdbEntry::PROTEIN &&
				ent.create_protein(i, prot)) {
				prot.assign_atom_types();
				print(prot, os);
			}
		}
	} else {
		SdfReader reader(ifs);
		SdfEntry ent;
		Molecule mol;
		while (reader.read(ent)) {
			if (ent.create_molecule(mol)) {
				mol.assign_atom_types();
				print(mol, os);
			}
		}
	}
}

// print:
void ShowAtomTypes::print(const Molecule& mol, std::ostream& os) const
{
	using std::endl;
	os << "ID: " << mol.id() << endl;
	os << "NAME: " << mol.name() << endl;
	os << "NATOMS: " << mol.natoms() << endl;
	int c;
	for (int i = 0; i < mol.natoms(); i++) {
		os << i+1 << '\t' << mol.atom(i)->name << '\t';
		c = mol.atom(i)->pc_class;
		if (c == Atom::CATION) os << "CATION     \t";
		else if (c == Atom::ANION) os << "ANION      \t";
		else if (c == Atom::DONOR) os << "DONOR      \t";
		else if (c == Atom::ACCEPTOR) os << "ACCEPTOR   \t";
		else if (c == Atom::POLAR) os << "POLAR      \t";
		else if (c == Atom::HYDROPHOBIC) os << "HYDROPHOBIC\t";
		else if (c == Atom::NONE) os << "NONE       \t";
		else os << "UNDEFINED  \t";
		os << mol.atom(i)->properties.str();
		os << endl;
	}
	os << endl;
}
