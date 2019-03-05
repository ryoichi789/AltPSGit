#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include "show_pdb_molinfo.h"
#include "pdb_entry.h"
#include "pdb_reader.h"
#include "molecule.h"
#include "protein.h"

using namespace Cba;

ShowPdbMolinfo the_app;

// ShowPdbMolinfo:
ShowPdbMolinfo::ShowPdbMolinfo() {}

// set_param_default:
void ShowPdbMolinfo::set_param_default()
{
	m_prog = "show_pdb_molinfo";
	m_outfile = 0;
}

// set_param:
void ShowPdbMolinfo::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (m_pdbfiles.size() == 0) {
		usage();
		std::exit(2);
	}
}

// run:
int ShowPdbMolinfo::run()
{
	using std::cout;
	using std::endl;

	std::ostream* os;
	if (m_outfile) os = new std::ofstream(m_outfile);
	else os = &std::cout;

	PdbEntry ent;
	Molecule mol;
	Protein prot;
	for (std::vector<const char*>::const_iterator
		p = m_pdbfiles.begin(); p != m_pdbfiles.end(); ++p) {
		std::ifstream ifs(*p);
		PdbReader::read(ent, ifs);
		*os << "------------------------------" << endl;
		*os << ent.nmolecules() << " molecules in " << *p << endl;
		*os << "------------------------------" << endl;
		for (int i = 0; i < ent.nmolecules(); i++) {
			*os << "molecule " << i+1 << ":" << endl;
			if (ent.type_of_molecule(i) == PdbEntry::PROTEIN) {
				ent.create_protein(i, prot);
				show_protein_info(prot, *os);
			} else {
				ent.create_molecule(i, mol);
				show_molecule_info(mol, *os);
			}
		}
	}

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void ShowPdbMolinfo::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tf: PDB file name [string] (multiple)" << endl;
	cerr << "\to: file for output [string] (unique; stdout by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void ShowPdbMolinfo::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "f") == 0) {
		m_obj->m_pdbfiles.push_back(param.value);
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_outfile = param.value;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

// show_molecule_info:
void ShowPdbMolinfo::show_molecule_info(const Molecule& mol, std::ostream& os)
{
	using std::endl;
	show_molecule_type(mol, os);
	os << "\tid: " << mol.id() << endl;
	os << "\tname: " << mol.name() << endl;
	os << "\tnatoms: " << mol.natoms() << endl;
}

// show_protein_info:
void ShowPdbMolinfo::show_protein_info(const Protein& prot, std::ostream& os)
{
	using std::endl;
	show_molecule_info(prot, os);
	os << "\tnresidues: " << prot.nresidues() << endl;

	// for each residue,
	// show name, its number, number of atoms,
	// name of its representative atom, and its position
	for (int i = 0; i < prot.nresidues(); i++) {
		const Protein::Residue& r = prot.residue(i);
		os << "\t\t" << i+1 << ": "
			<< r.name << " " << r.number << " " << r.natoms() << " ";
		const Atom* a = r.ca;
		if (!a) a = r.c;
		if (!a) a = r.n;
		os << a->name << " "
			<< a->pos.x() << " " << a->pos.y() << " " << a->pos.z() << endl;
	}
}

// show_molecule_type:
void ShowPdbMolinfo::show_molecule_type(const Molecule& mol, std::ostream& os)
{
	using std::endl;

	os << "\ttype: ";
	if (mol.type() == Molecule::PROTEIN) os << "protein" << endl;
	else if (mol.type() == Molecule::PEPTIDE) os << "peptide" << endl;
	else if (mol.type() == Molecule::DNA) os << "dna" << endl;
	else if (mol.type() == Molecule::RNA) os << "rna" << endl;
	else if (mol.type() == Molecule::SMALL) os << "small molecule" << endl;
	else if (mol.type() == Molecule::WATER) os << "water" << endl;
	else if (mol.type() == Molecule::ATOM) os << "atom" << endl;
	else os << "unknown" << endl;
}
