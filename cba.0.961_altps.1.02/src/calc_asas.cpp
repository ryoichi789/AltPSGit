#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <ios>
#include <iomanip>
#include "calc_asas.h"
#include "sdf_entry.h"
#include "sdf_reader.h"
#include "pdb_entry.h"
#include "pdb_reader.h"
#include "molecule.h"
#include "protein.h"

using namespace Cba;

CalcAsas the_app;

// CalcAsas:
CalcAsas::CalcAsas() {}

// set_param_default:
void CalcAsas::set_param_default()
{
	m_prog = "calc_asas";
	m_molfile = 0;
	m_sdf = false;
}

// set_param:
void CalcAsas::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!m_molfile) {
		usage();
		std::exit(2);
	}
}

// run:
int CalcAsas::run()
{
	std::ifstream ifs(m_molfile);
	std::ostream* os;
	if (m_outfile) os = new std::ofstream(m_outfile);
	else os = &std::cout;

	Molecule mol;
	if (m_sdf) {
		SdfEntry ent;
		SdfReader reader(ifs);
		Molecule mol;
		while (reader.read(ent)) {
			if (ent.create_molecule(mol)) {
				mol.assign_atom_types();
				mol.calc_atom_asas();
				print_asas(mol, *os);
			}
		}
	} else {
		PdbEntry ent;
		PdbReader::read(ent, ifs);
		Protein prot;
		for (int i = 0; i < ent.nmolecules(); i++) {
			if (ent.type_of_molecule(i) == PdbEntry::PROTEIN &&
				ent.create_protein(i, prot)) {
				prot.assign_atom_types();
				prot.calc_atom_asas();
				print_asas(prot, *os);
			}
		}
	}

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void CalcAsas::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tf: file of molecule (in SDF or PDB) [string] (unique; mandatory)" << endl;
	cerr << "\to: file for output [string] (unique; stdout by default)" << endl;
	cerr << "\tsdf: true if in SDF [t/f] (f by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void CalcAsas::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "f") == 0) {
		m_obj->m_molfile = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_outfile = param.value;
	} else if (strcmp(param.key, "sdf") == 0) {
		if (param.value[0] == 't') m_obj->m_sdf = true;
		else m_obj->m_sdf = false;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

// print_asas:
void CalcAsas::print_asas(const Molecule& mol, std::ostream& os) const
{
	os << mol.id() << ' ' << mol.natoms() << std::endl;
	const Atom* a;
	for (int i = 0; i < mol.natoms(); i++) {
		a = mol.atom(i);
		os << i+1 << '\t' << a->name << '\t'
			<< std::fixed << std::setw(8) << std::setprecision(3)
			<< a->asa << std::endl;
	}
	os << std::endl;
}
