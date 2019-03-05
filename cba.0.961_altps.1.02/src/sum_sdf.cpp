#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "sum_sdf.h"
#include "sdf_entry.h"
#include "sdf_reader.h"
#include "molecule.h"

using namespace Cba;

SumSdf the_app;

// SumSdf:
SumSdf::SumSdf() {}

// set_param_default:
void SumSdf::set_param_default()
{
	m_prog = "sum_sdf";
	m_sdffile = 0;
	m_outfile = 0;
	m_hydrogen = true;
	m_long = false;
}

// set_param:
void SumSdf::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!m_sdffile) {
		usage();
		std::exit(2);
	}
}

// run:
int SumSdf::run()
{
	std::ifstream ifs(m_sdffile);
	std::ostream* os;
	if (m_outfile) os = new std::ofstream(m_outfile);
	else os = &std::cout;

	SdfReader reader(ifs);
	SdfEntry ent;
	Molecule mol;
	int i = 0;
	while (reader.read(ent)) {
		ent.create_molecule(mol);
		if (!m_hydrogen) mol.remove_hydrogens();
		*os << "--------------------" <<std::endl;
		*os << "Molecule: " << ++i << std::endl;
		*os << "--------------------" <<std::endl;
		print_mol_info(mol, *os);
		if (m_long) print_sdf_annots(ent, *os);
	}

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void SumSdf::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tf: SDF file name [string] (unique; mandatory)" << endl;
	cerr << "\to: file for output [string] (unique; stdout by default)" << endl;
	cerr << "\th: include hydrogens [t/f] (true by default)" << endl;
	cerr << "\tl: show annotations [t/f] (false by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void SumSdf::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "f") == 0) {
		m_obj->m_sdffile = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_outfile = param.value;
	} else if (strcmp(param.key, "h") == 0) {
		if (param.value[0] == 't') m_obj->m_hydrogen = true;
		else m_obj->m_hydrogen = false;
	} else if (strcmp(param.key, "l") == 0) {
		if (param.value[0] == 't') m_obj->m_long = true;
		else m_obj->m_long = false;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

// print_mol_info:
void SumSdf::print_mol_info(const Molecule& mol, std::ostream& os) const
{
	using std::endl;
	os << "id: " << mol.id() << endl;
	os << "name: " << mol.name() << endl;
	if (mol.type() == Molecule::SMALL) os << "type: small" << endl;
	else os << "type: large" << endl;
	os << "natoms: " << mol.natoms() << endl;
	os << "nbonds: " << mol.nbonds() << endl;
}

// print_sdf_annots:
void SumSdf::print_sdf_annots(const SdfEntry& ent, std::ostream& os) const
{
	using std::endl;
	os << "nannots: " << ent.nannots() << endl;
	for (int i = 0; i < ent.nannots(); i++) {
		os << '\t' << ent.name_of_annot(i) << endl;
		const std::vector<std::string>& annot = ent.annot(i);
		for (std::vector<std::string>::const_iterator
			p = annot.begin(); p != annot.end(); ++p) {
			os << "\t\t" << *p << endl;
		}
	}
}
