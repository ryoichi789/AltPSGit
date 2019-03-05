#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "define_secondary_structures.h"
#include "pdb_entry.h"
#include "pdb_reader.h"
#include "protein.h"
#include "secondary_structure.h"

using namespace Cba;

DefineSecondaryStructures the_app;

// DefineSecondaryStructures:
DefineSecondaryStructures::DefineSecondaryStructures() {}

// set_param_default:
void DefineSecondaryStructures::set_param_default()
{
	m_prog = "define_secondary_structures";
	m_pdbfile = 0;
	m_outfile = 0;
	m_seq = false;
}

// set_param:
void DefineSecondaryStructures::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!m_pdbfile) {
		usage();
		std::exit(2);
	}
}

// run:
int DefineSecondaryStructures::run()
{
	std::ifstream ifs(m_pdbfile);
	std::ostream* os;
	if (m_outfile) os = new std::ofstream(m_outfile);
	else os = &std::cout;

	PdbEntry ent;
	Protein* protein;
	std::vector<Protein*> proteins;
	if (PdbReader::read(ent,ifs)) {
		for (int i = 0; i < ent.nmolecules(); i++) {
			if (ent.type_of_molecule(i) == PdbEntry::PROTEIN &&
				ent.nmr_model(i) <= 1) {
				protein = new Protein;
				ent.create_protein(i, *protein);
				proteins.push_back(protein);
			}
		}
		if (!proteins.empty()) {
			SecondaryStructure s(proteins);
			if (m_seq) s.print_seq(*os);
			else s.print(*os);
		}
	}

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void DefineSecondaryStructures::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tf: PDB file [string] (unique; mandatory)" << endl;
	cerr << "\to: file for output [string] (unique; stdout by default)" << endl;
	cerr << "\ts: output in sequence format [t/f] (f by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void DefineSecondaryStructures::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "f") == 0) {
		m_obj->m_pdbfile = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_outfile = param.value;
	} else if (strcmp(param.key, "s") == 0) {
		if (param.value[0] == 't') m_obj->m_seq = true;
		else m_obj->m_seq = false;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}
