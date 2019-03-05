#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "find_identical_mols.h"
#include "sdf_reader.h"
#include "sdf_entry.h"
#include "molecule.h"

using namespace Cba;
using std::vector;
using std::string;
using std::list;
using std::map;
using std::set;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;

FindIdenticalMols the_app;

//
FindIdenticalMols::FindIdenticalMols() {}

//
void FindIdenticalMols::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tq: SDF file for a query molecule [string]"
		<< " (unique; mandatory)" << endl;
	cerr << "\td: SDF file for molecule database to be searched [string]"
		<< " (unique; mandatory)" << endl;
	cerr << "\to: file for output [string]"
		<< " (unique; stdout by default)" << endl;
	cerr << "\th: remove hydrogens [t/f]"
		<< " (unique; t by default)" << endl;
	cerr << "\tv: output in a verbose form [t/f]"
		<< " (unique; f by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f]"
		<< " (f by default)" << endl;
}

//
void FindIdenticalMols::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "q") == 0) {
		m_obj->m_sdf_query = param.value;
	} else if (strcmp(param.key, "d") == 0) {
		m_obj->m_sdf_db = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_ofile = param.value;
	} else if (strcmp(param.key, "h") == 0) {
		if (param.value[0] == 't') m_obj->m_remove_hydrogens = true;
		else m_obj->m_remove_hydrogens = false;
	} else if (strcmp(param.key, "v") == 0) {
		if (param.value[0] == 't') m_obj->m_verbose = true;
		else m_obj->m_verbose = false;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter "
			<< param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

//
void FindIdenticalMols::set_param_default()
{
	m_prog = "find_identical_mols";
	m_sdf_query = 0;
	m_sdf_db = 0;
	m_ofile = 0;
	m_remove_hydrogens = true;
	m_verbose = false;
}

//
void FindIdenticalMols::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (m_sdf_query == 0 || m_sdf_db == 0) {
		usage();
		std::exit(2);
	}
}

//
int FindIdenticalMols::run()
{
	Molecule molq;
	get_query_molecule(molq);
	if (molq.natoms() == 0) {
		std::cerr << "query molecule empty" << std::endl;
		return 1;
	}
	if (m_remove_hydrogens) molq.remove_hydrogens();

	ostream* osp;
	if (m_ofile) {
		osp = new ofstream(m_ofile);
	} else {
		osp = &std::cout;
	}

	ifstream ifs(m_sdf_db);
	SdfReader reader(ifs);
	SdfEntry ent;
	Molecule moldb;
	vector<vector<int> > alignments;
	int n = 0;
	while (reader.read(ent)) {
		n++;
		if (ent.create_molecule(moldb)) {
			if (m_remove_hydrogens) moldb.remove_hydrogens();
			if (molq.is_identical_to(moldb, alignments)) {
				print_hit(*osp, n, moldb, alignments);
			}
		}
		if (m_debug && n%1000==0) {
			std::cerr << n << " molecules scanned" << std::endl;
		}
	}

	if (osp != & std::cout) delete osp;

	return 0;
}

//
void FindIdenticalMols::get_query_molecule(Molecule& mol) const
{
	mol.clear();
	ifstream ifs(m_sdf_query);
	SdfReader reader(ifs);
	SdfEntry ent;
	if (reader.read(ent)) {
		ent.create_molecule(mol);
	}
}

//
void FindIdenticalMols::print_hit(ostream& os, int n, const Molecule& mol, vector<vector<int> > alignments) const
{
	if (m_verbose) {
		os << '>' << n << '\t' << mol.id() << std::endl;
		os << alignments.size() << std::endl;
		for (int i = 0; i < alignments.size(); i++) {
			for (int j = 0; j < alignments[i].size(); j++) {
				os << ' ' << alignments[i][j]+1;
			}
			os << std::endl;
		}
	} else {
		os << n << '\t' << mol.id() << std::endl;
	}
}
