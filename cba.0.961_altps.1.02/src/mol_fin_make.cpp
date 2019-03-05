#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
//#include <vector>
//#include <string>
//#include <list>
//#include <map>
//#include <set>
#include <fstream>
//#include <iomanip>
//#include <cmath>
#include "mol_fin_make.h"
#include "sdf_reader.h"
#include "sdf_entry.h"
#include "molecule.h"
#include "mol_fin.h"

using namespace Cba;
//using std::vector;
//using std::string;
//using std::list;
//using std::map;
//using std::set;
//using std::istream;
//using std::ostream;
//using std::ifstream;
//using std::ofstream;

MolFinMake the_app;

//
MolFinMake::MolFinMake() {}

//
void MolFinMake::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\ti: SDF file [string]"
		<< " (unique; mandatory)" << endl;
	cerr << "\to: file for output [string]"
		<< " (unique; mandatory)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f]"
		<< " (f by default)" << endl;
}

//
void MolFinMake::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "i") == 0) {
		m_obj->m_fin = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_fout = param.value;
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
void MolFinMake::set_param_default()
{
	m_prog = "mol_fin_make";
	m_fin = 0;
	m_fout = 0;
}

//
void MolFinMake::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (m_fin == 0 || m_fout == 0) {
		usage();
		std::exit(2);
	}
}

//
int MolFinMake::run()
{
	std::ifstream ifs(m_fin);
	std::ofstream ofs(m_fout);

	SdfReader reader(ifs);
	SdfEntry ent;
	Molecule mol;
	MolFin mf;
	while (reader.read(ent)) {
		if (! ent.create_molecule(mol)) continue;
		mf.make(mol);
		mf.write(ofs);
	}

	return 0;
}
