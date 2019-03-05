#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <vector>
//#include <string>
//#include <list>
//#include <map>
//#include <set>
#include <fstream>
//#include <iomanip>
//#include <cmath>
#include "morgan_numbering.h"
#include "sdf_reader.h"
#include "sdf_entry.h"
#include "molecule.h"

using namespace Cba;
using std::vector;
//using std::string;
//using std::list;
//using std::map;
//using std::set;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::endl;

MorganNumbering the_app;

// NOTE: define static const member variables here if any

//
MorganNumbering::MorganNumbering() {}

//
void MorganNumbering::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\ti: SDF file [string]"
		<< " (unique; mandatory)" << endl;
	cerr << "\to: file for output [string]"
		<< " (unique; stdout by default)" << endl;
	cerr << "\th: count hydrogens [t/f]"
		<< " (unique; f by default)" << endl;
	cerr << "\tv: verbose [t/f]"
		<< " (unique; f by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f]"
		<< " (f by default)" << endl;
}

//
void MorganNumbering::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "i") == 0) {
		m_obj->m_fin = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_fout = param.value;
	} else if (strcmp(param.key, "h") == 0) {
		if (param.value[0] == 't') m_obj->m_count_hydrogens = true;
		else m_obj->m_count_hydrogens = false;
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
void MorganNumbering::set_param_default()
{
	m_prog = "morgan_numbering";
	m_fin = 0;
	m_fout = 0;
	m_count_hydrogens = false;
	m_verbose = false;
}

//
void MorganNumbering::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (m_fin == 0) {
		usage();
		std::exit(2);
	}
}

//
int MorganNumbering::run()
{
	ostream* os;
	if (m_fout) os = new ofstream(m_fout);
	else os = &std::cout;

	ifstream ifs(m_fin);
	SdfReader reader(ifs);
	SdfEntry ent;
	Molecule mol;
	vector<int> morgan_numbers;
	vector<vector<int> > morgan_hist;
	while (reader.read(ent)) {
		if (ent.create_molecule(mol)) {
			if (! m_count_hydrogens) mol.remove_hydrogens();
			mol.get_morgan_numbering(morgan_numbers, morgan_hist);
			*os << mol.id() << " natoms=" << mol.natoms();
			for (int i = 0; i < morgan_numbers.size(); i++) {
				*os << ' ' << morgan_numbers[i]+1;
			}
			*os << endl;
			if (m_verbose) {
				for (int k = 0; k < morgan_hist.size(); k++) {
					*os << '\t' << "step_" << k
						<< ' ' << "c=" << morgan_hist[k].back();
					for (int i = 0; i < morgan_hist[k].size()-1; i++) {
						*os << ' ' << morgan_hist[k][i];
					}
					*os << endl;
				}
			}
		} else {
			*os << "ERROR" << endl;
		}
	}

	if (os != &std::cout) delete os;
	return 0;
}
