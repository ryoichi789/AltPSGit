#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
//#include <vector>
#include <string>
#include <list>
//#include <map>
//#include <set>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cctype>
#include "mol_fin_cmp.h"
#include "mol_fin.h"

using namespace Cba;
//using std::vector;
using std::string;
using std::list;
//using std::map;
//using std::set;
//using std::istream;
//using std::ostream;
//using std::ifstream;
//using std::ofstream;

MolFinCmp the_app;

const int MolFinCmp::HITS_CAPACITY = 3;

//
MolFinCmp::MolFinCmp() {}

//
void MolFinCmp::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tq: query mol_fin file [string]"
		<< " (unique; mandatory)" << endl;
	cerr << "\tt: target mol_fin file [string]"
		<< " (multiple; optional;" << endl
		<< "\t\tbut at least one target file must be specified" << endl
		<< "\t\teither by this or 'l' option)" << endl;
	cerr << "\tl: file containing list of target mol_fin files [string]"
		<< " (optional)" << endl;
	cerr << "\tc: cutoff for Tanimoto coefficient [double]"
		<< " (unique; 0.0 by default)" << endl;
	cerr << "\tp: cutoff for p-value [double]"
		<< " (unique; 1.0 by default;" << endl
		<< "\t\tignored if 'c' option is specified)" << endl;
	cerr << "\ttop: output top-n hits for each query [int]" << endl
		<< "\t\t(unique; 0 by default; output all when 0;" << endl
		<< "\t\tignored if 'c' or 'p' option is specified)" << endl;
	cerr << "\to: file for output [string]"
		<< " (unique; stdout by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f]"
		<< " (f by default)" << endl;
}

//
void MolFinCmp::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "q") == 0) {
		m_obj->m_fquery = param.value;
	} else if (strcmp(param.key, "t") == 0) {
		m_obj->m_ftargets.push_back(string(param.value));
	} else if (strcmp(param.key, "l") == 0) {
		m_obj->m_ftargetlist = param.value;
	} else if (strcmp(param.key, "c") == 0) {
		m_obj->m_tc_cut = std::atof(param.value);
	} else if (strcmp(param.key, "p") == 0) {
		m_obj->m_pv_cut = std::atof(param.value);
	} else if (strcmp(param.key, "top") == 0) {
		m_obj->m_top_n = std::atoi(param.value);
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
void MolFinCmp::set_param_default()
{
	m_prog = "mol_fin_cmp";
	m_fquery = 0;
	m_ftargetlist = 0;
	m_tc_cut = 0.0;
	m_pv_cut = 1.0;
	m_top_n = 0;
	m_fout = 0;
}

//
void MolFinCmp::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (m_fquery == 0) {
		usage();
		std::exit(2);
	}
}

//
int MolFinCmp::run()
{
	set_targets();
	m_hits_capacity = m_top_n * HITS_CAPACITY;

	std::ostream* os;
	if (m_fout) os = new std::ofstream(m_fout);
	else os = &std::cout;

	print_search_info(*os);

	MolFin mfq;
	string idq;
	std::ifstream ifsq(m_fquery);
	for (int nq = 0; mfq.read(ifsq); nq++) {
		get_molid(m_fquery, nq, mfq, idq);
		scan_targets(nq, idq, mfq, *os);
	}

	if (os != &std::cout) delete os;
	return 0;
}

//
void MolFinCmp::print_search_info(std::ostream& os)
{
	using std::endl;
	os << "## Tc>= " << m_tc_cut << endl;
	os << "## P<=  " << m_pv_cut << endl;
	os << "## top  " << m_top_n << endl;
	os << "## query_file: " << m_fquery << endl;
	os << "## no_of_target_files: " << m_ftargets.size() << endl;
	for (int i = 0; i < m_ftargets.size(); i++) {
		os << "##  " << i+1 << " = " << m_ftargets[i] << endl;
	}
	os << endl;
}

//
void MolFinCmp::scan_targets(int nq, const string& idq, MolFin& mfq, std::ostream& os)
{
	m_hits.clear();

	MolFin mft;
	string idt;
	Hit hit;
	int nt;
	int ntall = 0;
	for (int i = 0; i < m_ftargets.size(); i++) {
		std::ifstream ifst(m_ftargets[i].c_str());
		hit.fileno = i+1;
		if (m_debug) { std::cerr << "query:" << idq << "\ttarget:" << m_ftargets[i] << "\tscanning..."; }
		for (nt = 0; mft.read(ifst); nt++) {
			get_molid(m_ftargets[i].c_str(), nt, mft, hit.idt);
			if ((hit.tc = mfq.tc(mft)) < m_tc_cut) continue;
			if ((hit.pv = MolFin::pvalue(hit.tc)) > m_pv_cut) continue;
			if (m_top_n) {
				if (m_hits.size() == m_hits_capacity) {
					m_hits.sort(CmpHit());
					m_hits.resize(m_top_n);
				}
				m_hits.push_back(hit);
			} else {
				os << idq << '\t' << hit.fileno << '\t' << hit.idt << '\t'
					<< std::fixed << std::setprecision(3) << std::setw(5) << hit.tc << '\t'
					<< std::fixed << std::setprecision(4) << std::setw(6) << hit.pv << std::endl;
			}
		}
		ntall += nt;
		if (m_debug) { std::cerr << "done." << std::endl; }
	}
	if (m_top_n) {
		m_hits.sort(CmpHit());
		print_hits(os, idq, nq+1, ntall);
	}
}

//
void MolFinCmp::print_hits(std::ostream& os, const string& idq, int nq, int nt)
{
	using std::endl;
	os << "#query_no: " << nq << endl;
	os << "#query_id: " << idq << endl;
	os << "#ntargets: " << nt << endl;
	os << "#top " << m_top_n << " hits..." << endl;
	list<Hit>::const_iterator lci = m_hits.begin();
	for (int i = 0; i < m_top_n && lci != m_hits.end(); i++, ++lci) {
		os << i+1 << '\t' << lci->fileno << '\t' << lci->idt << '\t'
			<< std::fixed << std::setprecision(3) << std::setw(5) << lci->tc << '\t'
			<< std::fixed << std::setprecision(4) << std::setw(6) << lci->pv << endl;
	}
	os << endl;
}

//
void MolFinCmp::set_targets()
{
	if (m_ftargetlist) {
		std::ifstream ifs(m_ftargetlist);
		string s;
		while (ifs >> s) {
			m_ftargets.push_back(s);
		}
	}

	if (m_ftargets.empty()) {
		usage();
		std::exit(2);
	}
}

//
void MolFinCmp::get_molid(const char* file, int i, const MolFin& mf, std::string& id) const
{
	const string& id0 = mf.id();
	if (id0.length() && (! std::isspace(id0.at(0)))) {
		id = id0;
		return;
	}

	bool no_id = false;
	if (id0.length() == 0) {
		no_id = true;
	} else {
		int k;
		for (k = 0; k < id0.length(); k++) {
			if (! std::isspace(id0.at(k))) break;
		}
		if (k == id0.length()) no_id = true;
	}
	if (no_id) {
		std::ostringstream oss;
		oss << file << '_' << i+1;
		id = oss.str();
		return;
	}
	id = id0;
}
