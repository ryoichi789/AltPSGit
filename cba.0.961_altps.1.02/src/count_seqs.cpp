#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "count_seqs.h"
#include "fasta_reader.h"
#include "bio_sequence.h"

using namespace Cba;

CountSeqs the_app;

// CountSeqs:
CountSeqs::CountSeqs() {}

// set_param_default:
void CountSeqs::set_param_default()
{
	m_progname.assign(m_prog = "count_seqs");
	m_seqfile = 0;
}

// set_param:
void CountSeqs::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	m_progname.assign(m_prog);
	int i = m_progname.find_last_of('/'); // to remove directory path
	if (i < m_progname.length()) m_progname = m_progname.substr(i+1);
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!m_seqfile) {
		usage();
		std::exit(2);
	}
}

// run:
int CountSeqs::run()
{
	std::ifstream ifs(m_seqfile);
	FastaReader reader(ifs);
	BioSequence seq;
	int nseqs = 0;
	int naas = 0;

	std::vector<std::string>::iterator vci;
	int i;
	while (reader.read(seq)) {
		nseqs++;
		naas += seq.length();
		if (m_progname.compare("get_seqs") == 0) {
			for (vci = m_ids.begin(); vci != m_ids.end(); ++vci)
				if ((i = seq.id().find(*vci)) < seq.id().length() &&
					(i == 0 || seq.id().at(i-1) == '|')) break;
			if (vci != m_ids.end()) {
				seq.print(std::cout);
				m_ids.erase(vci);
				if (m_ids.empty()) break;
			}
		}
	}
	if (m_progname.compare("count_seqs") == 0)
		std::cout << nseqs << '\t' << naas << std::endl;

	return 0;
}

// usage:
void CountSeqs::usage()
{
	using std::cerr;
	using std::endl;
	using std::strcmp;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tf: sequence file in FASTA format [string] (unique)" << endl;
	if (m_progname.compare("count_seqs") == 0) {
		cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
	} else {
		cerr << "\tc: ID code of a sequence entry in file [string] (multiple)" << endl;
	}
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void CountSeqs::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "f") == 0) {
		m_obj->m_seqfile = param.value;
	} else if (strcmp(param.key, "c") == 0) {
		m_obj->m_ids.push_back(param.value);
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}
