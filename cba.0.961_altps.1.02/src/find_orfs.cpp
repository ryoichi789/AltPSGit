#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "find_orfs.h"
#include "fasta_reader.h"
#include "dna_sequence.h"

using namespace Cba;

FindOrfs the_app;

// FindOrfs:
FindOrfs::FindOrfs() {}

// set_param_default:
void FindOrfs::set_param_default()
{
	m_prog = "find_orfs";
	m_seq_file = 0;
	m_out_file = 0;
	m_print_fasta = false;
	m_min_len = MIN_LEN;
}

// set_param:
void FindOrfs::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!m_seq_file) {
		usage();
		std::exit(2);
	}
}

// run:
int FindOrfs::run()
{
	std::ifstream ifs(m_seq_file);
	std::ostream* os;
	if (m_out_file) os = new std::ofstream(m_out_file);
	else os = &std::cout;

	FastaReader reader(ifs);
	DnaSequence dna;
	while (reader.read(dna))
		print_orfs(dna, *os);

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void FindOrfs::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\ti: sequence file in FASTA format [string] (mandatory)" << endl;
	cerr << "\to: file for output [string] (optional; stdout by default)" << endl;
	cerr << "\tfasta: output ORF sequence in FASTA format [t/f] (optional; f by default)" << endl;
	cerr << "\tl: minimum length of ORF [positive int] ("
		<< MIN_LEN << " by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void FindOrfs::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "i") == 0) {
		m_obj->m_seq_file = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_out_file = param.value;
	} else if (strcmp(param.key, "l") == 0) {
		m_obj->m_min_len = std::atoi(param.value);
	} else if (strcmp(param.key, "fasta") == 0) {
		if (param.value[0] == 't') m_obj->m_print_fasta = true;
		else m_obj->m_print_fasta = false;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

// print_orfs:
void FindOrfs::print_orfs(DnaSequence& dna, std::ostream& os) const
{
	int norfs = dna.find_orfs();
	if (m_print_fasta) {
		BioSequence seq;
		for (int i = 0; i < norfs && dna.translate_orf(i, seq, m_min_len); i++)
			seq.print(os);
	} else {
		os << '>' << dna.id() << ' ' << dna.length() << std::endl;
		dna.print_orfs(os, m_min_len);
	}
}
