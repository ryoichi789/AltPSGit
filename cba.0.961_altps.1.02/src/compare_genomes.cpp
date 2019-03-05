#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "compare_genomes.h"
#include "fasta_reader.h"
#include "bio_sequence.h"
#include "bio_seq_cmp_identical.h"
using namespace Cba;

CompareGenomes the_app;

// CompareGenomes:
CompareGenomes::CompareGenomes() {}

// set_param_default:
void CompareGenomes::set_param_default()
{
	m_prog = "compare_genomes";
	m_dna_file1 = 0;
	m_dna_file2 = 0;
	m_output_file = 0;
	m_min_len_equiv_subseq = MIN_LEN_EQUIV_SUBSEQ;
}

// set_param:
void CompareGenomes::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!(m_dna_file1 && m_dna_file2)) {
		usage();
		std::exit(2);
	}
}

// run:
int CompareGenomes::run()
{
	std::ostream* os;
	if (m_output_file) os = new std::ofstream(m_output_file);
	else os = &std::cout;

	BioSequence seq1;
	BioSequence seq2;
	BioSeqCmpIdentical cmp(m_min_len_equiv_subseq);

	std::ifstream ifs1(m_dna_file1);
	FastaReader fas1(ifs1);
	while (fas1.read(seq1)) {
		std::ifstream ifs2(m_dna_file2);
		FastaReader fas2(ifs2);
		while (fas2.read(seq2)) {
			cmp.compare(&seq1, &seq2, *os);
			*os << std::endl;
		}
	}

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void CompareGenomes::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\t1: first genome DNA sequence file [string] (mandatory)" << endl;
	cerr << "\t2: second genome DNA sequence file [string] (mandatory)" << endl;
	cerr << "\to: output file [string] (optional; stdout by default)" << endl;
	cerr << "\tm: minimum length of equivalent DNA subsequences [int] (unique; ";
	cerr << MIN_LEN_EQUIV_SUBSEQ << " by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void CompareGenomes::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "1") == 0) {
		m_obj->m_dna_file1 = param.value;
	} else if (strcmp(param.key, "2") == 0) {
		m_obj->m_dna_file2 = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_output_file = param.value;
	} else if (strcmp(param.key, "m") == 0) {
		m_obj->m_min_len_equiv_subseq = std::atoi(param.value);
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}
