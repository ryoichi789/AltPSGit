#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include "align_seq.h"
#include "bio_sequence.h"
#include "fasta_reader.h"
#include "bio_seq_score_matrix.h"
#include "bio_seq_cmp_dp.h"

using namespace Cba;

AlignSeq the_app;

// AlignSeq:
AlignSeq::AlignSeq() {}

// set_param_default:
void AlignSeq::set_param_default()
{
	m_prog = "align_seq";
	m_dna = false;
	m_fasta1 = 0;
	m_fasta2 = 0;
	m_output_file = 0;
	m_score_file = 0;
	m_max_nalignments = MAXNALIGN;
	m_min_score = MINSCORE;
}

// set_param:
void AlignSeq::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!(m_fasta1 && m_fasta2)) {
		usage();
		std::exit(2);
	}
}

// run:
int AlignSeq::run()
{
	std::ostream* os;
	if (m_output_file) os = new std::ofstream(m_output_file);
	else os = &std::cout;

	BioSequence seq1;
	BioSequence seq2;
	BioSeqCmpDp dp(m_max_nalignments, m_min_score, m_score_file, m_dna);

	std::ifstream ifs1(m_fasta1);
	FastaReader fas1(ifs1);
	while (fas1.read(seq1)) {
		std::ifstream ifs2(m_fasta2);
		FastaReader fas2(ifs2);
		while (fas2.read(seq2))
			dp.compare(&seq1, &seq2, *os);
	}

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void AlignSeq::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\t1: sequence file in FASTA format [string] (mandatory)" << endl;
	cerr << "\t2: sequence file in FASTA format [string] (mandatory)" << endl;
	cerr << "\to: file for output [string] (optional; stdout by default)" << endl;
	cerr << "\tn: number of alignments to find from a pair of sequences [int] (optional; "
		<< MAXNALIGN << " by default)" << endl;
	cerr << "\ts: acceptable alignment score [int] (optional; "
		<< MINSCORE << " by default)" << endl;
	cerr << "\tf: file for amino acid substitution matrix [string] (optional; BLOSUM62 by default)" << endl;
	cerr << "\tdna: to compare DNA sequences [t/f] (optional; false by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void AlignSeq::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "1") == 0) {
		m_obj->m_fasta1 = param.value;
	} else if (strcmp(param.key, "2") == 0) {
		m_obj->m_fasta2 = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_output_file = param.value;
	} else if (strcmp(param.key, "n") == 0) {
		m_obj->m_max_nalignments = std::atoi(param.value);
	} else if (strcmp(param.key, "s") == 0) {
		m_obj->m_min_score = std::atoi(param.value);
	} else if (strcmp(param.key, "f") == 0) {
		m_obj->m_score_file = param.value;
	} else if (strcmp(param.key, "dna") == 0) {
		if (param.value[0] == 't') m_obj->m_dna = true;
		else m_obj->m_dna = false;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}
