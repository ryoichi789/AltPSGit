#include <iostream>
#include <ios>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>
#include "predict_tm.h"
#include "trans_membrane_prediction.h"
#include "bio_sequence.h"
#include "fasta_reader.h"

using namespace Cba;

PredictTm the_app;

// PredictTm:
PredictTm::PredictTm() {}

// set_param_default:
void PredictTm::set_param_default()
{
	m_prog = "predict_tm";
	m_fasta_file = 0;
	m_out_file = 0;
}

// set_param:
void PredictTm::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!m_fasta_file) {
		usage();
		std::exit(2);
	}
}

// run:
int PredictTm::run()
{
	std::ifstream ifs(m_fasta_file);
	std::ostream* os;
	if (m_out_file) os = new std::ofstream(m_out_file);
	else os = &std::cout;

	FastaReader fas(ifs);
	BioSequence seq;
	TransMembranePrediction tm;
	std::vector<TransMembranePrediction::Region> rs;
	while (fas.read(seq)) {
		tm.predict(seq, rs);
		tm.sort(rs);
		print(seq, rs, *os);
	}

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void PredictTm::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\ti: sequence file in Fasta format [string] (mandatory)" << endl;
	cerr << "\to: file for output [string] (optional; stdout by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void PredictTm::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "i") == 0) {
		m_obj->m_fasta_file = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_out_file = param.value;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

// print:
void PredictTm::print(const BioSequence& seq, 
	const std::vector<TransMembranePrediction::Region>& rs,   
	std::ostream& os) const
{
	if (rs.size() == 0) return;
	os << '>' << seq.id() << ' ' << rs.size() << std::endl;
	for (std::vector<TransMembranePrediction::Region>::const_iterator
		i = rs.begin(); i != rs.end(); ++i) {
		os << i->s0+1 << '\t' << i->e0+1 << '\t'
			<< i->s+1 << '\t' << i->e+1 << '\t'
			<< std::fixed << std::setprecision(2) << i->h << std::endl;
	}
}
