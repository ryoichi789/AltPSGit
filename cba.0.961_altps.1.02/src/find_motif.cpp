#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <ios>
#include <iomanip>
#include "find_motif.h"
#include "fasta_reader.h"
#include "bio_sequence.h"
#include "region.h"

using namespace Cba;

FindMotif the_app;

// FindMotif:
FindMotif::FindMotif() {}

// set_param_default:
void FindMotif::set_param_default()
{
	m_prog = "find_motif";
	m_motif = 0;
	m_fasta_file = 0;
	m_out_file = 0;
}

// set_param:
void FindMotif::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!(m_motif && m_fasta_file)) {
		usage();
		std::exit(2);
	}
}

// run:
int FindMotif::run()
{
	std::ifstream ifs(m_fasta_file);
	std::ostream* os;
	if (m_out_file) os = new std::ofstream(m_out_file);
	else os = &std::cout;

	FastaReader fas(ifs);
	BioSequence seq;
	Regions regions;
	while (fas.read(seq)) {
		seq.find_motif(m_motif, regions);
		for (int i = 0; i < regions.size(); i++)
			print_match(seq, regions[i], *os);
	}

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void FindMotif::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tq: sequence motif in regular expression [string] (mandatory)" << endl;
	cerr << "\tt: sequence file in Fasta format [string] (mandatory)" << endl;
	cerr << "\to: file for output [string] (optional; stdout by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void FindMotif::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "q") == 0) {
		m_obj->m_motif = param.value;
	} else if (strcmp(param.key, "t") == 0) {
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

// print_match:
void FindMotif::print_match(
	const BioSequence& seq, const Region& region, std::ostream& os) const
{
	os << std::setw(40) << std::left << seq.id() << ' ';

	// N-terminal- (or 5'-) flanking region
	for (int i = region.first-L_FLANKING_REGION; i < region.first; i++)
		if (i < 0) os << '_';
		else os << seq.at(i);

	// matching region
	os << ' ' << std::setw(5) << std::right << region.first+1 << ' ';
	for (int i = region.first; i <= region.second; i++) os << seq.at(i);
	os << ' ' << std::setw(5) << std::left << region.second+1 << ' ';

	// C-terminal- (or 3'-) flanking region
	for (int i = region.second+1;
		i <= (region.second+L_FLANKING_REGION) && i < seq.length(); i++)
		os << seq.at(i);

	os << std::endl;
}
