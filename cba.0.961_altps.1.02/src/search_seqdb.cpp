#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include "search_seqdb.h"
#include "bio_seq_db_search.h"

using namespace Cba;

SearchSeqdb the_app;

// SearchSeqdb:
SearchSeqdb::SearchSeqdb() {}

// set_param_default:
void SearchSeqdb::set_param_default()
{
	m_prog = "search_seqdb";
	m_query_fasta_file = 0;
	m_target_fasta_file = 0;
	m_output_file = 0;
	m_score_min = SCORE_MIN;
	m_length_min = LENGTH_MIN;
	m_make_tables = false;
}

// set_param:
void SearchSeqdb::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!m_target_fasta_file) {
		usage();
		std::exit(2);
	}
}

// run:
int SearchSeqdb::run()
{
	BioSeqDbSearch dbsearch;
	if (m_make_tables) {
		dbsearch.create_tables(m_target_fasta_file);
	} else if (m_query_fasta_file) {
		dbsearch.search(m_query_fasta_file, m_target_fasta_file,
			m_score_min, m_length_min, m_output_file);
	}

	return 0;
}

// usage:
void SearchSeqdb::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tq: query FASTA file [string] (unique)" << endl;
	cerr << "\tt: target FASTA file [string] (unique; mandatory)" << endl;
	cerr << "\to: file for output [string] (unique; stdout by default)" << endl;
	cerr << "\ts: minimum acceptable score for matching region [int] (unique; "
		<< SCORE_MIN << " by default)" << endl;
	cerr << "\tl: minimum acceptable length of matching region [int] (unique; "
		<< LENGTH_MIN << " by default)" << endl;
	cerr << "\ttbl: to make segment tables from target FASTA file [t/f] (optional; f by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void SearchSeqdb::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "q") == 0) {
		m_obj->m_query_fasta_file = param.value;
	} else if (strcmp(param.key, "t") == 0) {
		m_obj->m_target_fasta_file = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_output_file = param.value;
	} else if (strcmp(param.key, "s") == 0) {
		m_obj->m_score_min = std::atoi(param.value);
	} else if (strcmp(param.key, "l") == 0) {
		m_obj->m_length_min = std::atoi(param.value);
	} else if (strcmp(param.key, "tbl") == 0) {
		if (param.value[0] == 't') m_obj->m_make_tables = true;
		else m_obj->m_make_tables = false;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}
