#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "pdb_reader.h"
#include "pdb_entry.h"
#include "protein.h"
#include "bio_sequence.h"
#include "get_pdb_atom_aas.h"

using namespace Cba;
using std::string;
using std::vector;

GetPdbAtomAas the_app;

//
GetPdbAtomAas::GetPdbAtomAas() {}

//
void GetPdbAtomAas::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\ti: list of PDB file paths (plus relevant chain id if necessary) [string]"
		<< " (unique; mandatory)" << endl;
	cerr << "\to: file for output of FASTA sequences [string]"
		<< " (unique; mandatory)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f]"
		<< " (f by default)" << endl;
}

//
void GetPdbAtomAas::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "i") == 0) {
		m_obj->m_f_pdbs = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_f_out = param.value;
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
void GetPdbAtomAas::set_param_default()
{
	m_prog = "get_pdb_atom_aas";
	m_f_pdbs = 0;
	m_f_out = 0;
}

//
void GetPdbAtomAas::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (m_f_pdbs == 0 || m_f_out == 0) {
		usage();
		std::exit(2);
	}
}

//
int GetPdbAtomAas::run()
{
	get_file_paths();
	std::ofstream ofs(m_f_out);
	for (int i = 0; i < m_pdb_files.size(); i++) {
		get_atom_aas(m_pdb_files[i], m_pdb_chains[m_pdb_files[i]], ofs);
	}

	return 0;
}

//
void GetPdbAtomAas::get_file_paths()
{
	/*
	std::ifstream ifs(m_f_pdbs);
	string file_path;
	while (ifs >> file_path) {
		m_pdb_files.push_back(file_path);
	}
	*/
	std::ifstream ifs(m_f_pdbs);
	string line;
	size_t p, q;
	char c;
	string file_path;
	while (std::getline(ifs, line)) {
		if (line.length() == 0 || line.at(0) == ' ' || line.at(0) == '\t' || line.at(0) == '\n') continue;
		c = ' ';
		if ((p = line.find_first_of(" \t")) != string::npos) {
			file_path = line.substr(0, p);
			if ((q = line.find_first_not_of(" \t", p+1)) != string::npos) {
				c = line.at(q);
			}
		} else {
			file_path = line;
		}
		m_pdb_files.push_back(file_path);
		if (c != ' ' && c != '-' && c != '_') m_pdb_chains[file_path] += c;
	}
}

//
void GetPdbAtomAas::get_atom_aas(const string& file, const string& chains, std::ostream& os)
{
	std::ifstream ifs(file.c_str());
	if (! PdbReader::read(m_pdb_ent, ifs)) {
		std::cerr << "cannt read: " << file << std::endl;
		return;
	}

	string fname;
	get_fname(file, fname);

	if (chains.length() == 0) {
		for (int i = 0; i < m_pdb_ent.nmolecules(); i++) {
			if (m_pdb_ent.is_protein(i) && m_pdb_ent.nmr_model(i) <= 1) {
				get_atom_aas_for_mol_i(m_pdb_ent, i, fname, os);
			}
		}
	} else {
		for (int k = 0; k < chains.length(); k++) {
			for (int i = 0; i < m_pdb_ent.nmolecules(); i++) {
				if (m_pdb_ent.is_protein(i) && m_pdb_ent.nmr_model(i) <= 1 && m_pdb_ent.chainid(i) == chains.at(k)) {
					get_atom_aas_for_mol_i(m_pdb_ent, i, fname, os);
				}
			}
		}
	}
}

//
void GetPdbAtomAas::get_atom_aas_for_mol_i(const PdbEntry& ent, int i, const string& fname, std::ostream& os)
{
	string seqid;
	string seqannot = "";
	if (ent.create_protein(i, m_prot)) {
		get_seqid(fname, ent.chainid(i), seqid);
		m_prot.get_seq(m_seq);
		m_seq.set_id(seqid);
		if (m_debug) {
			std::ostringstream oss;
			oss << m_prot.residue(0).resseq << ' ' << m_prot.residue(m_prot.nresidues()-1).resseq;
			seqannot = oss.str();
		}
		m_seq.set_annot(seqannot);
		m_seq.print(os, 80);
	}
}

//
void GetPdbAtomAas::get_fname(const string& file, string& fname)
{
	size_t p;
	size_t p0 = 0;
	if ((p = file.find_last_of('/')) != string::npos) {
		p0 = p + 1;
	}
	size_t p1 = file.length();
	if ((p = file.find_last_of('.')) != string::npos) {
		p1 = p;
	}
	fname = file.substr(p0, p1-p0);
}

//
void GetPdbAtomAas::get_seqid(const string& file, char chainid, string& seqid)
{
	std::ostringstream oss;
	//if (chainid == ' ') chainid = '_';
	oss << file << '_' << chainid;
	seqid = oss.str();
}
