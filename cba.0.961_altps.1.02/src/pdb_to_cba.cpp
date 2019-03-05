#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string>
#include "pdb_to_cba.h"
#include "pdb_reader.h"
#include "pdb_entry.h"
#include "protein.h"
#include "molecule.h"

using namespace Cba;

PdbToCba the_app;

//
PdbToCba::PdbToCba() {}

//
void PdbToCba::set_param_default()
{
	m_prog = "pdb_to_cba";
	m_pdb_file = 0;
	m_cba_file = 0;
	m_asa = false;
	m_protein_only = false;
}

//
void PdbToCba::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));

	if (!m_pdb_file) {
		usage();
		exit(2);
	}

	if (m_cba_file && std::strcmp(m_pdb_file, m_cba_file) == 0) {
		std::cerr << "output file must be different from input file" << std::endl;
		usage();
		exit(2);
	}
}

//
int PdbToCba::run()
{
	std::ifstream ifs(m_pdb_file);
	PdbReader reader;
	reader.read(m_pdb, ifs);
	Protein prot;
	Molecule mol;

	if (!m_cba_file) make_cba_file_name();
	std::ofstream ofs(m_cba_file);
	for (int i = 0; i < m_pdb.nmolecules(); i++) {
		if (! is_target(i)) continue;
		if (m_pdb.type_of_molecule(i) == PdbEntry::PROTEIN) {
			m_pdb.create_protein(i, prot);
			prot.assign_atom_types();
			if (m_asa) prot.calc_atom_asas();
			prot.print(ofs);
			prot.print_entry_terminator(ofs);
		} else if (m_pdb.type_of_molecule(i) == PdbEntry::SMALL ||
			m_pdb.type_of_molecule(i) == PdbEntry::PEPTIDE) {
			m_pdb.create_molecule(i, mol);
			mol.assign_atom_types();
			if (m_asa) mol.calc_atom_asas();
			mol.print(ofs);
			mol.print_entry_terminator(ofs);
		}
	}
}

//
void PdbToCba::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\ti: input file in PDB format [string] (unique; mandatory)" << endl;
	cerr << "\to: output file in CBA format [string] (unique; optional; 'pdb_code.cba' by default)" << endl;
	cerr << "\tc: target chain in PDB file [char] (optional; all chains by default)" << endl;
	cerr << "\tasa: calculate solvent-accessible surface areas of atoms "
		"if true [t/f] (optional; false by default)" << endl;
	cerr << "\tprot: process only proteins if true [t/f] (optional; false by default)" << endl;
}

//
void PdbToCba::SetParam::operator() (const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "i") == 0) {
		if (m_obj->m_pdb_file == 0) m_obj->m_pdb_file = param.value;
	} else if (strcmp(param.key, "o") == 0) {
		if (m_obj->m_cba_file == 0) m_obj->m_cba_file = param.value;
	} else if (strcmp(param.key, "c") == 0) {
		m_obj->m_chains.push_back(param.value[0]);
	} else if (strcmp(param.key, "asa") == 0) {
		if (param.value[0] == 't') m_obj->m_asa = true;
		else if (param.value[0] == 'f') m_obj->m_asa = false;
	} else if (strcmp(param.key, "prot") == 0) {
		if (param.value[0] == 't') m_obj->m_protein_only = true;
		else if (param.value[0] == 'f') m_obj->m_protein_only = false;
	}
}

//
bool PdbToCba::is_target(int imol) const
{
	if (m_protein_only && m_pdb.type_of_molecule(imol) != PdbEntry::PROTEIN) return false;
	if (m_chains.empty()) return true;
	if (std::find(m_chains.begin(), m_chains.end(), m_pdb.chainid(imol)) == m_chains.end()) return false;
	return true;
}

//
void PdbToCba::make_cba_file_name()
{
	if (m_cba_file) return;

	std::string fn;
	fn = m_pdb.id() + ".cba";
	m_cba_file = new char[fn.length() + 1];
	std::strcpy(m_cba_file, fn.c_str());
	for (int i = 0; i < fn.length(); i++)
		m_cba_file[i] = std::tolower(m_cba_file[i]);
}
