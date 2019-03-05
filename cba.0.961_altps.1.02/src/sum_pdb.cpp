#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cstdio>
#include <sys/types.h> // not C++ standard
#include <sys/stat.h> // not C++ standard
#include <dirent.h> // not C++ standard
#include "sum_pdb.h"
#include "pdb_entry.h"
#include "pdb_reader.h"

using namespace Cba;

SumPdb the_app;

const char* SumPdb::PDB_FILE_EXT = ".pdb";

// SumPdb:
SumPdb::SumPdb() {}

// set_param_default:
void SumPdb::set_param_default()
{
	m_prog = "sum_pdb";
	m_output_file = 0;
}

// set_param:
void SumPdb::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (m_pdb_files.size() == 0 && m_pdb_dirs.size() == 0) {
		usage();
		std::exit(2);
	}
}

// run:
int SumPdb::run()
{
	if (m_output_file) m_os = new std::ofstream(m_output_file);
	else m_os = &std::cout;
	m_num = 0;

	std::vector<std::string>::const_iterator p;
	for (p = m_pdb_files.begin(); p != m_pdb_files.end(); ++p)
		sum_pdb(*p);
		
	for (p = m_pdb_dirs.begin(); p != m_pdb_dirs.end(); ++p)
		scan_dir(*p);

	if (m_os != &std::cout) delete m_os;
	return 0;
}

// usage:
void SumPdb::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tf: PDB file name [string] (multiple)" << endl;
	cerr << "\td: PDB root directory [string] (multiple)" << endl;
	cerr << "\to: file for output [string] (optional; stdout by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void SumPdb::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "f") == 0) {
		m_obj->m_pdb_files.push_back(param.value);
	} else if (strcmp(param.key, "d") == 0) {
		m_obj->m_pdb_dirs.push_back(param.value);
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_output_file = param.value;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

// get_pdb_filepaths:
void SumPdb::scan_dir(const std::string& dir)
{
	DIR* dp;
	struct dirent* ent;
	struct stat statbuf;
	std::string current_ent;

	if ((dp = opendir(dir.c_str())) == NULL) return;
	while ((ent = readdir(dp)) != NULL) {
		current_ent = dir + "/" + ent->d_name;
		stat(current_ent.c_str(), &statbuf);
		if (S_ISDIR(statbuf.st_mode)) {
			if (std::strcmp(".", ent->d_name) == 0 ||
				std::strcmp("..", ent->d_name) == 0) continue;
			scan_dir(current_ent.c_str());
		} else {
			int l = std::strlen(ent->d_name);
			if (l > LEN_PDB_FILE_EXT) {
				char* s = ent->d_name;
				s += l - LEN_PDB_FILE_EXT;
				if (std::strcmp(s, PDB_FILE_EXT) == 0)
					sum_pdb(current_ent);
			}
		}
	}
	closedir(dp);
}

// sum_pdb:
void SumPdb::sum_pdb(const std::string& filepath)
{
	std::ifstream ifs(filepath.c_str());
	PdbReader::read(m_ent, ifs);
	print_sum(filepath, m_ent);
}

// print_sum:
void SumPdb::print_sum(const std::string& filepath, const PdbEntry& ent)
{
	using std::endl;
	*m_os << ++m_num << " =======" << endl;
	*m_os << "<FILE> " << filepath << endl;
	*m_os << "<ID> " << ent.id() << endl;
	*m_os << "<TITLE> " << ent.title() << endl;
	*m_os << "<COMPOUND> " << ent.compound() << endl;
	*m_os << "<RESOLUTION> " << ent.resolution();
	if (ent.nmolecules() && ent.nmr_model(0)) *m_os << " (model)";
	*m_os << endl;
	*m_os << "<NO_OF_MOLECULES> " << ent.nmolecules() << endl;
	for (int i = 0; i < ent.nmolecules(); i++)
		print_mol_info(ent, i);
	*m_os << endl;
}

// print_mol_info:
void SumPdb::print_mol_info(const PdbEntry& ent, int imol)
{
	*m_os << imol+1 << ":\t";
	*m_os << "c=" << ent.chainid(imol) << '\t';
	*m_os << "r=" << ent.nresidues(imol) << '\t';
	*m_os << "a=" << ent.natoms(imol) << '\t';

	PdbEntry::MoleculeType mt = ent.type_of_molecule(imol);
	if (mt == PdbEntry::UNDEFINED)
		*m_os << "undefined";
	else if (mt == PdbEntry::PROTEIN)
		*m_os << "protein";
	else if (mt == PdbEntry::PEPTIDE)
		*m_os << "peptide";
	else if (mt == PdbEntry::DNA)
		*m_os << "dna";
	else if (mt == PdbEntry::RNA)
		*m_os << "rna";
	else if (mt == PdbEntry::SMALL)
		*m_os << "small";
	else if (mt == PdbEntry::WATER)
		*m_os << "water";
	else if (mt == PdbEntry::ATOM)
		*m_os << "atom";

	*m_os << '\t';
	if (mt == PdbEntry::SMALL || mt == PdbEntry::WATER || mt == PdbEntry::ATOM)
		*m_os << ent.resname(imol, 0);

	*m_os << std::endl;
}
