#ifndef CBA_GET_PDB_ATOM_AAS_H
#define CBA_GET_PDB_ATOM_AAS_H

// get amino acids of a protein chain as appeared in PDB ATOM records;

#include <functional>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "app.h"
#include "pdb_entry.h"
#include "protein.h"
#include "bio_sequence.h"

namespace Cba {

class GetPdbAtomAas : public App {
public:
	GetPdbAtomAas();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();

	const char* m_f_pdbs;
	const char* m_f_out;

	std::vector<std::string> m_pdb_files;
	std::map<std::string, std::string> m_pdb_chains; // pdb file -> relevant chains
	PdbEntry m_pdb_ent;
	Protein m_prot;
	BioSequence m_seq;

	void get_file_paths();
	void get_atom_aas(const std::string& file, const std::string& chains, std::ostream& os);
	void get_atom_aas_for_mol_i(const PdbEntry& ent, int i, const std::string& fname, std::ostream& os);
	void get_fname(const std::string& file, std::string& fname);
	void get_seqid(const std::string& file, char chainid, std::string& seqid);

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(GetPdbAtomAas* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		GetPdbAtomAas* m_obj;
	};
	friend class SetParam;
};

}
#endif
