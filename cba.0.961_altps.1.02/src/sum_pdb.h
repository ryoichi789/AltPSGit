#ifndef CBA_SUM_PDB_H
#define CBA_SUM_PDB_H

#include <functional>
#include <vector>
#include <string>
#include <iostream>
#include "app.h"
#include "pdb_entry.h"

namespace Cba {

class SumPdb : public App {
public:
	SumPdb();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

	static const char* PDB_FILE_EXT; // ".pdb"
	static const int LEN_PDB_FILE_EXT = 4; // length of ".pdb"

private:
	void usage();
	std::vector<std::string> m_pdb_files;
	std::vector<std::string> m_pdb_dirs;
	const char* m_output_file;
	std::ostream* m_os; // for output
	PdbEntry m_ent;
	int m_num; // serial number of PDB files processed

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(SumPdb* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		SumPdb* m_obj;
	};
	friend class SetParam;

	void scan_dir(const std::string& dir);
	void sum_pdb(const std::string& filepath);
	void print_sum(const std::string& filepath, const PdbEntry& ent);
	void print_mol_info(const PdbEntry& ent, int imol);
};

}
#endif
