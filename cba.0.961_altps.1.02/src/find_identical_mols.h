#ifndef CBA_FIND_IDENTICAL_MOLS_H
#define CBA_FIND_IDENTICAL_MOLS_H

#include <functional>
#include <vector>
#include <iostream>
#include "app.h"
#include "molecule.h"

namespace Cba {

class FindIdenticalMols : public App {
public:
	FindIdenticalMols();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();

	const char* m_sdf_query;
	const char* m_sdf_db;
	const char* m_ofile;
	bool m_verbose;
	bool m_remove_hydrogens;

	void get_query_molecule(Molecule& mol) const;
	void print_hit(std::ostream& os, int n, const Molecule& mol, std::vector<std::vector<int> > alignments) const;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(FindIdenticalMols* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		FindIdenticalMols* m_obj;
	};
	friend class SetParam;
};

}
#endif
