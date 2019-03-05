#ifndef CBA_SHOW_ATOM_TYPES_H
#define CBA_SHOW_ATOM_TYPES_H

#include <functional>
#include <vector>
#include <iostream>
#include "app.h"
#include "molecule.h"

namespace Cba {

class ShowAtomTypes : public App {
public:
	ShowAtomTypes();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	std::vector<const char*> m_files;
	const char* m_outfile;
	bool m_pdb;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(ShowAtomTypes* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		ShowAtomTypes* m_obj;
	};
	friend class SetParam;

	void show_atom_types(const char* file, std::ostream& os) const;
	void print(const Molecule& mol, std::ostream& os) const;
};

}
#endif
