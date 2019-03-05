#ifndef CBA_SHOW_PDB_MOLINFO_H
#define CBA_SHOW_PDB_MOLINFO_H

#include <functional>
#include <vector>
#include <iostream>
#include "app.h"
#include "molecule.h"
#include "protein.h"

namespace Cba {

class ShowPdbMolinfo : public App {
public:
	ShowPdbMolinfo();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	std::vector<const char*> m_pdbfiles;
	const char* m_outfile;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(ShowPdbMolinfo* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		ShowPdbMolinfo* m_obj;
	};
	friend class SetParam;

	void show_molecule_info(const Molecule& mol, std::ostream& os);
	void show_protein_info(const Protein& mol, std::ostream& os);
	void show_molecule_type(const Molecule& mol, std::ostream& os);
};

}
#endif
