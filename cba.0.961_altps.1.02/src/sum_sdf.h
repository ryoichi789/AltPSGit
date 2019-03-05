#ifndef CBA_SUM_SDF_H
#define CBA_SUM_SDF_H

#include <functional>
#include <iostream>
#include "app.h"

namespace Cba {

class Molecule;
class SdfEntry;

class SumSdf : public App {
public:
	SumSdf();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	const char* m_sdffile;
	const char* m_outfile;
	bool m_hydrogen; // whether to include hydrogens
	bool m_long; // whether to show annotations

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(SumSdf* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		SumSdf* m_obj;
	};
	friend class SetParam;

	void print_mol_info(const Molecule& mol, std::ostream& os) const;
	void print_sdf_annots(const SdfEntry& ent, std::ostream& os) const;
};

}
#endif
