#ifndef CBA_CALC_ASAS_H
#define CBA_CALC_ASAS_H

#include <functional>
#include "app.h"
#include "molecule.h"

namespace Cba {

class CalcAsas : public App {
public:
	CalcAsas();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	const char* m_molfile;
	const char* m_outfile;
	bool m_sdf;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(CalcAsas* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		CalcAsas* m_obj;
	};
	friend class SetParam;

	void print_asas(const Molecule& mol, std::ostream& os) const;
};

}
#endif
