#ifndef CBA_DEFINE_SECONDARY_STRUCTURES_H
#define CBA_DEFINE_SECONDARY_STRUCTURES_H

#include <functional>
#include "app.h"

namespace Cba {

class DefineSecondaryStructures : public App {
public:
	DefineSecondaryStructures();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	const char* m_pdbfile;
	const char* m_outfile;
	bool m_seq;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(DefineSecondaryStructures* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		DefineSecondaryStructures* m_obj;
	};
	friend class SetParam;
};

}
#endif
