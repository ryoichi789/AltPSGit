#ifndef CBA_TEST_SUPERPOSER_H
#define CBA_TEST_SUPERPOSER_H

#include <functional>
#include "app.h"
#include "molecule.h"
#include "superposer.h"

namespace Cba {

class TestSuperposer : public App {
public:
	TestSuperposer();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	const char* m_molfile;
	bool m_sdf;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(TestSuperposer* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		TestSuperposer* m_obj;
	};
	friend class SetParam;

	void superpose_molecules(
		Molecule& mol1, Molecule& mol2, Superposer& spp) const;
};

}
#endif
