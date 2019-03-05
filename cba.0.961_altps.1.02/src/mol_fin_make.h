#ifndef CBA_MOL_FIN_MAKE_H
#define CBA_MOL_FIN_MAKE_H

#include <functional>
#include "app.h"

namespace Cba {

class MolFinMake : public App {
public:
	MolFinMake();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();

	const char* m_fin;
	const char* m_fout;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(MolFinMake* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		MolFinMake* m_obj;
	};
	friend class SetParam;
};

}
#endif
