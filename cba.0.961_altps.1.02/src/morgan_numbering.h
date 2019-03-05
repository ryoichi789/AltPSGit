#ifndef CBA_MORGAN_NUMBERING_H
#define CBA_MORGAN_NUMBERING_H

#include <functional>
#include "app.h"

namespace Cba {

class MorganNumbering : public App {
public:
	MorganNumbering();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();

	const char* m_fin;
	const char* m_fout;
	bool m_count_hydrogens;
	bool m_verbose;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(MorganNumbering* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		MorganNumbering* m_obj;
	};
	friend class SetParam;
};

}
#endif
