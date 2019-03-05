#ifndef CBA_LOCAL_TIME_H
#define CBA_LOCAL_TIME_H

#include <functional>
#include "app.h"

namespace Cba {

class LocalTime : public App {
public:
	LocalTime();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(LocalTime* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		LocalTime* m_obj;
	};
	friend class SetParam;
};

}
#endif
