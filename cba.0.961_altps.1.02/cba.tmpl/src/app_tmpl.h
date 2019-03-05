///
/// template for header (.h) file of a new application
///
/// NOTE: replace APP_TMPL by YOUR_APPLICATION
/// NOTE: replace 'AppTmpl' by 'YourApplication'
#ifndef CBA_APP_TMPL_H
#define CBA_APP_TMPL_H

#include <functional>
#include "app.h"
/// NOTE: include headers necessary for YourApplication

namespace Cba {

class AppTmpl : public App {
public:
	AppTmpl();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();

	/// NOTE:
	/// replace 'm_param1' and 'm_param2' below
	/// by actual parameters of YourApplication
	int m_param1;
	double m_param2;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(AppTmpl* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		AppTmpl* m_obj;
	};
	friend class SetParam;
};

}
#endif
