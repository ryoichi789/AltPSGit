#include "app.h"
using namespace Cba;

App* app;
App::App() : m_debug(false) { app = this; }
App::~App() { app = 0; }

// parse_argv:
void App::parse_argv(char* argv[]) throw(BadArgv)
{
	m_prog = *argv;

	Param param;
	while (*++argv) {
		if (**argv == '-') {
			param.key = ++(*argv);
			if (*++argv) {
				param.value = *argv;
			} else {
				throw BadArgv();
			}
		} else {
			throw BadArgv();
		}
		m_params.push_back(param);
	}
}
