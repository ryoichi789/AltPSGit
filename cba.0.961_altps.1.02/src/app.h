/*
 * class App: the base for all Cba applications
 */
#ifndef CBA_APP_H
#define CBA_APP_H

#include <list>
#include "exception.h"

namespace Cba {

class App {
public:
	App();
	virtual ~App();
	virtual void set_param_default() = 0;
	virtual void set_param(char* argv[]) = 0;
	virtual int run() = 0; // returns 0 if successful

protected:
	virtual void usage() = 0; // shows usage

	class BadArgv : public Exception {};
		// thrown by parse_argv() (see below) when argument list is bad

	struct Param { // each (key,value) pair from argument list
		char* key;
		char* value;
	};

	void parse_argv(char* argv[]) throw(BadArgv);
		// calls parse_argv() at the beginning of set_param()
		// when the argument list is assumed to follow the style:
		// 		-key value ...

	char* m_prog; // program name
	std::list<Param> m_params; // list of (key,value) pairs in argument list
	bool m_debug; // true if debug mode
};

}
#endif
