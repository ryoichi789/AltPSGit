#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include "local_time.h"
#include "date_time.h"

using namespace Cba;

LocalTime the_app;

//
LocalTime::LocalTime() {}

//
void LocalTime::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f]"
		<< " (f by default)" << endl;
}

//
void LocalTime::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter "
			<< param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

//
void LocalTime::set_param_default()
{
	m_prog = "local_time";
}

//
void LocalTime::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
}

//
int LocalTime::run()
{
	std::cout << DateTime::current_local_time_str() << std::endl;
	return 0;
}
