///
/// copy this to 'your_application.cpp';
//  follow the NOTEs below
///
/// NOTE: replace 'app_tmpl' by 'your_application'
/// NOTE: replace 'AppTmpl' by 'YourApplication'
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
//#include <vector>
//#include <string>
//#include <list>
//#include <map>
//#include <set>
//#include <fstream>
//#include <iomanip>
//#include <cmath>
#include "app_tmpl.h"
// NOTE: include necessary headers

using namespace Cba;
//using std::vector;
//using std::string;
//using std::list;
//using std::map;
//using std::set;
//using std::istream;
//using std::ostream;
//using std::ifstream;
//using std::ofstream;

AppTmpl the_app;

// NOTE: define static const member variables here if any

//
AppTmpl::AppTmpl() {}

//
// NOTE: STEP 1: define proper command line parameters
void AppTmpl::usage()
{
	using std::cerr;
	using std::endl;
	// NOTE: modify the following lines
	//         to show usage specific to YourApplication
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\ti: integer value [integer]"
		<< " (unique; 0 by default)" << endl;
	cerr << "\td: double-floating-point value [double]"
		<< " (unique; 0.0 by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f]"
		<< " (f by default)" << endl;
}

//
// NOTE: STEP 2:
// 2-1. in the header file, define proper member variables
//      to hold command line parameters;
// 2-2. set values to those member variables
void AppTmpl::SetParam::operator()(const Param& param)
{
	// NOTE: modify the following lines
	//         to set parameters specific to YourApplication
	using std::strcmp;
	if (strcmp(param.key, "i") == 0) {
		m_obj->m_param1 = std::atoi(param.value);
	} else if (strcmp(param.key, "d") == 0) {
		m_obj->m_param2 = std::atof(param.value);
	} else if (strcmp(param.key, "debug") == 0) {
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
// NOTE: STEP 3: define default values of parameters
void AppTmpl::set_param_default()
{
	// NOTE:
	// replace these lines by definition of actual parameters;
	m_prog = "app_tmpl";
	m_param1 = 0;
	m_param2 = 0;
}

//
// NOTE: STEP 4:
// define what to do if parameters are not set properly
void AppTmpl::set_param(char *argv[])
{
	// NOTE:
	// command line arguments must be listed in <-key value> format
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	// NOTE:
	// check parameters here;
	// if they are invalid, exit the program
	// by calling std::exit(2);
	// e.g.,
	//if (m_param1 == 0) {
	//	usage();
	//	std::exit(2);
	//}
}

//
// NOTE: STEP 5: define the main routine
int AppTmpl::run()
{
	// NOTE: sample:
	//
	// == DO SOMETHING HERE ==
	//
	// e.g.,
	//std::cout << m_param1 << " " << m_param2 << std::endl;
	//if (m_debug)
	//	std::cerr << "running in debug mode..." << std::endl;
	//
	// In case of error, terminate the program by:
	// return 1;

	return 0;
}

// NOTE: define other member functions below
