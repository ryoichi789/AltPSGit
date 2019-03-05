#include "app_sample.h"
#include <cstdlib>
#include <numeric>
#include <iostream>

using namespace Cba;

AppSample the_app;

// set_param_default:
void AppSample::set_param_default() {}

// set_param:
void AppSample::set_param(char *argv[])
{
	while (*++argv)
		m_numbers.push_back(std::atoi(*argv));
}

// run:
int AppSample::run()
{
	std::cout << std::accumulate(m_numbers.begin(), m_numbers.end(), 0) << std::endl;
	return 0;
}

// usage:
void AppSample::usage() {}
