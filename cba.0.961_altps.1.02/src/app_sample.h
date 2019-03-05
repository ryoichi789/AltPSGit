/*
 * app_sample:
 *	- sums up numbers (integers) in command line arguments
 */
#ifndef CBA_APP_SAMPLE_H
#define CBA_APP_SAMPLE_H

#include "app.h"
#include <vector>

namespace Cba {

class AppSample : public App {
public:
	void set_param_default();
	void set_param(char* argv[]);
	int run();
private:
	void usage();
	std::vector<int> m_numbers;
};

}
#endif
