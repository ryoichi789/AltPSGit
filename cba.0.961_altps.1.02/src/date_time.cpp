//#include <vector>
#include <string>
//#include <iostream>
//#include <fstream>
//#include <algorithm>
#include <ctime>
#include "date_time.h"

using namespace Cba;
using std::string;

const std::string DateTime::LOCAL_TIME_FORMAT = "YYYYMMDDhhmmss";

//
DateTime::DateTime() {}

//
DateTime::~DateTime() {}

//
string DateTime::current_local_time_str()
{
	char buf[LOCAL_TIME_FORMAT.length() + 1];
	time_t t = std::time(0);
	std::strftime(buf, LOCAL_TIME_FORMAT.length()+1, "%4Y%2m%2d%2H%2M%2S", localtime(&t));
	return string(buf);
}
