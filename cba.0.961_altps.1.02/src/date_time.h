#ifndef CBA_DATE_TIME_H
#define CBA_DATE_TIME_H

//#include <vector>
#include <string>

namespace Cba {

class DateTime {
public:
	DateTime();
	~DateTime();

	static std::string current_local_time_str();
	static const std::string LOCAL_TIME_FORMAT;

private:

};

}
#endif
