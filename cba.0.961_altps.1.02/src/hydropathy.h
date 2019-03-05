#ifndef CBA_HYDROPATHY_H
#define CBA_HYDROPATHY_H

#include <map>

namespace Cba {

class Hydropathy {
public:
	Hydropathy();
	double value_for(char a) const;
private:
	std::map<char, double> m_values;
	void set_values();
	static const double VALUE_FOR_UNKNOWN_AMINO_ACID = 0.0;
};

}
#endif
