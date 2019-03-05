#include <map>
#include <cctype>
#include "hydropathy.h"

using namespace Cba;

//
Hydropathy::Hydropathy() { set_values(); }

//
double Hydropathy::value_for(char a) const
{
	std::map<char, double>::const_iterator p = m_values.find(std::toupper(a));
	if (p != m_values.end()) return p->second;
	return VALUE_FOR_UNKNOWN_AMINO_ACID;
}

//
void Hydropathy::set_values()
{
	// Kyte and Doolittle (1982) J.Mol.Biol. 157:105-132.
	m_values['A'] =  1.8;
	m_values['C'] =  2.5;
	m_values['D'] = -3.5;
	m_values['E'] = -3.5;
	m_values['F'] =  2.8;
	m_values['G'] = -0.4;
	m_values['H'] = -3.2;
	m_values['I'] =  4.5;
	m_values['K'] = -3.9;
	m_values['L'] =  3.8;
	m_values['M'] =  1.9;
	m_values['N'] = -3.5;
	m_values['P'] = -1.6;
	m_values['Q'] = -3.5;
	m_values['R'] = -4.5;
	m_values['S'] = -0.8;
	m_values['T'] = -0.7;
	m_values['V'] =  4.2;
	m_values['W'] = -0.9;
	m_values['Y'] = -1.3;
	m_values['B'] = -3.5;
	m_values['Z'] = -3.5;
}
