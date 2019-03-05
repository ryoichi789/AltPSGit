#include <cstdlib>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "position.h"

using namespace Cba;

//
const char* Atom::symbols[] = {
	"x",
	"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
	"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
	"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
	"Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
	"Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
	"Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
	"Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
	"Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
	"T1", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
	"Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
	"Md", "No", "Lr", 0 
};

//////////////////////////////////////////////////////////////////
// ** van der Waals radii of atom groups **
// hydrogens are ignored;
// according to:
// Linus Pauling "The Nature of the Chemical Bond, 3rd ed.",
//	Cornell University Press, Ithaca, N.Y., 1960.
//
// see p.529 of Gordon M. Barrow "Physical Chemistry, 4th ed.",
//	McGraw-Hill, 1979.
// cf. p.26 of Ooi, Oobatake, Nemethy, Scheraga (1987)
//	Proc Natl Acad Sci USA 84:3086-3090.
//////////////////////////////////////////////////////////////////
// vdw_radii[]:
const double Atom::vdw_radii[] = {
	1.7, // unknown

	1.95, // Br_sp3_6

	1.8, // Cl_sp3_6

	1.7, // C_sp1_6
	1.7, // C_sp1_7
	1.7, // C_sp2_6
	1.7, // C_sp2_7
	2.0, // C_sp3_6
	2.0, // C_sp3_7

	1.35, // F_sp3_6

	2.15, // I_sp3_6
	
	1.5, // N_sp1_4
	1.5, // N_sp1_7
	1.5, // N_sp2_1
	1.5, // N_sp2_2
	1.5, // N_sp2_3
	1.5, // N_sp2_4
	1.5, // N_sp2_5
	1.5, // N_sp2_7
	1.5, // N_sp3_1
	1.5, // N_sp3_2
	1.5, // N_sp3_3
	1.5, // N_sp3_5
	1.5, // N_sp3_7

	1.4, // O_sp2_1
	1.4, // O_sp2_2
	1.4, // O_sp2_4
	1.4, // O_sp2_5
	1.4, // O_sp3_1
	1.4, // O_sp3_2
	1.4, // O_sp3_4
	1.4, // O_sp3_5

	1.9, // P_sp2_7
	1.9, // P_sp3_1
	1.9, // P_sp3_7

	1.85, // S_sp1_7
	1.85, // S_sp2_2
	1.85, // S_sp2_4
	1.85, // S_sp2_6
	1.85, // S_sp2_7
	1.85, // S_sp3_1
	1.85, // S_sp3_2
	1.85, // S_sp3_6
	1.85 // S_sp3_7
};

//
Atom::Atom() : element(0), charge(0), asa(0.0), number(0)
{
	properties.init(Atom::number_of_properties);
	pc_class = Atom::UNDEFINED;
}

// get_element_number:
unsigned char Atom::get_element_number(const char* symbol)
{
	for (unsigned char i = 1; symbols[i]; i++)
		if (strcmp(symbol, symbols[i]) == 0) return i;
	return 0;
}

// type:
int Atom::type() const
{
	bool _sp1 = properties.get_bit(sp);
	bool _sp2 = properties.get_bit(sp2);
	bool _sp3 = properties.get_bit(sp3);
	if (element == get_element_number("Br")) {
		if (_sp3 && pc_class == 6) return Br_sp3_6;
	} else if (element == get_element_number("Cl")) {
		if (_sp3 && pc_class == 6) return Cl_sp3_6;
	} else if (element == get_element_number("C")) {
		if (_sp1 && pc_class == 6) return C_sp1_6;
		else if (_sp1 && pc_class == 7) return C_sp1_7;
		else if (_sp2 && pc_class == 6) return C_sp2_6;
		else if (_sp2 && pc_class == 7) return C_sp2_7;
		else if (_sp3 && pc_class == 6) return C_sp3_6;
		else if (_sp3 && pc_class == 7) return C_sp3_7;
	} else if (element == get_element_number("F")) {
		if (_sp3 && pc_class == 6) return F_sp3_6;
	} else if (element == get_element_number("I")) {
		if (_sp3 && pc_class == 6) return I_sp3_6;
	} else if (element == get_element_number("N")) {
		if (_sp1 && pc_class == 4) return N_sp1_4;
		else if (_sp1 && pc_class == 7) return N_sp1_7;
		else if (_sp2 && pc_class == 1) return N_sp2_1;
		else if (_sp2 && pc_class == 2) return N_sp2_2;
		else if (_sp2 && pc_class == 3) return N_sp2_3;
		else if (_sp2 && pc_class == 4) return N_sp2_4;
		else if (_sp2 && pc_class == 5) return N_sp2_5;
		else if (_sp2 && pc_class == 7) return N_sp2_7;
		else if (_sp3 && pc_class == 1) return N_sp3_1;
		else if (_sp3 && pc_class == 2) return N_sp3_2;
		else if (_sp3 && pc_class == 3) return N_sp3_3;
		else if (_sp3 && pc_class == 5) return N_sp3_5;
		else if (_sp3 && pc_class == 7) return N_sp3_7;
	} else if (element == get_element_number("O")) {
		if (_sp2 && pc_class == 1) return O_sp2_1;
		else if (_sp2 && pc_class == 2) return O_sp2_2;
		else if (_sp2 && pc_class == 4) return O_sp2_4;
		else if (_sp2 && pc_class == 5) return O_sp2_5;
		else if (_sp3 && pc_class == 1) return O_sp3_1;
		else if (_sp3 && pc_class == 2) return O_sp3_2;
		else if (_sp3 && pc_class == 4) return O_sp3_4;
		else if (_sp3 && pc_class == 5) return O_sp3_5;
	} else if (element == get_element_number("P")) {
		if (_sp2 && pc_class == 7) return P_sp2_7;
		else if (_sp3 && pc_class == 1) return P_sp3_1;
		else if (_sp3 && pc_class == 7) return P_sp3_7;
	} else if (element == get_element_number("S")) {
		if (_sp1 && pc_class == 7) return S_sp1_7;
		else if (_sp2 && pc_class == 2) return S_sp2_2;
		else if (_sp2 && pc_class == 4) return S_sp2_4;
		else if (_sp2 && pc_class == 6) return S_sp2_6;
		else if (_sp2 && pc_class == 7) return S_sp2_7;
		else if (_sp3 && pc_class == 1) return S_sp3_1;
		else if (_sp3 && pc_class == 2) return S_sp3_2;
		else if (_sp3 && pc_class == 6) return S_sp3_6;
		else if (_sp3 && pc_class == 7) return S_sp3_7;
	}
	return UNKNOWN;
}
