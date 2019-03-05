#ifndef CBA_ATOM_H
#define CBA_ATOM_H

#include <string>
#include "position.h"
#include "bit_string.h"

namespace Cba {

struct Atom {
	Atom();
	std::string name;
		// atom name as defined by the data source;
		// e.g., " CA ", " CB ", ..., from PDB
	unsigned char element; // element number (H=1, He=2,...,C=6,N=7,O=8,...)
	Position pos; // position (coordinates) in 3D space
	int charge;
	BitString properties; // physico-chemical properties of Bush & Sheridan (1993)
	unsigned char pc_class; // physico-chemical class of Bush & Sheridan (1993)
	int type() const; // // atom type of Miller et al. (1999)
	double vdwr() const; // Van der Waals radius
	double asa; // solvent-accessible surface area

	static const char* symbols[]; // atomic symbols in order of element numbers
	static unsigned char get_element_number(const char* symbol);
	static const double vdw_radii[]; // Van der Waals radii of atoms according to their types

	//////////////////////////////////////////////////////////////////////
	// physico-chemical properties of non-hydrogen atoms
	//	using Bush & Sheridan method (1993):
	//	Bush BL, Sheridan RP (1993) J Chem Inf Comput Sci 33:756-762.
	// uses the following sixteen physico-chemical properties:
	//   # on hybridization
	//       0=sp, 1=sp2, 2=sp3, 3=conj, 4=res, 
	//   # on number of neighbors (bonded atoms)
	//       5=x0, 6=x1, 7=x2, 8=x3, 9=x4,
	//   # on aromaticity (5-, 6-membered)
	//       10=ar, 11=ar5, 12=ar6, 
	//   # electronegativity 
	//       13=neg,
	//   # for amide nitrogen
	//       14=namide,
	//   # for carboxylate carbon
	//       15=cx
	//////////////////////////////////////////////////////////////////////
	// each atom has 1 or 0 for each of the 16 properties in 'properties'
	enum Properties {
		sp = 0, sp2 = 1, sp3 = 2, conj = 3, res = 4,
		x0 = 5, x1 = 6, x2 = 7, x3 = 8, x4 = 9,
		ar = 10, ar5 = 11, ar6 = 12,
		neg = 13, namide = 14, cx = 15, number_of_properties
	};

	///////////////////////////////////////////////////////////////////////////
	// physicochemical classes of atoms:
	//	1=cation, 2=anion, 3=donor, 4=acceptor, 5=polar, 6=hydrophobic, 7=none
	//
	//	according to:
	//	Bush BL, Sheridan RP (1993) J Chem Inf Comput Sci 33:756-762.
	///////////////////////////////////////////////////////////////////////////
	enum Classes {
		UNDEFINED = 0,
		CATION = 1, ANION = 2, DONOR = 3, ACCEPTOR = 4,
		POLAR = 5, HYDROPHOBIC = 6, NONE = 7, NUMBER_OF_ATOM_PC_CLASSES
	};

	//////////////////////////////////////////////////////////////////////
	// atom types:
	// according to:
	// Miller MD, Sheridan RP, Kearsley SK (1999) J Med Chem 42: 1505-1514.
	//
	// 44 = 43 types + 1 (unknown)
	//////////////////////////////////////////////////////////////////////
	enum Types {
		UNKNOWN  = 0,
		Br_sp3_6 = 1,
		Cl_sp3_6 = 2,
		C_sp1_6  = 3, C_sp1_7  = 4,
		C_sp2_6  = 5, C_sp2_7  = 6,
		C_sp3_6  = 7, C_sp3_7  = 8,
		F_sp3_6  = 9,
		I_sp3_6  = 10,
		N_sp1_4  = 11, N_sp1_7  = 12,
		N_sp2_1  = 13, N_sp2_2  = 14, N_sp2_3  = 15,
		N_sp2_4  = 16, N_sp2_5  = 17, N_sp2_7  = 18,
		N_sp3_1  = 19, N_sp3_2  = 20, N_sp3_3  = 21,
		N_sp3_5  = 22, N_sp3_7  = 23,
		O_sp2_1  = 24, O_sp2_2  = 25, O_sp2_4  = 26, O_sp2_5  = 27,
		O_sp3_1  = 28, O_sp3_2  = 29, O_sp3_4  = 30, O_sp3_5  = 31,
		P_sp2_7  = 32, P_sp3_1  = 33, P_sp3_7  = 34,
		S_sp1_7  = 35,
		S_sp2_2  = 36, S_sp2_4  = 37, S_sp2_6  = 38, S_sp2_7  = 39,
		S_sp3_1  = 40, S_sp3_2  = 41, S_sp3_6  = 42, S_sp3_7  = 43
	};

	int number; // serial number in a molecule
};

inline double Atom::vdwr() const { return vdw_radii[type()]; }

}
#endif
