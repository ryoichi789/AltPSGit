#ifndef CBA_SDF_ENTRY_H
#define CBA_SDF_ENTRY_H

#include <string>
#include <vector>
#include "molecule.h"

namespace Cba {

class SdfReader;

// SdfEntry:
// represents an entry from SD format (MOL format + annotations) file
// containing information on the structure of a (mainly small) molecule
//      see "CTfile formats" (August 2002) by MDL Information Systems, Inc.

class SdfEntry {
public:
	SdfEntry();
	const std::string& id() const;
	const std::string& name() const;
	int natoms() const;
	int nbonds() const;
	int nannots() const;
	const std::string& name_of_annot(size_t i) const;
	const std::vector<std::string>& annot(size_t i) const;
	bool create_molecule(Molecule& mol) const; // returns false in failure;

private:
	///
	/// Molfile block
	///
	struct Header {
		// line 1
		std::string molecule_name;

		// line 2
		std::string user_initials; // columns 0-1
		std::string program_name; // columns 2-9
		std::string date; // columns 10-19 MMDDYYHHmm
		unsigned char dimension; // columns 20-21 (e.g., '2D')

		// line 3
		std::string comments;

		Header() : dimension(0) {}
		void clear();
	};

	// Counts line contains
	//	number of atoms (columns 0-2) and
	//	number of bonds (columns 3-5).

	struct AtomRecord {
		double x; // columns 0-9
		double y; // columns 10-19
		double z; // columns 20-29
		std::string symbol; // columns 31-33; e.g., C, O, N, Cl, ...
		//char mass_difference; // columns 34-35; either of -3,-2,-1,0,1,2,3,4
		char charge; // columns 36-38;
		//char atom_stereo_parity; // columns 39-41
		//char hydrogen_count; // columns 42-44
		//char stereo_care_box; // columns 45-47
		//char valence; // columns 48-50
	};

	struct BondRecord {
		int atom1; // columns 0-2
		int atom2; // columns 3-5
		int bond_type; // columns 6-8;
			// 1=single, 2=double, 3=triple, 4=aromatic,
			// 5=single or double, 6=single or aromatic,
			// 7=double or aromatic, 8=any
		//char bond_stereo; // columns 9-11
	};

	///
	/// Annotation from SDfile
	///
	struct Annot {
		// header part
		std::string name; // enclosed in brackets; e.g., <melting.point>
		int registry_number;
			// external registry numbers are enclosed in parentheses;
			// e.g.,
			// >  <molecular.weight> (5452)
		std::vector<std::string> annot;
		// may extend over multiple lines
		// containing up to 200 characters each, terminated by a blank line 
		Annot() : registry_number(0) {}
	};

	Header m_header;
	std::string m_id;
	std::string m_name;
	int m_natoms;
	int m_nbonds;
	std::vector<AtomRecord> m_atoms;
	std::vector<BondRecord> m_bonds;
	std::vector<Annot> m_annots;

	void clear();
	friend class SdfReader;
};

inline const std::string& SdfEntry::id() const { return m_id; }
inline const std::string& SdfEntry::name() const { return m_name; }
inline int SdfEntry::natoms() const { return m_natoms; }
inline int SdfEntry::nbonds() const { return m_nbonds; }
inline int SdfEntry::nannots() const { return m_annots.size(); }

}
#endif
