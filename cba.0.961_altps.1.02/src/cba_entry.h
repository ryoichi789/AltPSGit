#ifndef CBA_CBA_ENTRY_H
#define CBA_CBA_ENTRY_H

#include <string>
#include <iostream>
#include "atom.h"
#include "molecule.h"
#include "protein.h"

namespace Cba
{
////////////////////////////////////////////////////////////////////////
// CBA format
////////////////////////////////////////////////////////////////////////

// * An entry consists of records of several types described below.
// * An entry starts with the "I" record (entity id).
// * A record is a single line.
// * A record (line) starts with a single-letter record identifier.
// * Fields in a record are delimited by tabs.
// * When two or more entries are contained in a file, the end of an entry
//  is indicated by a blank line or a line beginning with the character '/'.

// --------------------------------
// record "I": (1 record per entry; mandatory)
// --------------------------------
//	field: entity id
// e.g.,
// I	4HHBA_lbr_0

// --------------------------------
// record "N": (1)
// --------------------------------
//	field: entity name
// e.g.,
// N	HEMOGLOBIN (DEOXY)

// --------------------------------
// record "T": (1)
// --------------------------------
//	field: type of entity; value = {protein, peptide, dna, rna, small, water, atom, undefined}
// e.g.,
// T	protein

// --------------------------------
// record "A": (number of atoms)
// --------------------------------
//	field:
//		serial number,
//		serial number in the original database,
//		atom name in the original database,
//		element symbol,
//		x coordinate,
//		y coordinate,
//		z coordinate,
//		solvent-accessible surface area,
//		Van der Waals radius,
//		physico-chemical class (see atom.h),
//		physico-chemical type (see atom.h),
//		physico-chemical properties (see atom.h)
// e.g.,
// A	7   	307 	 CA 	C	    4.3980	   -1.0150	  -17.8610	    6.2708	2.00	6	7	0000000000000000

// --------------------------------
// record "B": (number of bonds)
// --------------------------------
//	field:
//		serial number of one atom,
//		serial number of the other atom,
//		valence of the bond (1=single, 2=double, ...)
// e.g.,
// B	2	7	1

// --------------------------------
// record "R": (number of residues)
// --------------------------------
//	field:
//		serial number,
//		serial number (of type string) in the original database,
//		residue name in the original database,
//		serial number of the last atom of the residue
// e.g.,
// R	7   	58  	HIS	20

// --------------------------------
// record "C":
// --------------------------------
//	field:
//		comments in any format in a line

class CbaEntry {
public:
	CbaEntry();
	~CbaEntry();

	const std::string& id() const;
	const std::string& name() const;
	const std::string& type() const;
	bool is_protein() const;

	void clear();
	bool read(std::istream& is);
	bool create_molecule(Molecule& m) const;
	bool create_protein(Protein& p) const;

	int ncomments() const;
	const std::string& get_comment(size_t i) const;

protected:
	struct BRecord {
		size_t a1;
		size_t a2;
		unsigned int v;
	};

	struct RRecord {
		size_t serial_num;
		std::string serial_num_orig;
		std::string name;
		size_t end_atom;
	};

private:
	std::string m_id;
	std::string m_name;
	std::string m_type;
	std::vector<Atom> m_atoms;
	std::vector<BRecord> m_bonds;
	std::vector<RRecord> m_residues;
	std::vector<std::string> m_comments;

	void read_i(const std::string& line);
	void read_n(const std::string& line);
	void read_t(const std::string& line);
	void read_a(const std::string& line);
	void read_b(const std::string& line);
	void read_r(const std::string& line);
	void read_c(const std::string& line);
};

inline const std::string& CbaEntry::id() const { return m_id; }
inline const std::string& CbaEntry::name() const { return m_name; }
inline const std::string& CbaEntry::type() const { return m_type; }
inline bool CbaEntry::is_protein() const { return m_type == "protein"; }
inline int CbaEntry::ncomments() const { return m_comments.size(); }
inline const std::string& CbaEntry::get_comment(size_t i) const { return m_comments.at(i); }
}
#endif
