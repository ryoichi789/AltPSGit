#ifndef CBA_MOL2_READER_H
#define CBA_MOL2_READER_H

#include <iostream>
#include <string>
#include "mol2_entry.h"

namespace Cba {

// assumes 1 entry in 1 file in MOL2 format

class Mol2Reader {
public:
	static bool read(Mol2Entry& ent, std::istream& from);

private:
	enum RecordType {
		ALT_TYPE = 0,
		ANCHOR_ATOM,
		ASSOCIATED_ANNOTATION,
		ATOM,
		BOND,
		CENTER_OF_MASS,
		CENTROID,
		COMMENT,
		CRYSIN,
		DICT,
		DATA_FILE,
		EXTENSION_POINT,
		FF_PBC,
		FFCON_ANGLE,
		FFCON_DIST,
		FFCON_MULTI,
		FFCON_RANGE,
		FFCON_TORSION,
		LINE,
		LSPLANE,
		MOLECULE,
		NORMAL,
		QSAR_ALIGN_RULE,
		RING_CLOSURE,
		ROTATABLE_BOND,
		SEARCH_DIST,
		SEARCH_OPTIONS,
		SET,
		SUBSTRUCTURE,
		U_FEAT,
		UNITY_ATOM_ATTR,
		UNITY_BOND_ATTR,
		NO_OF_RECORD_TYPES
	};

	static const char* record_names[];
	static int get_first_record_type(std::istream& from);
	static int get_record_type(const std::string& line);

	static int read_molecule_records(Mol2Entry& ent, std::istream& from);
	static int read_atom_records(Mol2Entry& ent, std::istream& from);
	static int read_bond_records(Mol2Entry& ent, std::istream& from);
	static int skip_records(std::istream& from);
		// these functions return next record type
	static void read_atom_record(Mol2Entry::AtomRecord& atom_record, const std::string& from);
	static void read_bond_record(Mol2Entry::BondRecord& bond_record, const std::string& from);
};

}
#endif
