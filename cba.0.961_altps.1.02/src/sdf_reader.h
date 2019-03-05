#ifndef CBA_SDF_READER_H
#define CBA_SDF_READER_H

#include <string>
#include <vector>
#include <iostream>
#include "sdf_entry.h"

namespace Cba {

class SdfReader {
public:
	SdfReader(std::istream& is);
	bool read(SdfEntry& ent);

private:
	std::istream& m_is;
	int m_serial_number; // used for id of molecule if no other id available

	bool read_header(SdfEntry& ent);
	bool read_counts_line(SdfEntry& ent);
	bool read_atoms(SdfEntry& ent);
	bool read_bonds(SdfEntry& ent);
	bool skip_properties_block();
	void read_annots(SdfEntry& ent);
	void set_id(SdfEntry& ent);
	void set_name(SdfEntry& ent);
};

}
#endif
