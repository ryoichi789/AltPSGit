#ifndef CBA_PDB_READER_H
#define CBA_PDB_READER_H

#include <iostream>
#include <string>
#include "pdb_record.h"
#include "pdb_entry.h"

namespace Cba {

class PdbReader {
public:
	static bool read(PdbEntry& ent, std::istream& from);

private:
	static void add_record(PdbEntry& ent, const std::string& line);
};

}
#endif
