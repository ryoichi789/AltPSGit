#ifndef CBA_PDB_FILE_PROCESSOR_H
#define CBA_PDB_FILE_PROCESSOR_H

// define a derived class from PdbFileProcessor to do some specific
// operation to each pdb file;
// define a specific operation by overloading the function process_file()
// in this way:
//
//      YourPdbFileProcessor::YourPdbFileProcessor(const char* pdb_file_ext)
//          : PdbFileProcessor(pdb_file_ext) {}
//
//      void YourPdbFileProcessor::process_file(const std::string& file_path)
//      {
//          if (! is_pdb_file(file_path)) return;
//          // do some operations ...
//          // ... ... ...
//      }
//
// apply an object of the derived class to a DirScan object (see dir_scan.h);

#include <string>
#include "dir_scan.h"

namespace Cba {

class PdbFileProcessor : public FileProcessor {
public:
	PdbFileProcessor(const char* pdb_file_ext = "pdb");
	virtual ~PdbFileProcessor();
	virtual void process_file(const std::string& file_path) = 0;
	bool is_pdb_file(const std::string& file_path) const;
private:
	std::string m_pdb_file_ext;
	int m_len_pdb_file_ext;
};

}
#endif
