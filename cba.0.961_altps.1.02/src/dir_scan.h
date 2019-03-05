#ifndef CBA_DIR_SCAN_H
#define CBA_DIR_SCAN_H

// class DirScan:
// tranverses all files and directories in a given directory (root) and
// applies to each file the process_file() function of a FileProcessor object.

#include <string>

namespace Cba {

//
class FileProcessor {
public:
	FileProcessor();
	virtual ~FileProcessor();
	virtual void process_file(const std::string& file_path) = 0;
private:
};

//
class DirScan {
public:
	DirScan(const char* root, FileProcessor& file_processor);
	~DirScan();
	void scan();
private:
	std::string m_root;
	FileProcessor& m_file_processor;
	void scan(const std::string& dir);
};

inline void DirScan::scan() { scan(m_root); }

}
#endif
