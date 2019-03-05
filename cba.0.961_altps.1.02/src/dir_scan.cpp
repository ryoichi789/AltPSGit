#include <string>
#include <cstring> //For CentOS6.2
#include <sys/types.h> // not C++ standard
#include <sys/stat.h> // not C++ standard
#include <dirent.h> // not C++ standard
#include "dir_scan.h"
 
using namespace Cba;
using std::string;

//
FileProcessor::FileProcessor() {}

//
FileProcessor::~FileProcessor() {}

//
DirScan::DirScan(const char* root, FileProcessor& file_processor)
	: m_root(root), m_file_processor(file_processor)
{}

//
DirScan::~DirScan() {}

//
void DirScan::scan(const string& dir)
{
	DIR* dp;
	struct dirent* ent;
	struct stat statbuf;
	std::string current_ent;

	if ((dp = opendir(dir.c_str())) == NULL) return;
	while ((ent = readdir(dp)) != NULL) {
		current_ent = dir + "/" + ent->d_name;
		stat(current_ent.c_str(), &statbuf);
		if (S_ISDIR(statbuf.st_mode)) {
			if (std::strcmp(".", ent->d_name) == 0 ||
				std::strcmp("..", ent->d_name) == 0) continue;
			scan(current_ent);
		} else {
			m_file_processor.process_file(current_ent);
		}
	}
	closedir(dp);
}
