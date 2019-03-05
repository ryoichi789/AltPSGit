#include <string>
#include "pdb_file_processor.h"

using namespace Cba;

//
PdbFileProcessor::PdbFileProcessor(const char* pdb_file_ext)
{
	m_pdb_file_ext = ".";
	m_pdb_file_ext += pdb_file_ext;
	m_len_pdb_file_ext = m_pdb_file_ext.length();
}

//
PdbFileProcessor::~PdbFileProcessor() {}

//
bool PdbFileProcessor::is_pdb_file(const std::string& file_path) const
{
	int l = file_path.length();
	if (l <= m_len_pdb_file_ext) return false;
	if (0 == file_path.compare(l-m_len_pdb_file_ext, m_len_pdb_file_ext,
		m_pdb_file_ext, 0, m_len_pdb_file_ext)) return true;
	return false;
}
