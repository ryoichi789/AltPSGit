#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "ncbi_nrpdb.h"
#include "exception.h"

using namespace Cba;
using std::vector;
using std::string;
using std::map;

//
NcbiNrpdb::NcbiNrpdb() {}

//
NcbiNrpdb::~NcbiNrpdb() {}

//
int NcbiNrpdb::read(std::istream& is)
{
	m_entries.clear();
	m_entry_map.clear();

	Entry ent;
	char c;
	string s;
	while (is.get(c)) {
		if (c == '#') {
			getline(is, s);
		} else {
			ent.pdb_id.erase();
			ent.pdb_id += c;
			for (int i = 1; i < 4; i++) {
				is.get(c);
				ent.pdb_id += c;
			}
			is.get(c); // skip one space
			is.get(ent.chain_id);
			is >> ent.mmdb_id
				>> ent.group_7.id >> ent.group_7.rank >> ent.group_7.is_representative
				>> ent.group_40.id >> ent.group_40.rank >> ent.group_40.is_representative
				>> ent.group_80.id >> ent.group_80.rank >> ent.group_80.is_representative
				>> ent.group_ni.id >> ent.group_ni.rank >> ent.group_ni.is_representative
				>> ent.unknown_res >> ent.incomplete_res >> ent.missing_res
				>> ent.incomplete_sc >> ent.resolution
				>> ent.nchains >> ent.nheterogens >> ent.nhettypes >> ent.nresidues
				>> ent.method >> c;
			if (c == 'a') ent.acceptable = true;
			else ent.acceptable = false;

			is.get(c); // end of line in Unix
			if (c == '\r') is.get(c); // end of line in Windows
			if (c != '\n') throw Exception("input failed");
			m_entries.push_back(ent);
			s = ent.pdb_id;
			s += ent.chain_id;
			m_entry_map[s] = &(m_entries.back());
		}
	}

	return m_entries.size();
}

//
const NcbiNrpdb::Entry* NcbiNrpdb::get_entry(int i) const
{
	if (i < 0 || i > nentries()) return 0;
	return &(m_entries[i]);
}

//
const NcbiNrpdb::Entry* NcbiNrpdb::get_entry(const string& id)
{
	if (id.length() < 5) return 0;
	m_id_work = id;
	m_id_work.erase(5);

	map<string, Entry*>::const_iterator p = m_entry_map.find(m_id_work);
	if (p != m_entry_map.end()) {
		return p->second;
	} else if (m_id_work.at(4) == '_') {
		m_id_work.at(4) = ' ';
		p = m_entry_map.find(m_id_work);
		if (p != m_entry_map.end()) return p->second;
	}
	return 0;
}

//
const NcbiNrpdb::Entry* NcbiNrpdb::get_entry(const char* pdb_code, char chain_id)
{
	m_id_work = pdb_code;
	m_id_work += chain_id;
	return get_entry(m_id_work);
}
