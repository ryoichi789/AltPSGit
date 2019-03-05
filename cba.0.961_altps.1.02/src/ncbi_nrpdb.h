#ifndef CBA_NCBI_NRPDB_H
#define CBA_NCBI_NRPDB_H

#include <vector>
#include <map>
#include <string>
#include <iostream>

namespace Cba {

// Data contents of NCBI nrpdb
/*
#---------------------------------------------------------------------
# Non-redundant PDB chain set
#---------------------------------------------------------------------
#
# 1: PDB code
# 2: Chain ID
# 3: MMDB ID
#
# 4: Group ID (BLAST pvalue 10e-7)
# 5: Rank (BLAST pvalue 10e-7)
# 6: Representative (1) or not (0) (BLAST pvalue 10e-7)
#
# 7: Group ID (BLAST pvalue 10e-40)
# 8: Rank (BLAST pvalue 10e-40)
# 9: Representative (1) or not (0) (BLAST pvalue 10e-40)
#
# A: Group ID (BLAST pvalue 10e-80)
# B: Rank (BLAST pvalue 10e-80)
# C: Representative (1) or not (0) (BLAST pvalue 10e-40)
#
# D: Group ID (Non-identical sequences)
# E: Rank (Non-identical sequence)
# F: Representative (1) or not (0) (Non-identical sequence)
#
# G: Percentage of Unknown residues
# H: Percentage of Incomplete residues
# I: Percentage of Missing residues
# J: Percentage of Incomplete side-chain residues
# K: Resolution (0.0 if NMR)
# L: No. of chains (subunits) in the PDB entry
# M: No. of heterogens in the PDB entry
# N: No. of different heterogen types in the PDB entry
# O: No. of residues in the chain
# P: Method of coordinate determination (X: X-ray, N: NMR, M: theoretical)
# Q: Acceptable (a) in structural quality or not (n)
*/

class NcbiNrpdb {
public:
	NcbiNrpdb();
	~NcbiNrpdb();

	int read(std::istream& is); // returns number of etries
	int nentries() const;

	struct Entry;
	const Entry* get_entry(const std::string& id);
		// id = pdb-code + chain-id
		// returns 0 if no such entry
	const Entry* get_entry(const char* pdb_code, char chain_id);
		// returns 0 if no such entry
	const Entry* get_entry(int i) const;
		// returns i-th entry; 0 if no such entry

	struct Entry {
		std::string pdb_id;
		char chain_id;
		int mmdb_id;
		struct Group {
			int id;
			int rank;
			bool is_representative;
		};
		Group group_7; // group from p-value 10e-7
		Group group_40; // group from p-value 10e-40
		Group group_80; // group from p-value 10e-80
		Group group_ni; // group from non-identical
		float unknown_res; // Percentage of Unknown residues
		float incomplete_res; // Percentage of Incomplete residues
		float missing_res; // Percentage of Missing residues
		float incomplete_sc; // Percentage of Incomplete side-chain residues
		float resolution; // Resolution (0.0 if NMR)
		int nchains; // No. of chains (subunits) in the PDB entry
		int nheterogens; // No. of heterogens in the PDB entry
		int nhettypes; // No. of different heterogen types in the PDB entry
		int nresidues; // No. of residues in the chain
		char method; // Method of coordinate determination (X: X-ray, N: NMR, M: theoretical)
		bool acceptable; // Acceptable (a) in structural quality or not (n)
	};

private:
	std::vector<Entry> m_entries;
	std::map<std::string, Entry*> m_entry_map;
	std::string m_id_work;
};

inline int NcbiNrpdb::nentries() const { return m_entries.size(); }
}
#endif
