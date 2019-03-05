#ifndef CBA_PDB_ENTRY_H
#define CBA_PDB_ENTRY_H

#include <string>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <list>
#include <vector>
#include <map>
#include <istream>
#include <fstream>
#include "pdb_record.h"
#include "molecule.h"
#include "protein.h"
#include "dir_scan.h"
#include "bio_sequence.h"

namespace Cba
{
class PdbReader;

class PdbEntry {
public:
	PdbEntry();
	~PdbEntry();
	PdbEntry(const PdbEntry& ent);
	PdbEntry& operator=(const PdbEntry& ent);
	PdbEntry* clone(PdbEntry* ent) const;
		// deep copy;
		// if argument 'ent' == 0, create a new PdbEntry object;
		// returns 0 in error;

	int nrecords() const;
	void set_id(const std::string& id);
	const std::string& id() const;
	std::string title() const;
	std::string compound() const;
	int nmolecules() const;
	int nresidues(int imol) const; // no. of residues with coordinates
	int natoms(int imol) const;
	int natoms(int imol, int ires) const;
	char chainid(int imol) const;
	const std::string& resname(int imol, int ires) const;
	std::string resnum(int imol, int ires) const;
		// returns string of residue serial number + insertion code
	int resseq(int imol, int ires) const;
		// returns residue number as in PDB file
	enum MoleculeType {
		UNDEFINED = 0,
		PROTEIN,
		PEPTIDE,
		DNA,
		RNA,
		SMALL,
		WATER,
		ATOM,
		NO_OF_MOLECULE_TYPES
	};
	MoleculeType type_of_molecule(int imol) const;
	bool is_protein(int imol) const;
	bool is_peptide(int imol) const;
	bool is_small_molecule(int imol) const;
	int nmr_model(int imol) const;
		// returns NMR model number (0 if not an NMR model)
	float resolution() const;

	static const int MAX_LEN_PEPTIDE;
		// polypeptides with length <= MAX_LEN_PEPTIDE
		// are considered "PEPTIDE" (oligopeptide).

	const std::string& get_small_molecule_id(int imol) const;
	const std::string& get_small_molecule_name(const std::string& id) const;

	bool create_molecule(int imol, Molecule& mol) const;
	bool create_protein(int imol, Protein& protein) const;
		// returns true if successfully created

	int get_depdate() const;
		// returns date of deposition in the form of YYYYMMDD

	int nseqs() const; // number of SEQRES records
	char get_seq_chain_id(int i) const;
	int length_of_seqres(int i) const;
	bool get_protein_sequence(BioSequence& seq, char chain = ' ') const;
	bool get_ith_sequence(BioSequence& seq, int i) const;

private:
	std::vector<PdbRecord*> m_records;
	PdbHeader* m_header;
	PdbObslte* m_obslte;
	PdbTitle* m_title;
	PdbCaveat* m_caveat;
	PdbCompnd* m_compnd;
	PdbSource* m_source;
	PdbKeywds* m_keywds;
	PdbExpdta* m_expdta;
	PdbAuthor* m_author;
	std::vector<PdbRevdat*> m_revdats;
	PdbSprsde* m_sprsde;
	PdbJrnl* m_jrnl;
	PdbRemark1* m_remark1;
	PdbRemark2* m_remark2;
	PdbRemark3* m_remark3;
	std::vector<PdbRemarkN*> m_remarkns;
	std::vector<PdbDbref*> m_dbrefs;
	std::vector<PdbSeqadv*> m_seqadvs;
	std::vector<PdbSeqres*> m_seqreses;
	std::vector<PdbModres*> m_modreses;
	std::vector<PdbHet*> m_hets;
	std::vector<PdbHetnam*> m_hetnams;
	std::vector<PdbHetsyn*> m_hetsyns;
	std::vector<PdbFormul*> m_formuls;
	std::vector<PdbHelix*> m_helices;
	std::vector<PdbSheet*> m_sheets;
	std::vector<PdbTurn*> m_turns;
	std::vector<PdbSsbond*> m_ssbonds;
	std::vector<PdbLink*> m_links;
	std::vector<PdbHydbnd*> m_hydbnds;
	std::vector<PdbSltbrg*> m_sltbrgs;
	std::vector<PdbCispep*> m_cispeps;
	std::vector<PdbSite*> m_sites;
	PdbCryst1* m_cryst1;
	PdbOrigx* m_origx;
	PdbScale* m_scale;
	std::vector<PdbMtrix*> m_mtrices;
	std::vector<PdbTvect*> m_tvects;
	std::vector<PdbModel*> m_models;
	std::vector<PdbEndmdl*> m_endmdls;
	std::vector<PdbAtom*> m_atoms;
	std::vector<PdbSigatm*> m_sigatms;
	std::vector<PdbAnisou*> m_anisous;
	std::vector<PdbSiguij*> m_siguijs;
	std::vector<PdbTer*> m_ters;
	std::vector<PdbConect*> m_conects;
	PdbMaster* m_master;
	PdbEnd* m_end;
	std::vector<PdbUnknownRecord*> m_unknowns;

	std::string m_id;
	std::vector<int> m_residues; // serial number of last atom of each residue
	std::vector<int> m_molecules; // serial number of last residue of each molecule
	std::map<std::string, std::string> m_het_names;

	void clear();
	//unsigned char get_element_number(const std::string& atom_name) const;
	unsigned char get_element_number(const PdbAtom* atom) const;

	// these functions are used by PdbReader
	void make_entry();
	void assign_id();
	void define_residues();
	void define_molecules();
	typedef std::vector<PdbRecord*>::const_iterator RCI;
	bool end_of_molecule(int ia, RCI p) const;
	int residue_begin(int imol) const; // serial number of first residue of a molecule
	int atom_begin(int imol) const; // serial number of first atom of a molecule
	void make_het_name_list();
	bool covalently_bonded_to_next(int ir) const; // used in define_molecules()
	static const double MAX_COVALENT_BOND_LENGTH; // used in covalently_bonded_to_next()
	static const double OVERLAPPED_DISTANCE; // used in covalently_bonded_to_next()
	bool peptide_bonded_to_next(int ir) const; // used in define_molecules()
	static const double MAX_PEPTIDE_BOND_LENGTH; // used in peptide_bonded_to_next()
	bool is_het_molecule(int ia) const; // used in end_of_molecule()
	bool is_metal(const PdbAtom* atom) const; // used in covalently_bonded_to_next()
	void remove_altloc(); // used in make_entry()
	friend class PdbReader;
};

inline void PdbEntry::set_id(const std::string& id) { m_id = id; }
inline const std::string& PdbEntry::id() const { return m_id; }
inline int PdbEntry::nrecords() const { return m_records.size(); }
inline int PdbEntry::nmolecules() const { return m_molecules.size(); }
inline bool PdbEntry::is_protein(int imol) const { return PdbEntry::PROTEIN == type_of_molecule(imol); }
inline bool PdbEntry::is_peptide(int imol) const { return PdbEntry::PEPTIDE == type_of_molecule(imol); }
inline bool PdbEntry::is_small_molecule(int imol) const { return PdbEntry::SMALL == type_of_molecule(imol); }
inline const std::string& PdbEntry::get_small_molecule_id(int imol) const { return resname(imol, 0); }
inline int PdbEntry::nseqs() const { return m_seqreses.size(); }

}
#endif
