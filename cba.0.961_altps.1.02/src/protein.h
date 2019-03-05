#ifndef CBA_PROTEIN_H
#define CBA_PROTEIN_H

#include <string>
#include <vector>
#include "molecule.h"
#include "atom.h"
#include "bio_sequence.h"
#include "exception.h"

namespace Cba {

class SecondaryStructure;

class Protein : public Molecule {
public:
	Protein();
	~Protein();
	void clear();
	//int nbonds() const { return 0; }
	//int bond_between(const Atom* a1, const Atom* a2) const { return 0; }
	Protein* clone(Molecule* protein) const;

	Protein(const Protein& prot);
	Protein& operator=(const Protein& prot);

	struct Residue {
		Residue() : n(0), h(0), ca(0), c(0), o(0), oxt(0), resseq(0) {}
		Atom* n; // main-chain N
		Atom* h; // main-chain H of NH
		Atom* ca; // main-chain C-alpha
		Atom* c; // main-chain C'
		Atom* o; // main-chain O
		Atom* oxt; // C-terminal main-chain O
		Atom* cb; // side-chain C-beta
		std::vector<Atom*> side_chain_atoms;
			// side_chain_atoms except for C-beta (including hydrogens)
		std::string name; // PDB residue name
		std::string number; // PDB residue number with insertion code as in PDB file
		int resseq; // PDB residue number as in PDB file
		int natoms() const;
		void clear();
	};

	int nresidues() const;
	const Residue& residue(size_t i) const;
	Residue& residue(size_t i);
	void add_residue(const Residue& residue);
	char residue_type(size_t i) const;
	void get_seq(BioSequence& seq) const;

	void remove_hydrogens();
	void assign_atom_types();
	char secondary_structure_at(size_t i) const;

	void print(std::ostream& os) const;
	void print_in_sdf(std::ostream& os) const {}

	int residue_containing_atom(const Atom* a) const;
		// returns residue number (serial) containing this atom;
		// returns -1 if no residue has the atom;

	//
	// ligand binding regions
	//
	struct LigandBindingRegion {
		std::vector<Atom*> atoms;
		bool predicted; // true if this is a prediction
		void clear() { atoms.clear(); }
	};
	int nlbrs() const; // number of ligand binding regions
	const LigandBindingRegion& get_lbr(size_t i) const;

private:
	std::vector<Residue> m_residues;
	std::string m_secondary_structures;
	std::vector<LigandBindingRegion> m_lbrs;

	void clone_residue(const Residue& r, Protein& p) const;
	void print_residues(std::ostream& os) const;

	// functions used by assign_atom_types()
	void set_type_for_n(Atom* a);
	void set_type_for_ca(Atom* a);
	void set_type_for_c(Atom* a);
	void set_type_for_o(Atom* a);
	void set_type_for_oxt(Atom* a);
	void set_type_for_cb(Atom* a);
	void set_cys_atom_property_type(std::vector<Atom*>& atoms);
	void set_asp_atom_property_type(std::vector<Atom*>& atoms);
	void set_glu_atom_property_type(std::vector<Atom*>& atoms);
	void set_phe_atom_property_type(std::vector<Atom*>& atoms);
	void set_his_atom_property_type(std::vector<Atom*>& atoms);
	void set_ile_atom_property_type(std::vector<Atom*>& atoms);
	void set_lys_atom_property_type(std::vector<Atom*>& atoms);
	void set_leu_atom_property_type(std::vector<Atom*>& atoms);
	void set_met_atom_property_type(std::vector<Atom*>& atoms);
	void set_asn_atom_property_type(std::vector<Atom*>& atoms);
	void set_pro_atom_property_type(std::vector<Atom*>& atoms);
	void set_gln_atom_property_type(std::vector<Atom*>& atoms);
	void set_arg_atom_property_type(std::vector<Atom*>& atoms);
	void set_ser_atom_property_type(std::vector<Atom*>& atoms);
	void set_thr_atom_property_type(std::vector<Atom*>& atoms);
	void set_val_atom_property_type(std::vector<Atom*>& atoms);
	void set_trp_atom_property_type(std::vector<Atom*>& atoms);
	void set_tyr_atom_property_type(std::vector<Atom*>& atoms);
	void set_ptr_atom_property_type(std::vector<Atom*>& atoms);
	void set_ace_atom_property_type(std::vector<Atom*>& atoms);
	void set_nh2_atom_property_type(std::vector<Atom*>& atoms);
	void set_unknown_property_type(std::vector<Atom*>& atoms);

	friend class SecondaryStructure;
};

inline int Protein::nresidues() const { return m_residues.size(); }
inline const Protein::Residue& Protein::residue(size_t i) const
{ return m_residues.at(i); }
inline Protein::Residue& Protein::residue(size_t i)
{ return m_residues.at(i); }
inline void Protein::add_residue(const Residue& residue)
{ m_residues.push_back(residue); }
inline char Protein::secondary_structure_at(size_t i) const
{ return m_secondary_structures.at(i); }
inline int Protein::nlbrs() const { return m_lbrs.size(); }
inline const Protein::LigandBindingRegion&
Protein::get_lbr(size_t i) const { return m_lbrs[i]; }

}
#endif
