#ifndef CBA_MOLECULE_H
#define CBA_MOLECULE_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "atom.h"
#include "matrix.h"
#include "asa_calculator.h"

namespace Cba {

class Molecule {
public:
	Molecule();
	virtual ~Molecule();
	virtual void clear();

	Molecule(const Molecule& mol);
	Molecule& operator=(const Molecule& mol);

	std::string id() const;
	std::string name() const;

	int natoms() const;
	Position position_of_atom(size_t i) const;
	const Atom* atom(size_t i) const;
	Atom* atom(size_t i);
	void get_atoms_bonded_to(const Atom* a, std::vector<Atom*>& as);
	void get_atoms_bonded_to(size_t i, std::vector<int>& as);

	int get_graph_distance(int a1, int a2);
	int get_graph_distance(const Atom* a1, const Atom* a2);
		// graph_distance = number of bonds in the shortest path between the two atoms;
		// returns -1 if the atoms are not connected;

	bool is_identical_to(Molecule& mol);
		// returns true if this molecule is identical to mol;
	bool is_identical_to(Molecule& mol, std::vector<std::vector<int> >& alignments);
		// alignments[i][j] is the serial number of the atom of the molecule 'mol'
		//  equivalent to the j-th atom of this molecule in the i-th alignment

	bool is_a_complex();
		// returns true if this is actually a complex of two or more molecules
	int get_component_molecules(std::vector<Molecule>& mols);
		// NOTE: use this function only if is_a_complex() is true;
		// creates Molecule objects for component molecules;
		// returns the number of component molecules;

	virtual int nbonds() const;
	virtual int bond_between(int a1, int a2) const;
	virtual int bond_between(const Atom* a1, const Atom* a2) const;
		// returns 0 (= no bond), 1 (= single bond), 2 (= double bond), etc.;
		// nbonds() and bond_between() are implemented differently
		// for protein and nucleic acid (DNA and RNA).
		// a1, a2 should be values from function atom(size_t i)
	struct BondBetweenAtoms {
		int atom1;
		int atom2;
		int v; // 1 = single, 2 = double, 3 = triple bonds
	};
	virtual void get_bond(size_t i, BondBetweenAtoms& bond) const;
	virtual int get_bonds(std::vector<BondBetweenAtoms>& bonds) const;
		// returns number of bonds
	virtual int get_rotatable_bonds(std::vector<int>& rotatable_bonds);
		// returns number of rotatable bonds;
		// bonds are identified by their serial numbers (0,1,2,...)
		//   as stored in the vector 'bonds' returned from get_bonds()
	virtual int n_rotatable_bonds();

	virtual int get_rings(std::vector<std::vector<int> >& rings);
	virtual int get_rings(std::vector<std::vector<Atom*> >& rings);
		// returns number of 3- to 8-membered rings
	virtual int get_aromatic_rings(std::vector<std::vector<int> >& rings);
	virtual int get_aromatic_rings(std::vector<std::vector<Atom*> >& rings);
		// returns number of aromatic rings

	virtual int nfragments();
	virtual void get_fragment(size_t i, std::vector<int>& fragment);
	virtual void get_fragment(size_t i, std::vector<Atom*>& fragment);
	virtual void get_fragment(size_t i, Molecule& fragment);
		// store the atoms and bonds of the fragment into a Molecule object
	virtual int get_fragments(std::vector<std::vector<int> >& fragments);
	virtual int get_fragments(std::vector<std::vector<Atom*> >& fragments);
		// fragment is a part of molecule that contains no rotatable bond.

	virtual Molecule* clone(Molecule* mol) const;
	// deep copy;
	// if molecule (argument) == 0, create new Molecule object;
	// returns 0 in error;

	virtual bool clone(Molecule& mol) const;
		// deep copy;
		// returns false in error;

	void shift_by(double x, double y, double z);
	void shift_by(const Position& p);
	void rotate_by(const Matrix& rot); // rot must be a 3x3 rotation matrix

	enum Type {
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
	Type type() const;

	void set_id(const std::string& id);
	void set_name(const std::string& name);
	void set_type(Type t);
	void add_atom(const Atom& a);
	void add_bond(Atom* a1, Atom* a2, int v);
		// a1, a2 should be values from function atom(size_t i)

	virtual void remove_hydrogens();
	virtual void assign_atom_types();
		// using Bush & Sheridan (1993) J.Chem.Inf.Comput.Sci 33:756-762;
		// hydrogen atoms are removed from the molecule;
		// Van der Waals radii are also assigned according to the types;
		// for macromolecules from PDB, this function is implemented differently; 
	virtual void remove_atom(const Atom* a);
		// NOTE: number of atoms will change after remove_atom();
		// do not use, or be careful in using a 'for' loop to remove atoms

	void calc_atom_asas();
	void calc_atom_asas(double rsolv, double zstep);

	// atom numbering by Morgan algorithm
	void get_morgan_numbering(std::vector<int>& ns);
		// ns: atom numbering (order of atoms)
	void get_morgan_numbering(std::vector<int>& ns, std::vector<std::vector<int> >& hist);
		// hist: numbers assigned to atoms incrementally during process + number of numbers

	/*
	 * for printing
	 */
	virtual void print(std::ostream& os) const; // print in CBA format
	void print_id(std::ostream& os) const;
	void print_name(std::ostream& os) const;
	void print_type(std::ostream& os) const;
	void print_atoms(std::ostream& os) const;
	void print_atom(std::ostream& os, size_t i) const;
	void print_bonds(std::ostream& os) const;
	void print_bond(std::ostream& os, size_t i) const;
	void print_entry_terminator(std::ostream& os) const;

	virtual void print_in_sdf(std::ostream& os, bool annot = false); // print in SDF format
		// NOTE: if annot = true, assign_atom_types() are executed, and all hydrogens are removed.
	void print_in_sdf_header_line2(std::ostream& os) const;
	void print_in_sdf_count_line(std::ostream& os) const;
	void print_in_sdf_atom(std::ostream& os, const Atom* atom) const;
	void print_in_sdf_bond(std::ostream& os, int i, int j, int t) const;
	void print_in_sdf_annot(std::ostream& os);

protected:
	struct Bond { // represents a covalent bond between two atoms
		Atom* atom1;
		Atom* atom2;
		int v; // bond type as defined by SDF
	};
	void remove_bonds_to(const Atom* a);
	int atom_number_of(const Atom* a) const;
	int bond_number_of(const Bond* b) const;

	virtual int find_rings();
		// returns number of 3- to 8-membered rings
		// rings are stored in the class object for later use
	virtual int get_n_membered_rings(std::vector<std::vector<Atom*> >& rings, int n);
		// returns number of n-membered rings

	virtual int find_rotatable_bonds();
		// returns number of rotatable bonds
		// rotatable bonds are stored in the class object for later use

	virtual int find_fragments();

	bool m_atom_types_assigned;

private:
	std::string m_id;
	std::string m_name;
	Type m_type;
	std::vector<Atom*> m_atoms; // atoms to count
	std::vector<Bond*> m_bonds; // bonds to count
	std::vector<Atom*> m_all_atoms;
	std::vector<Bond*> m_all_bonds;

	std::vector<std::vector<int> > m_rings; // 3- to 8-membered rings
	std::vector<int> m_rotatable_bonds;
	bool is_ring_bond(const Bond* b) const;
	bool is_ring_bond(int b) const;

	std::vector<std::vector<int> > m_fragments;
	void set_fragment_id(int i, Molecule& fragment);
	void set_fragment_name(int i, Molecule& fragment);
	void set_fragment_atoms(const std::vector<int>& atoms, Molecule& fragment);
	void set_fragment_bonds(const std::vector<int>& atoms, Molecule& fragment);

	std::vector<std::vector<int> > m_graph_distances;
	void calc_graph_distances();
	void initialize_graph_distances();
	void get_linked_atoms(std::vector<std::vector<int> >& linked_atoms);
	void calc_graph_distances_from(int atomi, const std::vector<std::vector<int> >& linked_atoms);

	void assign_atom_properties();
	void assign_atom_pc_classes();

	// data and functions used by assign_atom_properties()
	std::map<Atom*, std::vector<Bond*> > m_bonds_to_atom;
	void define_bonds_to_atom();
	void clear_bonds_to_atom();
	void get_atoms_bonded_to_atom(Atom* a, std::vector<Atom*>& as);
	void clear_atom_properties();

	void assign_atom_sp_sp2_res();
	void assign_atom_conj();
	void assign_atom_sp3_res();
	void assign_atom_nneighbors();
	void assign_atom_aromaticity();
	void get_reachable_atoms_from(Atom* a, int nsteps,
		std::vector<std::vector<Atom*> >& paths);
	void visit_atom(
		Atom* a, int nsteps, int msteps,
		std::vector<std::vector<Atom*> >& paths, std::vector<Atom*>& path);
	void assign_atom_electronegativity();
	void assign_atom_amide();
	void assign_atom_carboxylate();

	void assign_atom_pc_class_default();
	void assign_atom_special_nitrogen_types();
	void assign_atom_quaternary_nitrogen();
	void assign_atom_sp3_amines();
	void assign_atom_conjugated_ns();
	void assign_atom_basic_conjugated_sp3_ns();
	void assign_atom_trisubstituted_sp2_amines();
	void assign_atom_disubstituted_sp_amines();
	void assign_atom_trisubstituted_os();
	void assign_atom_miscellaneous_anions();
	void assign_atom_ma_conjugated_sulhydryls();
	void assign_atom_ma_carboxylates_sequiv();
	void assign_atom_ma_phosphate_arsenate();
	void assign_atom_ma_muscimol();
	void assign_atom_ma_tetrazole();
	void assign_atom_diaminopyrimidine();
	void get_sb_nx1s(Atom* a0, std::vector<Atom*>& nx1s);
	void assign_atom_ns_in_imidazole();
	void assign_atom_sulfonamides_phosphonamides();
	void get_bonded_osx1s(Atom* a0, std::vector<Atom*>& osx1s);
	void assign_atom_carbonyl_hydroxide();
	void assign_atom_adjacent_keto_enols();
	void assign_atom_ake1();
	void assign_atom_ake2();
	void assign_atom_misc_buried();

	// used by is_identical_to()
	bool atom_counts_equal(const Molecule& mol) const;
	bool bond_counts_equal(const Molecule& mol) const;
	void get_equivalent_atom_candidates(Molecule& mol, std::vector<std::vector<int> >& equiv_candidates);
	void search_equivalent_atom_alignments(int iatom,
		const Molecule& mol, const std::vector<std::vector<int> >& equiv_candidates,
		std::vector<int>& alignment, std::vector<std::vector<int> >& alignments);

	// used by get_morgan_numbering()
	int get_morgan_hist(std::vector<std::vector<int> >& hist, const std::vector<std::vector<int> >& atom_neighbors);
	void do_morgan_numbering(std::vector<int>& ns, int nsteps, std::vector<std::vector<int> >& hist, std::vector<std::vector<int> >& atom_neighbors);
	void sort_atoms_in_morgan(int a, std::vector<int>& as, const std::vector<int>& ns);
};

inline std::string Molecule::id() const { return m_id; }
inline std::string Molecule::name() const { return m_name; }
inline Molecule::Type Molecule::type() const { return m_type; }
inline int Molecule::natoms() const { return m_atoms.size(); }
inline int Molecule::nbonds() const { return m_bonds.size(); }
inline int Molecule::nfragments() { return find_fragments(); }
inline int Molecule::n_rotatable_bonds() { return find_rotatable_bonds(); }
inline Position Molecule::position_of_atom(size_t i) const
{ return m_atoms[i]->pos; }
inline const Atom* Molecule::atom(size_t i) const { return m_atoms[i]; }
inline Atom* Molecule::atom(size_t i) { return m_atoms[i]; }
inline void Molecule::set_id(const std::string& id) { m_id = id; }
inline void Molecule::set_name(const std::string& name) { m_name = name; }
inline void Molecule::set_type(Type t) { m_type = t; }

}
#endif
