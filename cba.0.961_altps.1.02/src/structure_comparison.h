#ifndef CBA_STRUCTURE_COMPARISON_H
#define CBA_STRUCTURE_COMPARISON_H

#include <vector>
#include <utility>
#include <functional>
#include <iostream>
#include "molecule.h"
#include "protein.h"
#include "atom.h"
#include "superposer.h"
#include "position.h"

namespace Cba {

class StructureComparison {
public:
	StructureComparison();
	~StructureComparison();
	int compare(Molecule* mol1, Molecule* mol2);
		// assumes assign_atom_types() is already done for mol1 and mol2;
		// returns number of atom (residue) alignments detected;

	struct AtomPair {
		int a1; // atom number in molecule 1
		int a2; // atom number in molecule 2
		AtomPair() {}
		AtomPair(int i, int j) : a1(i), a2(j) {}
	};

	struct AtomAlignment {
		std::vector<AtomPair> pairs;
		Superposer spp; // to superpose molecule2 onto molecule1
		double rmsd;
		int score;
	};

	struct ProteinAlignmentSegment {
		int r1;
		int r2;
		int length;
		ProteinAlignmentSegment() {}
		ProteinAlignmentSegment(int i, int j, int l)
			: r1(i), r2(j), length(l) {}
		bool operator==(const ProteinAlignmentSegment& s) const
		{ return r1 == s.r1 && r2 == s.r2 && length == s.length; }
		bool operator!=(const ProteinAlignmentSegment& s) const
		{ return r1 != s.r1 || r2 != s.r2 || length != s.length; }
		bool contains(const ProteinAlignmentSegment& s) const;
		bool overlaps(const ProteinAlignmentSegment& s) const;
		bool merges(const ProteinAlignmentSegment& s);
	};

	struct ProteinAlignment {
		std::vector<ProteinAlignmentSegment> segments;
		Superposer spp; // to superpose protein2 onto protein1
		double rmsd;
		void add_segment(int i, int j, int l);
		int length() const;
		void print(std::ostream& os) const;
		bool operator==(const ProteinAlignment& a) const;
		bool contains(const ProteinAlignment& a) const;
		bool overlaps(const ProteinAlignment& a) const;
		bool merges(const ProteinAlignment& a);
	};

	int n_atom_alignments() const;
	void get_atom_alignment(size_t i, AtomAlignment& a) const;

	int n_protein_alignments() const;
	void get_protein_alignment(size_t i, ProteinAlignment& a) const;

private:
	Molecule* m_mol1;
	Molecule* m_mol2;
	Protein* m_prot1;
	Protein* m_prot2;

	///////////////////////////////////
	// for small molecule comparison
	///////////////////////////////////
	static const int MIN_ATOM_PAIRS;
	static const double MAX_DIST_ATOM_PAIR;
	static const int MIN_ATOM_ALIGNMENT_SCORE;

	struct AtomTriad {
		Atom* a0; // center of triad
		Atom* a1;
		Atom* a2;
		AtomTriad() {}
		AtomTriad(Atom* a, Atom* b, Atom* c)
			: a0(a), a1(b), a2(c) {}
	};

	std::vector<AtomAlignment> m_atom_alignments;
	int compare();
	void get_atom_positions(const Molecule* mol,
		std::vector<Position>& ps) const;
	void get_atom_triads(std::vector<AtomTriad>& triads1,
		std::vector<AtomTriad>& triads2);
	void get_atom_triads(Molecule& mol,
		std::vector<AtomTriad>& triads);
	bool atom_triads_already_aligned(const AtomTriad& t1,
		const AtomTriad& t2) const;
	bool match_atom_triads(const AtomTriad& t1,
		const AtomTriad& t2) const;
	bool match_atoms(const Atom* a1, const Atom* a2) const;
	void superpose_atom_triads(
		const AtomTriad& t1, const AtomTriad& t2,
		AtomAlignment& a);
	void align_atoms(std::vector<Position>& ps1,
		std::vector<Position>& ps2, AtomAlignment& a);
	int get_equiv_atom(const Position& p,
		const std::vector<Position>& ps) const;
	void get_new_superposition(
		const std::vector<Position>& ps1,
		const std::vector<Position>& ps2, 
	    AtomAlignment& a);
	void add_atom_alignment(
		const std::vector<Position>& ps1,
		const std::vector<Position>& ps2, 
	    AtomAlignment& a);
	void get_atom_alignment_score(AtomAlignment& a) const;

	//////////////////////////////
	// for protein comparison
	//////////////////////////////
	static const int LENGTH_OF_SEGMENT;
	static const int MIN_LEN_ALIGNMENT_SEGMENT;
	static const int MIN_LEN_ALIGNMENT;
	static const double MAX_INIT_RMSD_PROTEIN;
	static const double MAX_DIST_RESIDUE_PAIR;

	std::vector<ProteinAlignment> m_protein_alignments;

	struct CmpProteinAlignmentSegment : std::binary_function
			<ProteinAlignmentSegment, ProteinAlignmentSegment, bool> {
		bool operator()
			(const ProteinAlignmentSegment& s1,
				const ProteinAlignmentSegment& s2) const
		{ return s1.r1 < s2.r1; }
	};

	struct CmpProteinAlignments
		: std::binary_function<ProteinAlignment, ProteinAlignment, bool> {
		bool operator()
			(const ProteinAlignment& a1, const ProteinAlignment& a2) const
		{ return a1.length() > a2.length(); }
	};

	int compare_proteins();
	bool protein_segment_aligned(int ri, int rj) const;
	void get_ca_positions(const Protein* prot,
		size_t res0, size_t len,
		std::vector<Position>& pcas) const;
	void align_proteins(const std::vector<Position>& ps1,
		std::vector<Position>& ps2, ProteinAlignment& a);
	int find_closest_residue(
		const Position& p, const std::vector<Position>& ps);
	void get_new_superposition(
		const std::vector<Position>& ps1,
		const std::vector<Position>& ps2, 
	    ProteinAlignment& a);
	void add_protein_alignment(
		const std::vector<Position>& ps1,
		const std::vector<Position>& ps2, 
		ProteinAlignment& a);
	void refine_protein_alignments(
		const std::vector<Position>& ps1,
		const std::vector<Position>& ps2);
	void merge_protein_alignments();
};

inline int
StructureComparison::n_atom_alignments() const
{ return m_atom_alignments.size(); }

inline int
StructureComparison::n_protein_alignments() const
{ return m_protein_alignments.size(); }

inline void
StructureComparison::ProteinAlignment::add_segment(int i, int j, int l)
{ segments.push_back(ProteinAlignmentSegment(i,j,l)); }

}
#endif
