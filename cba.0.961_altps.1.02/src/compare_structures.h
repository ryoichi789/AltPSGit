#ifndef CBA_COMPARE_STRUCTURES_H
#define CBA_COMPARE_STRUCTURES_H

#include <iostream>
#include <functional>
#include "app.h"
#include "molecule.h"
#include "protein.h"
#include "bio_sequence.h"
#include "structure_comparison.h"

namespace Cba {

class CompareStructures : public App {
public:
	CompareStructures();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();

	const char* m_infile1;
	const char* m_infile2;
	const char* m_outfile;
	bool m_pdb_vs_pdb;
	char m_chain1; // chain id in a PDB entry
	char m_chain2;

	void compare_proteins(std::ostream& os) const;
	bool get_protein(std::istream& is,
		char chain, Protein& p) const;
	void print_protein_alignments(
		const Protein& p1, const Protein& p2,
		const StructureComparison& cmp,
		std::ostream& os) const;
	void print_protein_alignment_in_seq(
		const BioSequence& seq1, const BioSequence& seq2,
		const StructureComparison::ProteinAlignment& a,
		std::ostream& os) const;

	void compare_molecules(std::ostream& os) const;
	void print_atom_alignments(
		const Molecule& m1, const Molecule& m2,
		const StructureComparison& cmp,
		std::ostream& os) const;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(CompareStructures* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		CompareStructures* m_obj;
	};
	friend class SetParam;
};

}
#endif
