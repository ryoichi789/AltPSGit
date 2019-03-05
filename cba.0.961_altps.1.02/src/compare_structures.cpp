#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include "compare_structures.h"
#include "structure_comparison.h"
#include "molecule.h"
#include "sdf_entry.h"
#include "sdf_reader.h"
#include "protein.h"
#include "pdb_entry.h"
#include "pdb_reader.h"
#include "bio_sequence.h"

using namespace Cba;

CompareStructures the_app;


// CompareStructures:
CompareStructures::CompareStructures() {}

// set_param_default:
void CompareStructures::set_param_default()
{
	m_prog = "compare_structures";
	m_infile1 = 0;
	m_infile2 = 0;
	m_chain1 = ' ';
	m_chain2 = ' ';
	m_outfile = 0;
	m_pdb_vs_pdb = false;
}

// set_param:
void CompareStructures::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!(m_infile1 && m_infile2)) {
		usage();
		std::exit(2);
	}
}

// run:
int CompareStructures::run()
{
	std::ostream* os;
	if (m_outfile) os = new std::ofstream(m_outfile);
	else os = &std::cout;

	if (m_pdb_vs_pdb) compare_proteins(*os);
	else compare_molecules(*os);

	if (os != &std::cout) delete os;
	return 0;
}

// usage:
void CompareStructures::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\t1: molecule structure 1 in SDF or PDB format [string] (unique; mandatory)" << endl;
	cerr << "\t2: molecule structure 2 in SDF or PDB format [string] (unique; mandatory)" << endl;
	cerr << "\tc1: chain id in protein 1 (valid only for PDB-vs-PDB comparison) [char] (unique; ' ' by default)" << endl;
	cerr << "\tc2: chain id in protein 2 (valid only for PDB-vs-PDB comparison) [char] (unique; ' ' by default)" << endl;
	cerr << "\tm: 0 if SDF vs SDF; 1 if PDB vs PDB [integer] (unique; 0 by default)" << endl;
	cerr << "\to: file for output [string] (unique; stdout by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void CompareStructures::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "1") == 0) {
		m_obj->m_infile1 = param.value;
	} else if (strcmp(param.key, "2") == 0) {
		m_obj->m_infile2 = param.value;
	} else if (strcmp(param.key, "c1") == 0) {
		m_obj->m_chain1 = param.value[0];
	} else if (strcmp(param.key, "c2") == 0) {
		m_obj->m_chain2 = param.value[0];
	} else if (strcmp(param.key, "m") == 0) {
		if (param.value[0] == '1') m_obj->m_pdb_vs_pdb = true;
	} else if (strcmp(param.key, "o") == 0) {
		m_obj->m_outfile = param.value;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

// compare_proteins:
void CompareStructures::compare_proteins(
	std::ostream& os) const
{
	StructureComparison cmp;
	Protein p1, p2;
	std::ifstream ifs1(m_infile1);
	std::ifstream ifs2(m_infile2);
	if (get_protein(ifs1, m_chain1, p1) &&
		get_protein(ifs2, m_chain2, p2)) {
		cmp.compare(&p1, &p2);
		print_protein_alignments(p1, p2, cmp, os);
	}
}

// get_protein:
bool CompareStructures::get_protein(
	std::istream& is, char chain, Protein& p) const
{
	PdbEntry ent;
	if (PdbReader::read(ent, is)) {
		for (int i = 0; i < ent.nmolecules(); i++) {
			if (ent.type_of_molecule(i) == PdbEntry::PROTEIN &&
				ent.nmr_model(i) <= 1 && ent.chainid(i) == chain) {
				ent.create_protein(i, p);
				return true;
			}
		}
	}
	return false;
}

// print_protein_alignment:
void CompareStructures::print_protein_alignments(
	const Protein& p1, const Protein& p2, 
	const StructureComparison& cmp,
	std::ostream& os) const
{
	using std::endl;
	os << '>' << p1.id() << ' ' << p1.nresidues() << endl;
	os << '>' << p2.id() << ' ' << p2.nresidues() << endl;
	os << "* " << cmp.n_protein_alignments()
		<< " alignments" << endl;

	BioSequence seq1, seq2;
	p1.get_seq(seq1);
	p2.get_seq(seq2);

	StructureComparison::ProteinAlignment a;
	for (int i = 0; i < cmp.n_protein_alignments(); i++) {
		os << "# " << i+1 << std::endl;
		cmp.get_protein_alignment(i, a);
		a.print(os);
		print_protein_alignment_in_seq(seq1, seq2, a, os);
		os << '/' << endl;
	}
	for (int i = 0; i < cmp.n_protein_alignments(); i++) {
		cmp.get_protein_alignment(i, a);
	}
}

// print_protein_alignment_in_seq:
void CompareStructures::print_protein_alignment_in_seq(
	const BioSequence& seq1, const BioSequence& seq2,
	const StructureComparison::ProteinAlignment& a,
	std::ostream& os) const
{
	using std::endl;

	for (int i = 0; i < a.segments.size(); i++) {
		const StructureComparison::
			ProteinAlignmentSegment& s
				= a.segments[i];

		// segment of sequence 1
		os << endl;
		os << std::fixed << std::setw(5) << std::left
			<< s.r1+1 << ' ';
		for (int j = s.r1; j < s.r1+s.length; j++)
			os << seq1.at(j);
		os << ' ' << s.r1+s.length << endl;;

		// segment of sequence 2
		os << std::fixed << std::setw(5) << std::left
			<< s.r2+1 << ' ';
		for (int j = s.r2; j < s.r2+s.length; j++)
			os << seq2.at(j);
		os << ' ' << s.r2+s.length << endl;;

		// residue match
		os << "      ";
		for (int j = 0; j < s.length; j++)
			if (seq1.at(s.r1+j) == seq2.at(s.r2+j))
				os << '*';
			else os << ' ';
		os << endl;
	}
}

// compare_molecules:
void CompareStructures::compare_molecules(
	std::ostream& os) const
{
	StructureComparison cmp;
	SdfEntry e1, e2;
	Molecule m1, m2;

	std::ifstream ifs1(m_infile1);
	SdfReader reader1(ifs1);
	while (reader1.read(e1)) {
		e1.create_molecule(m1);
		m1.assign_atom_types();

		std::ifstream ifs2(m_infile2);
		SdfReader reader2(ifs2);
		while (reader2.read(e2)) {
			e2.create_molecule(m2);
			m2.assign_atom_types();
			cmp.compare(&m1, &m2);
			print_atom_alignments(m1, m2, cmp, os);
		}
	}
}

// print_atom_alignment:
void CompareStructures::print_atom_alignments(
	const Molecule& m1, const Molecule& m2,
	const StructureComparison& cmp,
	std::ostream& os) const
{
	using std::endl;
	if (cmp.n_atom_alignments() == 0) return;

	os << '>' << m1.id() << ' ' << m1.natoms() << endl;
	os << '>' << m2.id() << ' ' << m2.natoms() << endl;
	os << "* " << cmp.n_atom_alignments()
		<< " alignments" << endl;
	StructureComparison::AtomAlignment a;
	for (int i = 0; i < cmp.n_atom_alignments(); i++) {
		cmp.get_atom_alignment(i, a);
		os << "# " << i+1 << ' '
			<< a.pairs.size() << ' '
			<< std::fixed << std::setprecision(2)
			<< std::setw(5) << a.rmsd << ' '
			<< a.score << std::endl;
		for (int j = 0; j < a.pairs.size(); j++) {
			int a1 = a.pairs[j].a1;
			int a2 = a.pairs[j].a2;
			os << a1+1 << ' '
				<< m1.atom(a1)->name << ' '
				<< (int)m1.atom(a1)->pc_class << ' '
				<< " : "
				<< a2+1 << ' '
				<< m2.atom(a2)->name << ' '
				<< (int)m2.atom(a2)->pc_class << std::endl;
		}
	}
	os << '/' << std::endl;
}
