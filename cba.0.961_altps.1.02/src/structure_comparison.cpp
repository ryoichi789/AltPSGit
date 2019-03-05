#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "structure_comparison.h"
#include "superposer.h"

using namespace Cba;
using std::vector;

const int StructureComparison::MIN_ATOM_PAIRS = 10;
const double StructureComparison::MAX_DIST_ATOM_PAIR = 1.0;
const int StructureComparison::MIN_ATOM_ALIGNMENT_SCORE = 4;

const int StructureComparison::LENGTH_OF_SEGMENT = 15;
const int StructureComparison::MIN_LEN_ALIGNMENT_SEGMENT = 10;
const int StructureComparison::MIN_LEN_ALIGNMENT = 40;
const double StructureComparison::MAX_INIT_RMSD_PROTEIN = 1.0;
const double StructureComparison::MAX_DIST_RESIDUE_PAIR = 2.5;

//
StructureComparison::StructureComparison() {}

//
StructureComparison::~StructureComparison() {}

// compare:
int StructureComparison::compare(Molecule* mol1, Molecule* mol2)
{
	m_mol1 = mol1;
	m_mol2 = mol2;
	m_prot1 = m_prot2 = 0;
	m_prot1 = dynamic_cast<Protein*>(m_mol1);
	m_prot2 = dynamic_cast<Protein*>(m_mol2);
	if (m_prot1 && m_prot2)
		return compare_proteins();
	return compare();
}

// get_atom_alignment:
void StructureComparison::get_atom_alignment(
	size_t i, AtomAlignment& a) const
{
	if (i < m_atom_alignments.size())
		a = m_atom_alignments[i];
}

// get_protein_alignment:
void StructureComparison::get_protein_alignment(
	size_t i, ProteinAlignment& a) const
{
	if (i < m_protein_alignments.size())
		a = m_protein_alignments[i];
}

// compare:
int StructureComparison::compare()
{
	m_atom_alignments.clear();

	vector<Position> ps1;
	get_atom_positions(m_mol1, ps1);
	vector<Position> ps2;

	vector<AtomTriad> triads1, triads2;
	get_atom_triads(triads1, triads2);
	AtomAlignment a;
	for (int i = 0; i < triads1.size(); i++) {
		for (int j = 0; j < triads2.size(); j++) {
			if (atom_triads_already_aligned(
				triads1[i], triads2[j])) continue;
			if (match_atom_triads(
				triads1[i], triads2[j])) {
				superpose_atom_triads(
					triads1[i], triads2[j], a);
				get_atom_positions(m_mol2, ps2);
				align_atoms(ps1, ps2, a);

				get_atom_positions(m_mol2, ps2);
				add_atom_alignment(ps1, ps2, a);
			}
		}
	}

	return m_atom_alignments.size();
}

// get_atom_positions:
void StructureComparison::get_atom_positions(
	const Molecule* mol, vector<Position>& ps) const
{
	ps.clear();
	for (int i = 0; i < mol->natoms(); i++)
		ps.push_back(mol->atom(i)->pos);
}

// get_atom_triads:
void StructureComparison::get_atom_triads(
	vector<AtomTriad>& triads1, vector<AtomTriad>& triads2)
{
	get_atom_triads(*m_mol1, triads1);
	get_atom_triads(*m_mol2, triads2);

	vector<AtomTriad> triads3;
	for (int i = 0; i < triads2.size(); i++)
		triads3.push_back(AtomTriad(triads2[i].a0,
			triads2[i].a2, triads2[i].a1));
	triads2.insert(triads2.end(),
		triads3.begin(), triads3.end());
}

// get_atom_triads:
void StructureComparison::get_atom_triads(
	Molecule& mol, vector<AtomTriad>& triads)
{
	triads.clear();
	Atom* a0;
	vector<Atom*> as; // atoms bonded to an atom
	for (int i = 0; i < mol.natoms(); i++) {
		a0 = mol.atom(i);
		int c0 = a0->pc_class;
		mol.get_atoms_bonded_to(a0, as);
		for (int j = 0; j < (int)as.size()-1; j++) {
			int c1 = as[j]->pc_class;
			for (int k = j+1; k < as.size(); k++) {
				int c2 = as[k]->pc_class;
				if ((c0 != Atom::HYDROPHOBIC &&
						c0 != Atom::NONE) ||
					(c1 != Atom::HYDROPHOBIC &&
						c1 != Atom::NONE) ||
					(c2 != Atom::HYDROPHOBIC &&
						c2 != Atom::NONE))
					triads.push_back(
						AtomTriad(a0,as[j],as[k]));
			}
		}
	}
}

// atom_triads_already_aligned:
bool StructureComparison::atom_triads_already_aligned
	(const AtomTriad& t1, const AtomTriad& t2) const
{
	for (int i = 0; i < m_atom_alignments.size(); i++) {
		const AtomAlignment& a = m_atom_alignments[i];
		int k = 0;
		for (int j = 0; j < a.pairs.size(); j++) {
			const Atom* a1 = m_mol1->atom(a.pairs[j].a1);
			const Atom* a2 = m_mol2->atom(a.pairs[j].a2);
			if ((a1 == t1.a0 && a2 == t2.a0) ||
				(a1 == t1.a1 && a2 == t2.a1) ||
				(a1 == t1.a2 && a2 == t2.a2))
				k++;
		}
		if (k == 3) return true;
	}
	return false;
}

// match_atom_triads:
bool StructureComparison::match_atom_triads(
	const AtomTriad& t1, const AtomTriad& t2) const
{
	if (match_atoms(t1.a0, t2.a0) &&
		match_atoms(t1.a1, t2.a1) &&
		match_atoms(t1.a2, t2.a2))
		return true;
	return false;
}

// match_atoms:
bool StructureComparison::match_atoms(
	const Atom* a1, const Atom* a2) const
{
	if (a1->pc_class == Atom::CATION) {
		if (a2->pc_class == Atom::CATION)
			return true;
	} else if (a1->pc_class == Atom::ANION) {
		if (a2->pc_class == Atom::ANION)
			return true;
	} else if (a1->pc_class == Atom::DONOR) {
		if (a2->pc_class == Atom::DONOR ||
			a2->pc_class == Atom::POLAR)
			return true;
	} else if (a1->pc_class == Atom::ACCEPTOR) {
		if (a2->pc_class == Atom::ACCEPTOR ||
			a2->pc_class == Atom::POLAR)
			return true;
	} else if (a1->pc_class == Atom::POLAR) {
		if (a2->pc_class == Atom::DONOR ||
			a2->pc_class == Atom::ACCEPTOR ||
			a2->pc_class == Atom::POLAR)
			return true;
	} else if (a1->pc_class == Atom::HYDROPHOBIC) {
		if (a2->pc_class == Atom::HYDROPHOBIC)
			return true;
	} else if (a1->pc_class == Atom::NONE) {
		if (a2->pc_class == Atom::NONE)
			return true;
	}
	return false;
}

// superpose_atom_triads:
void StructureComparison::superpose_atom_triads(
	const AtomTriad& t1, const AtomTriad& t2,
	AtomAlignment& a)
{
	vector<Position> ps1;
	ps1.push_back(t1.a0->pos);
	ps1.push_back(t1.a1->pos);
	ps1.push_back(t1.a2->pos);
	vector<Position> ps2;
	ps2.push_back(t2.a0->pos);
	ps2.push_back(t2.a1->pos);
	ps2.push_back(t2.a2->pos);
	a.spp.superpose(ps2, ps1);
}

// align_atoms:
void StructureComparison::align_atoms(
	vector<Position>& ps1, vector<Position>& ps2,
	AtomAlignment& a)
{
	a.pairs.clear();
	int dn, n0;
	for (;;) {
		n0 = a.pairs.size();
		a.pairs.clear();
		for (int i = 0; i < ps2.size(); i++)
			a.spp.transform(ps2[i]);
		for (int i = 0; i < ps1.size(); i++) {
			int j = get_equiv_atom(ps1[i], ps2);
			if (j < ps2.size())
				a.pairs.push_back(AtomPair(i,j));
		}
		dn = a.pairs.size() - n0;
		if (dn > 0)
			get_new_superposition(ps1, ps2, a);
		else
			break;
	}
}

// get_equiv_atom:
int StructureComparison::get_equiv_atom(
	const Position& p, const vector<Position>& ps) const
{
	int imin = ps.size();
	double dmin = MAX_DIST_ATOM_PAIR + 1.0;
	double d;
	for (int i = 0; i < ps.size(); i++) {
		if (p.is_within_range(ps[i], MAX_DIST_ATOM_PAIR)
			&& (d = p.distance_from(ps[i])) < dmin) {
			dmin = d;
			imin = i;
		}
	}
	return imin;
}

// get_new_superposition:
void StructureComparison::get_new_superposition(
	const vector<Position>& ps1,
	const vector<Position>& ps2,
	AtomAlignment& a)
{
	vector<Position> qs1;
	vector<Position> qs2;
	for (int i = 0; i < a.pairs.size(); i++) {
		qs1.push_back(ps1[a.pairs[i].a1]);
		qs2.push_back(ps2[a.pairs[i].a2]);
	}
	a.rmsd = a.spp.superpose(qs2, qs1);
}

// add_atom_alignment:
void StructureComparison::add_atom_alignment(
	const vector<Position>& ps1,
	const vector<Position>& ps2,
	AtomAlignment& a)
{
	if (a.pairs.size() < MIN_ATOM_PAIRS) return;
	get_atom_alignment_score(a);
	if (a.score < MIN_ATOM_ALIGNMENT_SCORE) return;

	get_new_superposition(ps1, ps2, a);
	m_atom_alignments.push_back(a);
}

// get_atom_alignment_score:
void StructureComparison::get_atom_alignment_score(
	AtomAlignment& a) const
{
	a.score = 0;
	int c1, c2; // physico-chemical classes of atoms
	for (int i = 0; i < a.pairs.size(); i++) {
		c1 = m_mol1->atom(a.pairs[i].a1)->pc_class;
		c2 = m_mol2->atom(a.pairs[i].a2)->pc_class;
		if (c1 == Atom::CATION) {
			if (c2 == Atom::CATION) a.score++;
			else a.score--;
		} else if (c1 == Atom::ANION) {
			if (c2 == Atom::ANION) a.score++;
			else a.score--;
		} else if (c1 == Atom::DONOR) {
			if (c2 == Atom::DONOR ||
				c2 == Atom::POLAR) a.score++;
			else a.score--;
		} else if (c1 == Atom::ACCEPTOR) {
			if (c2 == Atom::ACCEPTOR ||
				c2 == Atom::POLAR) a.score++;
			else a.score--;
		} else if (c1 == Atom::HYDROPHOBIC) {
			if (c2 != Atom::HYDROPHOBIC) a.score--;
		} else if (c1 == Atom::NONE) {
			if (c2 != Atom::NONE) a.score--;
		}
	}
}

// compare_proteins:
int StructureComparison::compare_proteins()
{
	m_protein_alignments.clear();
	vector<Position> pcas1, pcas2; // positions of Ca's
	get_ca_positions(m_prot1, 0, m_prot1->nresidues(), pcas1);
	get_ca_positions(m_prot2, 0, m_prot2->nresidues(), pcas2);
	// distance between ends of each segment
	double* ds1 = new double[m_prot1->nresidues()];
	for (int i = 0; i <= m_prot1->nresidues()-LENGTH_OF_SEGMENT; i++)
		ds1[i] = pcas1[i].distance_from(pcas1[i+LENGTH_OF_SEGMENT-1]);
	double* ds2 = new double[m_prot2->nresidues()];
	for (int i = 0; i <= m_prot2->nresidues()-LENGTH_OF_SEGMENT; i++)
		ds2[i] = pcas2[i].distance_from(pcas2[i+LENGTH_OF_SEGMENT-1]);

	ProteinAlignment a;
	vector<Position> pseg1, pseg2; // positions of Ca's
	vector<Position> pcas2_work;
	for (int i = 0; i <= m_prot1->nresidues()-LENGTH_OF_SEGMENT; i++) {
		get_ca_positions(m_prot1, i, LENGTH_OF_SEGMENT, pseg1);
		for (int j = 0; j <= m_prot2->nresidues()-LENGTH_OF_SEGMENT; j++) {
			if (std::fabs(ds1[i]-ds2[j]) > MAX_INIT_RMSD_PROTEIN) continue;
			if (protein_segment_aligned(i, j)) continue;
			get_ca_positions(m_prot2, j, LENGTH_OF_SEGMENT, pseg2);
			if (a.spp.superpose(pseg2, pseg1) > MAX_INIT_RMSD_PROTEIN) continue;
			get_ca_positions(m_prot2, 0, m_prot2->nresidues(), pcas2_work);
			align_proteins(pcas1, pcas2_work, a);
			add_protein_alignment(pcas1, pcas2, a);
		}
	}
	delete [] ds1;
	delete [] ds2;

	refine_protein_alignments(pcas1, pcas2);
	return m_protein_alignments.size();
}

// refine_protein_alignments:
void StructureComparison::refine_protein_alignments(
	const std::vector<Position>& ps1, 
	const std::vector<Position>& ps2)
{
	merge_protein_alignments();
	for (vector<ProteinAlignment>::iterator
		i = m_protein_alignments.begin();
		i != m_protein_alignments.end(); ++i)
		get_new_superposition(ps1, ps2, *i);
	std::sort(m_protein_alignments.begin(),
		m_protein_alignments.end(),
		CmpProteinAlignments());
}

// residue_pair_aligned:
bool StructureComparison::protein_segment_aligned(int ri, int rj) const
{
	for (int i = 0; i < m_protein_alignments.size(); i++) {
		const ProteinAlignment& a = m_protein_alignments[i];
		for (int j = 0; j < a.segments.size(); j++) {
			const ProteinAlignmentSegment& seg = a.segments[j];
			if (seg.r1 <= ri &&
				ri+LENGTH_OF_SEGMENT <= seg.r1+seg.length &&
				ri-rj == seg.r1-seg.r2)
				return true;
		}
	}
	return false;
}

// get_ca_positions:
void StructureComparison::get_ca_positions(
	const Protein* prot, size_t res0, size_t len, vector<Position>& pcas) const
{
	Position p0;
	pcas.clear();
	for (int i = res0; i < res0+len; i++) {
		const Protein::Residue& r = prot->residue(i);
		if (r.ca) pcas.push_back(r.ca->pos);
		else pcas.push_back(p0);
	}
}

// align_proteins:
void StructureComparison::align_proteins(
	const vector<Position>& ps1, vector<Position>& ps2, ProteinAlignment& a)
{
	a.segments.clear();
	int l0, dl;
	do {
		for (int i = 0; i < ps2.size(); i++) a.spp.transform(ps2[i]);
		l0 = a.length();
		a.segments.clear();
		for (int i0 = 0; i0 <= ps1.size()-MIN_LEN_ALIGNMENT_SEGMENT; ) {
			int j0 = find_closest_residue(ps1[i0], ps2);
			if (j0 < ps2.size()) {
				int i1, j1, l;
				for (i1=i0+1, j1=j0+1; i1<ps1.size() && j1<ps2.size() &&
					ps1[i1].is_within_range(ps2[j1], MAX_DIST_RESIDUE_PAIR);
					i1++, j1++) ;
				i1--; j1--;
				if ((l = i1-i0+1) >= MIN_LEN_ALIGNMENT_SEGMENT)
					a.add_segment(i0, j0, l);
				i0 += l;
			} else {
				i0++;
			}
		}
		get_new_superposition(ps1, ps2, a);
		dl = a.length() - l0;
	} while (a.length() >= MIN_LEN_ALIGNMENT && dl > 0);
}

// find_closest_residue:
int StructureComparison::find_closest_residue(
	const Position& p, const vector<Position>& ps)
{
	int imin = ps.size();
	double dmin = MAX_DIST_RESIDUE_PAIR*MAX_DIST_RESIDUE_PAIR+1.0;
	double d;
	for (int i = 0; i < ps.size(); i++) {
		if (p.is_within_range(ps[i], MAX_DIST_RESIDUE_PAIR) &&
			(d = p.square_distance_from(ps[i])) < dmin) {
			dmin = d;
			imin = i;
		}
	}
	return imin;
}

// get_new_superposition:
void StructureComparison::get_new_superposition(
	const vector<Position>& ps1, const vector<Position>& ps2,
	ProteinAlignment& a)
{
	vector<Position> qs1;
	vector<Position> qs2;
	for (int i = 0; i < a.segments.size(); i++) {
		ProteinAlignmentSegment& seg = a.segments[i];
		for (int k = seg.r1; k < seg.r1+seg.length; k++)
			qs1.push_back(ps1[k]);
		for (int k = seg.r2; k < seg.r2+seg.length; k++)
			qs2.push_back(ps2[k]);
	}
	a.rmsd = a.spp.superpose(qs2, qs1); // superpose protein 2 onto 1
}

// add_protein_alignment:
void StructureComparison::add_protein_alignment(
	const vector<Position>& ps1, const vector<Position>& ps2,
	ProteinAlignment& a)
{
	if (a.length() < MIN_LEN_ALIGNMENT) return;
	for (vector<ProteinAlignment>::iterator
		i = m_protein_alignments.begin(); i != m_protein_alignments.end(); ++i) {
		if (a == *i) return;
		if (i->contains(a)) return;
		if (a.contains(*i)) {
			m_protein_alignments.erase(i);
			break;
		}
	}
	get_new_superposition(ps1, ps2, a);
	m_protein_alignments.push_back(a);
}

// merge_protein_alignments:
void StructureComparison::merge_protein_alignments()
{
	bool merged;
	do {
		merged =false;
		for (vector<ProteinAlignment>::iterator
			i = m_protein_alignments.begin();
			i != m_protein_alignments.end(); ++i) {
			vector<ProteinAlignment>::iterator j = i;
			for (++j; j != m_protein_alignments.end(); ++j) {
				if (i->merges(*j)) {
					m_protein_alignments.erase(j);
					merged = true;
					break;
				}
			}
			if (merged) break;
		}
	} while (merged);
}

// ProteinAlignmentSegment::contains:
bool StructureComparison::ProteinAlignmentSegment::contains(
	const ProteinAlignmentSegment& s) const
{
	if (r1-r2 != s.r1-s.r2) return false;
	if (r1 <= s.r1 && s.r1+s.length <= r1+length) return true;
	return false;
}

// ProteinAlignmentSegment::overlaps:
bool StructureComparison::ProteinAlignmentSegment::overlaps(
	const ProteinAlignmentSegment& s) const
{
	if (r1-r2 != s.r1-s.r2) return false;
	if (r1+length <= s.r1 || s.r1+s.length <= r1) return false;
	return true;
}

// ProteinAlignmentSegment::merges:
bool StructureComparison::ProteinAlignmentSegment::merges(
	const ProteinAlignmentSegment& s)
{
	if (!overlaps(s)) return false;

	// starts and ends of this alignment segment
	int r1s = r1;
	int r1e = r1 + length -1;
	// starts and ends of alignment segment s
	int q1s = s.r1;
	int q1e = s.r1 + s.length -1;

	int dr1 = r1s - q1s;
	if (dr1 > 0) {
		r1 = q1s;
		r2 -= dr1;
	}
	dr1 = r1e - q1e;
	if (dr1 < 0) length = q1e - r1 + 1;

	return true;
}

// ProteinAlignment::length:
int StructureComparison::ProteinAlignment::length() const
{
	int l = 0;
	for (int i = 0; i < segments.size(); i++)
		l += segments[i].length;
	return l;
}

// ProteinAlignment::print:
void StructureComparison::ProteinAlignment::print(std::ostream& os) const
{
	using std::endl;
	os << "length: " << length() << endl;
	os << "rmsd: "
		<< std::fixed << std::setprecision(2) << std::setw(5)
		<< rmsd << endl;
	os << "nsegments: " << segments.size() << endl;
	for (int i = 0; i < segments.size(); i++)
		os << segments[i].r1+1 << ' ' << segments[i].r1+segments[i].length
			<< " : "
			<< segments[i].r2+1 << ' ' << segments[i].r2+segments[i].length
			<< " = "
			<< segments[i].length << endl;
	os << '/' << endl;
}

// ProteinAlignment::operator==:
bool StructureComparison::ProteinAlignment::operator==(
	const ProteinAlignment& a) const
{
	int n = segments.size();
	if (n != a.segments.size()) return false;
	for (int i = 0; i < n; i++)
		if (segments[i] != a.segments[i]) return false;
	return true;
}

// ProteinAlignment::contains:
bool StructureComparison::ProteinAlignment::contains(
	const ProteinAlignment& a) const
{
	int k = 0;
	for (int i = 0; i < a.segments.size(); i++) {
		for ( ; k < segments.size(); k++)
			if (segments[k].contains(a.segments[i])) break;
		if (k == segments.size()) return false;
		k++;
	}
	return true;
}

// ProteinAlignment::overlaps:
bool StructureComparison::ProteinAlignment::overlaps(
	const ProteinAlignment& a) const
{
	for (int i = 0; i < segments.size(); i++)
		for (int j = 0; j < a.segments.size(); j++)
			if (segments[i].overlaps(a.segments[j]))
				return true;
	return false;
}

// ProteinAlignment::merges:
bool StructureComparison::ProteinAlignment::merges(
	const ProteinAlignment& a)
{
	if (!overlaps(a)) return false;

	segments.insert(segments.end(),
		a.segments.begin(), a.segments.end());
	std::sort(segments.begin(), segments.end(),
		CmpProteinAlignmentSegment());
	bool merged;
	do {
		merged = false;
		for (vector<ProteinAlignmentSegment>::iterator
			i = segments.begin(); i != segments.end(); ++i) {
			vector<ProteinAlignmentSegment>::iterator j = i;
			if (++j == segments.end()) break;
			if (i->merges(*j)) {
				segments.erase(j);
				merged = true;
				break;
			};
		}
	} while(merged);

	return true;
}
