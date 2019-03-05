#include <iostream>
#include <list>
#include <algorithm>
#include "bio_seq_cmp_identical.h"
#include "bio_seq_equiv_region.h"
#include "bio_sequence.h"
using namespace Cba;

// BioSeqCmpIdentical:
BioSeqCmpIdentical::BioSeqCmpIdentical(size_t min_len_region)
	: m_min_len_region(min_len_region) {}

// ~BioSeqCmp:
BioSeqCmpIdentical::~BioSeqCmpIdentical() {}

// compare:
void BioSeqCmpIdentical::compare(
	const BioSequence* seq1, const BioSequence* seq2, std::ostream& os)
{
	m_seq1 = seq1;
	m_seq2 = seq2;
	int l1 = m_seq1->length();
	int l2 = m_seq2->length();
	int j1 = l1 - m_min_len_region;
	int j2 = l2 - m_min_len_region;
	m_equiv_regions.clear();
	for (int i = 0; i <= j1; i++) compare_from(i, 0, j1, j2, l1, l2);
	for (int i = 1; i <= j2; i++) compare_from(0, i, j1, j2, l1, l2);
	std::sort(m_equiv_regions.begin(), m_equiv_regions.end(), BioSeqEquivRegionCmp());

	print(os);
}

// print:
void BioSeqCmpIdentical::print(std::ostream& os)
{
	using std::endl;

	os << ">" << m_seq1->id() << " " << m_seq1->length() << endl;
	os << ">" << m_seq2->id() << " " << m_seq2->length() << endl;
	for (int i = 0; i < m_equiv_regions.size(); i++) {
		os << (m_equiv_regions[i].pos1+1) << " "
			<< (m_equiv_regions[i].pos1+m_equiv_regions[i].len) << " "
			<< (m_equiv_regions[i].pos2+1) << " "
			<< (m_equiv_regions[i].pos2+m_equiv_regions[i].len) << " "
			<< m_equiv_regions[i].len << endl;
	}
}

// compare_from:
void BioSeqCmpIdentical::compare_from(int s1, int s2, int e1, int e2, int l1, int l2)
{
	int p1, p2, l;
	while (s1<=e1 && s2<=e2) {
		for (p1=s1, p2=s2; p1<l1 && p2<l2; p1++, p2++)
			if (m_seq1->at(p1) != m_seq2->at(p2)) break;
		// l = length of identical subsequence
		if ((l=p1-s1) >= m_min_len_region)
			m_equiv_regions.push_back(BioSeqEquivRegion(s1,s2,l));
		s1 = p1+1;
		s2 = p2+1;
	}
}
