#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "dna_sequence.h"
#include "bio_sequence.h"
#include "genetic_code.h"
using namespace Cba;

// ~DnaSequence():
DnaSequence::~DnaSequence() {}

// clear:
void DnaSequence::clear()
{
	BioSequence::clear();
	clear_data();
}

// complement_at:
char DnaSequence::complement_at(size_t i)
{
	if (i >= length()) return 0;
	char c = at(i);
	if (c == 'A' || c == 'a') return 'T';
	if (c == 'C' || c == 'c') return 'G';
	if (c == 'G' || c == 'g') return 'C';
	if (c == 'T' || c == 't') return 'A';
	if (c == 'U' || c == 'u') return 'A';
	if (c == 'R' || c == 'r') return 'Y';
	if (c == 'Y' || c == 'y') return 'R';
	return 'N';
}

// make_c_seq:
void DnaSequence::make_c_seq()
{
	int l = length();
	if (m_c_seq.length() == l) return; // already made
	clear_data();
	for (int i = l-1; i >= 0; i--)
		m_c_seq.push_back(complement_at(i));
}

// find_orfs:
int DnaSequence::find_orfs()
{
	make_c_seq();
	find_stop_codons();
	for (int i = 1; i <= 3; i++) {
		find_orfs_in_reading_frame(i);
		find_orfs_in_reading_frame(-i);
	}
	sort_orfs();
	return m_orfs.size();
}

// find_stop_codons:
void DnaSequence::find_stop_codons()
{
	int k = length() - 3; // position of last codon
	for (int j = 0; j < 3; j++) { // for each reading frame
		for (int i = j; i <= k; i += 3) { // for each codon
			if (GeneticCode::is_stop_codon(at(i), at(i+1), at(i+2)))
				m_stop_codons[j].push_back(i);
			if (GeneticCode::is_stop_codon(m_c_seq[i], m_c_seq[i+1], m_c_seq[i+2]))
				m_stop_codons[j+3].push_back(i);
		}
	}
}

// find_orf_in_reading_frame:
void DnaSequence::find_orfs_in_reading_frame(int frame_i)
{
	Orf orf;
	orf.m_reading_frame = frame_i;
	// frame +1->0, +2->1, +3->2, -1->3, -2->4, -3->5
	if (frame_i > 0) --frame_i;
	else frame_i = -frame_i + 2;

	orf.m_start = frame_i % 3; // start of current ORF
	int p;
	for (int i = 0; i < m_stop_codons[frame_i].size(); i++) {
		p = m_stop_codons[frame_i].at(i);
		if (p > orf.m_start) {
			orf.m_end = p - 1;
			m_orfs.push_back(orf);
		}
		orf.m_start = p + 3; // start of next ORF
	}

	// ORF from last stop codon to end of sequence
	if ((p = length()-orf.m_start) >= 3) {
		orf.m_end = (length()-1) - (p%3);
		m_orfs.push_back(orf);
	}
}

// sort_orfs:
void DnaSequence::sort_orfs()
{
	std::sort(m_orfs.begin(), m_orfs.end(), OrfCmp());
}

// clear_data:
void DnaSequence::clear_data()
{
	m_orfs.clear();
	m_c_seq.erase();
	for (int i = 0; i < 6; i++)
		m_stop_codons[i].clear();
}

// print_orfs:
void DnaSequence::print_orfs(std::ostream& os, int min_len) const
{
	std::for_each(m_orfs.begin(), m_orfs.end(), OrfPrint(os, min_len));
}

// OrfPrint::operator():
void DnaSequence::OrfPrint::operator()(const Orf& orf)
{
	if (orf.length() < m_min_len) return;
	m_os << orf.m_reading_frame << '\t'
		<< (orf.m_start+1) << '\t' << (orf.m_end+1) << '\t'
		<< orf.length() << std::endl;
}

// translate_orf:
bool DnaSequence::translate_orf(size_t orf_i, BioSequence& seq, int min_len) const
{
	seq.clear();
	if (orf_i >= m_orfs.size()) return false;
	const Orf& orf = m_orfs.at(orf_i);
	if (orf.length() < min_len) return false;

	std::ostringstream oss;
	oss << id() << '_' << orf.m_reading_frame << '_'
		<< (orf.m_start+1) << '_' << (orf.m_end+1) << '_' << orf.length();
	seq.set_id(oss.str());

	if (orf.m_reading_frame > 0) {
		for (int i = orf.m_start; i < orf.m_end; i += 3)
			seq.add(GeneticCode::decode(at(i), at(i+1), at(i+2)));
	} else {
		for (int i = orf.m_start; i < orf.m_end; i += 3)
			seq.add(GeneticCode::decode(m_c_seq.at(i), m_c_seq.at(i+1), m_c_seq.at(i+2)));
	}

	return true;
}
