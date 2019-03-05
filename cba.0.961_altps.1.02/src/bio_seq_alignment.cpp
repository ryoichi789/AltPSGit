#include <algorithm>
#include <ios>
#include <iomanip>
#include "bio_seq_alignment.h"
//#include "bio_seq_equiv_region.h"

using namespace Cba;

// BioSeqAlignment:
BioSeqAlignment:: BioSeqAlignment(
	const BioSequence* seq1, const BioSequence* seq2)
	: m_seq1(seq1), m_seq2(seq2) {}

// nmatches:
int BioSeqAlignment::nmatches() const
{
	int n = 0;
	for (int i = 0; i < length(); i++)
		if (m_alignment.at(i).first != GAP &&
			m_alignment.at(i).second != GAP) n++;
	return n;
}

// get_seq_identity:
double BioSeqAlignment::get_seq_identity() const
{
	int n = 0; // number of identical residues
	int pos1, pos2;
	for (int i = 0; i < length(); i++) {
		if ((pos1 = m_alignment.at(i).first) == GAP) continue;
		if ((pos2 = m_alignment.at(i).second) == GAP) continue;
		if (m_seq1->at(pos1) == m_seq2->at(pos2)) n++;
	}

	if (length() == 0) return 0.0;
	return 100.0*(double)n/length();
}

// print:
void BioSeqAlignment::print(std::ostream& os) const
{
	using std::endl;

	int l = length();
	os << m_seq1->id() << " " << m_seq1->length() << endl;
	os << m_seq2->id() << " " << m_seq2->length() << endl;
	os << "length: " << l << "  ";
	os << "ngaps: " << ngaps() << "  ";
	os << "score: " << get_score() << "  ";
	os << "id: " << std::fixed << std::setprecision(1) <<
		get_seq_identity() << endl;

	for (int i = 0; i < length(); i += LINE_SEQ_FOLD)
		print_range(os, i, ((i+LINE_SEQ_FOLD) > l ? (l-i) : LINE_SEQ_FOLD));

	os << '/' << endl; // line beginning with '/' indicates end of alignment
}

// print_range:
void BioSeqAlignment::print_range(std::ostream& os, size_t pos, size_t len) const
{
	print_scale(os, len);
	print_seq(0, os, pos, len);
	print_seq(1, os, pos, len);
	print_star(os, pos, len);
}

// print_scale:
void BioSeqAlignment::print_scale(std::ostream& os, size_t len) const
{
	os << "       "; // 7 spaces at line head
	for (int i = 1; i <= len; i++)
		if (i%10) os << ' ';
		else os << '+'; // scale ('+' for every 10 characters)
	os << std::endl;
}

// print_seq:
void BioSeqAlignment::print_seq(int seq_i, std::ostream& os, size_t pos, size_t len) const
{
	const BioSequence* seq = (seq_i == 0 ? m_seq1 : m_seq2);
	int pa; // position in alignment
	int ps; // position in sequence

	// line head
	for (pa = pos; pa < length() && residue_number_at(pa, seq_i) == GAP; pa++)
		;
	if (pa == length())
		os << std::setw(6) << std::right << seq->length() << ' ';
	else
		os << std::setw(6) << std::right
			<< residue_number_at(pa, seq_i)+1 << ' ';

	// alignment
	for (int pa = pos; pa < (pos+len); pa++)
		if ((ps = residue_number_at(pa, seq_i)) == GAP)
			os << '-';
		else
			os << seq->at(ps);

	// line tail
	for (pa = (pos+len-1); pa >= 0 && residue_number_at(pa, seq_i) == GAP; pa--)
		;
	if (pa < 0)
		os << " 1";
	else
		os << ' ' << std::setw(6) << std::left << residue_number_at(pa, seq_i)+1;

	os << std::endl;
}

// residue_number_at:
int BioSeqAlignment::residue_number_at(size_t pos, int seq_i) const
{
	if (seq_i == 0) return m_alignment.at(pos).first;
	return m_alignment.at(pos).second;
}

// PrintAlignment::print_star:
void BioSeqAlignment::print_star(std::ostream& os, size_t pos, size_t len) const
{
	os << "       "; // 7 spaces at line head
	int ps1, ps2; // positions in sequences
	for (int i = pos; i < (pos+len); i++)
		if ((ps1 = m_alignment.at(i).first) != GAP &&
			(ps2 = m_alignment.at(i).second) != GAP &&
			m_seq1->at(ps1) == m_seq2->at(ps2))
			os << '*'; // identical residues
		else
			os << ' ';
	os << std::endl;
}
