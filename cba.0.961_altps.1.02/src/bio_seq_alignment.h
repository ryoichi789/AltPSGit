#ifndef CBA_BIO_SEQ_ALIGNMENT_H
#define CBA_BIO_SEQ_ALIGNMENT_H

#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>//For CentOS6.2
#include "bio_sequence.h"

namespace Cba {

class BioSeqAlignment {
public:
	BioSeqAlignment(const BioSequence* seq1, const BioSequence* seq2);

	static const int GAP = -1; // -1 indicates a gap in alignment
	typedef std::pair<int,int> ResiduePair;

	void add(int pos1, int pos2); // adds a new pair of residue pair
	void reverse(); // reverses the order of pairs in alignment
	void set_score(int s);

	int length() const; // length of alignment (including gaps)
	int nmatches() const; // number of matches in alignment
	int ngaps() const; // number of gaps in alignment
	const ResiduePair& at(size_t i) const; // positions of first and second sequences
	int get_score() const;
	double get_seq_identity() const;

	void print(std::ostream& os) const;
	static const int LINE_SEQ_FOLD = 70;

private:
	const BioSequence* m_seq1;
	const BioSequence* m_seq2;
	std::vector<ResiduePair> m_alignment;
	int m_score;

	void print_range(std::ostream& os, size_t pos, size_t len) const;
	void print_scale(std::ostream& os, size_t len) const;
	void print_seq(int seq_i, std::ostream& os, size_t pos, size_t len) const;
	int residue_number_at(size_t pos, int seq_i) const;
	void print_star(std::ostream& os, size_t pos, size_t len) const;
		// '*' indicates a pair of identical residues
};

inline void BioSeqAlignment::add(int pos1, int pos2)
{ m_alignment.push_back(std::make_pair<int,int>(pos1, pos2)); }
inline void BioSeqAlignment::reverse()
{ std::reverse(m_alignment.begin(), m_alignment.end()); }
inline void BioSeqAlignment::set_score(int s) { m_score = s; }
inline int BioSeqAlignment::length() const { return m_alignment.size(); }
inline int BioSeqAlignment::ngaps() const { return length() - nmatches(); }
inline const std::pair<int,int>& BioSeqAlignment::at(size_t i) const
{ return m_alignment.at(i); }
inline int BioSeqAlignment::get_score() const { return m_score; }

}
#endif
