#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include "bio_seq_cmp_dp.h"
#include "bio_sequence.h"
#include "bio_seq_score_matrix.h"
using namespace Cba;

// BioSeqCmpDp:
BioSeqCmpDp::BioSeqCmpDp(int n, int min_score, const char* score_file, bool dna)
	: m_max_nalignments(n), m_score_cutoff(min_score),
		m_score_file(score_file), m_dna(dna)
{
	if (m_score_file) {
		std::ifstream ifs(score_file);
		m_score_matrix.read(ifs);
	} else if (m_dna) {
		m_score_matrix.for_nucleotides();
	}
}

// ~BioSeqCmpDp:
BioSeqCmpDp::~BioSeqCmpDp() {}

// compare:
void BioSeqCmpDp::compare(
	const BioSequence* seq1, const BioSequence* seq2, std::ostream& os)
{
	init_dp(seq1, seq2);
	for (int i = 0; i < m_max_nalignments; i++) {
		prep_dp();
		if (scan_path_matrix() < m_score_cutoff) break;
		generate_alignment();
	}

	print(os);
}

// init_dp:
void BioSeqCmpDp::init_dp(const BioSequence* seq1, const BioSequence* seq2)
{
	m_seq1 = seq1;
	m_seq2 = seq2;
	m_lseq1 = seq1->length();
	m_lseq2 = seq2->length();

	m_path_matrix.resize(m_lseq1+1, m_lseq2+1);
	for (int i = 0; i <= m_lseq1; i++)
		for (int j = 0; j <= m_lseq2; j++)
			m_path_matrix.at(i, j) = ' ';
	for (int i = 0; i < 2; i++)
		m_scores[i].resize(m_lseq2+1);

	m_path_matrix.at(0,0) = USED; // origin
	m_scores[0][0] = 0;
	m_alignments.clear();
}

// prep_dp:
void BioSeqCmpDp::prep_dp()
{
	for (int i = 1; i <= m_lseq1; i++)
		for (int j = 1; j <= m_lseq2; j++)
			if (m_path_matrix.at(i,j) != USED) m_path_matrix.at(i,j) = ' ';

	for (int j = 1; j <= m_lseq2; j++) { // for first row
		m_path_matrix.at(0,j) = USED;
		m_scores[0][j] = 0;
	}

	m_max_score = 0;
	m_max_score_pos.first = m_max_score_pos.second = 0;

}

// scan_path_matrix:
int BioSeqCmpDp::scan_path_matrix()
{
	for (size_t i = 1; i <= m_lseq1; i++) {
		// go to row i
		step_to_row(i);
		// before going to next row...
		for (int j = 0; j <= m_lseq2; j++)
			m_scores[0][j] = m_scores[1][j];
	}

	return m_max_score;
}

// step_to_row:
void BioSeqCmpDp::step_to_row(size_t i)
{
	// for (i,0)
	m_path_matrix.at(i,0) = USED;
	m_scores[1][0] = 0;

	// for (i,j) (j>=1)
	for (size_t j = 1; j <= m_lseq2; j++) {
		step_to_pos(i, j);

		// update position of max score
		if (m_scores[1][j] > m_max_score) {
			m_max_score_pos.first = i;
			m_max_score_pos.second = j;
			m_max_score = m_scores[1][j];
		}
	}
}

// step_to_pos:
void BioSeqCmpDp::step_to_pos(size_t i, size_t j)
{
	if (m_path_matrix.at(i,j) == USED) {
		m_scores[1][j] = 0; // score for (i,j)
		return;
	}

	// calculate scores for STEP_MATCH, STEP_DOWN, STEP_LEFT
	int score_m = m_scores[0][j-1] +
		m_score_matrix.score(m_seq1->at(i-1), m_seq2->at(j-1));
	int score_d = m_scores[0][j] + insert_gap_down(i,j);
	int score_l = m_scores[1][j-1] + insert_gap_left(i,j);

	// find best step of STEP_MATCH, STEP_DOWN, STEP_LEFT
	int score_max = score_m;
	char istep = STEP_MATCH;
	if (score_d > score_max) {
		score_max = score_d;
		istep = STEP_DOWN;
	}
	if (score_l > score_max) {
		score_max = score_l;
		istep = STEP_LEFT;
	}
	m_scores[1][j] = score_max; // score for (i,j)
	m_path_matrix.at(i,j) = istep; // step to (i,j)

	if (m_scores[1][j] <= 0) {
		m_scores[1][j] = 0;
		m_path_matrix.at(i,j) = USED;
	}
}

// insert_gap_down:
int BioSeqCmpDp::insert_gap_down(size_t i, size_t j)
{
	if (m_path_matrix.at(i-1,j) == STEP_DOWN) return m_score_matrix.gap_ext();
	return m_score_matrix.gap_open();
}

// insert_gap_left:
int BioSeqCmpDp::insert_gap_left(size_t i, size_t j)
{
	if (m_path_matrix.at(i,j-1) == STEP_LEFT) return m_score_matrix.gap_ext();
	return m_score_matrix.gap_open();
}

// generate_alignment:
void BioSeqCmpDp::generate_alignment()
{
	std::pair<int,int> p; // position in path matrix
	BioSeqAlignment alignment(m_seq1, m_seq2);

	// start backtracking from position with max score
	p.first = m_max_score_pos.first;
	p.second = m_max_score_pos.second;
	// backtrack to the first 'USED' position
	char step;
	while ((step = m_path_matrix.at(p.first,p.second)) != USED) {
		m_path_matrix.at(p.first,p.second) = USED;
			// this position will not be used in the next alignment
		if (step == STEP_MATCH)
			alignment.add((p.first--)-1, (p.second--)-1);
		else if (step == STEP_DOWN)
			alignment.add((p.first--)-1, BioSeqAlignment::GAP);
		else if (step == STEP_LEFT)
			alignment.add(BioSeqAlignment::GAP, (p.second--)-1);
			// position (i,j) of path matrix corresponds to
			// residue i-1 of sequence 1 and residue j-1 of sequence 2
	}

	alignment.reverse();
	alignment.set_score(m_max_score);
	m_alignments.push_back(alignment);
}

//
void BioSeqCmpDp::print(std::ostream& os)
{
	for (std::vector<BioSeqAlignment>::const_iterator i = m_alignments.begin();
		i != m_alignments.end(); ++i)
		i->print(os);
}

// PathMatrix::resize:
void BioSeqCmpDp::PathMatrix::resize(size_t nrows, size_t ncolumns)
{
	m_nrows = nrows;
	m_ncolumns = ncolumns;
	m_matrix.resize(m_nrows * m_ncolumns);
}
