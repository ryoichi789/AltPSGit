#ifndef CBA_BIO_SEQ_CMP_DP_H
#define CBA_BIO_SEQ_CMP_DP_H

#include <iostream>
#include <vector>
#include <utility>
#include "bio_seq_cmp.h"
#include "bio_sequence.h"
#include "bio_seq_score_matrix.h"
#include "bio_seq_alignment.h"

namespace Cba {

class BioSeqCmpDp : public BioSeqCmp {
public:
	BioSeqCmpDp(int n = 1, int min_score = 0,
		const char* score_file = 0, bool dna = false);
		// n = number of alignments to find
		// min_score = to find alignment with score >= min_score
		// score_file = file_containing amino acid substitution matrix (score matrix)
		//	(if 0, use default matrix)
		// dna = set true when comparing DNA sequences
	~BioSeqCmpDp();
	void compare(
		const BioSequence* seq1, const BioSequence* seq2, std::ostream& os);

private:
	int m_max_nalignments; // max no. of alignments to find from a pair of sequences
	int m_score_cutoff; // min acceptable alignment score
	const char* m_score_file;
	bool m_dna; // true when comparing DNA sequences
	BioSeqScoreMatrix m_score_matrix;

	const BioSequence* m_seq1;
	const BioSequence* m_seq2;
	int m_lseq1; // length of first sequence
	int m_lseq2; // length of second sequence

	class PathMatrix {
	public:
		PathMatrix() : m_nrows(0), m_ncolumns(0) {}
		void resize(size_t nrows, size_t ncolumns);
		char& at(size_t i, size_t j)
			{ return m_matrix[i*m_ncolumns+j]; }
		const char& at(size_t i, size_t j) const
			{ return m_matrix[i*m_ncolumns+j]; }
	private:
		size_t m_nrows; // number of rows
		size_t m_ncolumns; // number of columns
		std::vector<char> m_matrix;
	};

	PathMatrix m_path_matrix; // matrix of dimension (m_lseq1+1)x(m_lseq2+1)
	std::vector<int> m_scores[2];
		// alignment scores for two consecutive rows of path matrix

	int m_max_score;
	std::pair<int, int> m_max_score_pos;
	std::vector<BioSeqAlignment> m_alignments;

	// signs assigned to each position in path matrix
	static const char STEP_MATCH = 'm'; // match
	static const char STEP_DOWN = 'd'; // gap in second sequence
	static const char STEP_LEFT = 'l'; // gap in first sequence
	static const char USED = 'x'; // position already used by a previous alignment

	// functions used in cmopare()
	void init_dp(const BioSequence* seq1, const BioSequence* seq2);
	void prep_dp();
	int scan_path_matrix();
	void generate_alignment();
	void print(std::ostream& os);

	// functions used in scan_path_matrix()
	void step_to_row(size_t i); // i>= 1
	void step_to_pos(size_t i, size_t j); // i>=1, j>=1
	int insert_gap_down(size_t i, size_t j); // i>=1, j>=1
	int insert_gap_left(size_t i, size_t j); // i>=1, j>=1
};

}
#endif
