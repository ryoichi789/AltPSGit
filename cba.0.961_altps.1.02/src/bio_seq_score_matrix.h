#ifndef CBA_BIO_SEQ_SCORE_MATRIX_H
#define CBA_BIO_SEQ_SCORE_MATRIX_H

#include <vector>
#include <istream>//For CentOS6.2

namespace Cba {

class BioSeqScoreMatrix {
public:
	BioSeqScoreMatrix();
	~BioSeqScoreMatrix();

	int score(char a1, char a2) const; // returns score for a1 versus a2
	void read(std::istream& is); // reads a score matrix
	int gap_open() const; // returns score for opening a gap
	int gap_ext() const; // returns score for extending a gap
	void for_nucleotides(); // to use a default score matrix for nucleotides

	static const size_t N_AMINO_TYPES = 23; // 20 amino acids + B,Z,X
	static const char* amino_acid_code; // 1-letter code of amino acids
	static const int blosum62[N_AMINO_TYPES][N_AMINO_TYPES];
		// score matrix used by default for amino acids

	static const size_t N_NUCLEOTIDE_TYPES = 15; // ACGTMRWSYKVHDBN (U = T)
	static const char* nucleotide_code; // 1-letter code of nucleotides
	static const int nuc_matrix[N_NUCLEOTIDE_TYPES][N_NUCLEOTIDE_TYPES];
		// score matrix used by default for nucleotides

private:
	static const size_t MAX_N_TYPES = N_AMINO_TYPES;
		// largest of N_AMINO_TYPES and N_NUCLEOTIDE_TYPES;
		// obviously, it is always N_AMINO_TYPES;
	int m_score_matrix[MAX_N_TYPES][MAX_N_TYPES];
	int m_gap_open; // score for opening a gap
	int m_gap_ext; // score for extending a gap
	bool m_for_amino_acids; // true if matrix is for amino acids (default)

	size_t amino_nucleotide_number(char a) const;
};

inline int BioSeqScoreMatrix::gap_open() const { return m_gap_open; }
inline int BioSeqScoreMatrix::gap_ext() const { return m_gap_ext; }

}
#endif
