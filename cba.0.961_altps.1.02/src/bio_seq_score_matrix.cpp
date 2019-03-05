#include "bio_seq_score_matrix.h"

using namespace Cba;

const char* BioSeqScoreMatrix::amino_acid_code = "ARNDCQEGHILKMFPSTWYVBZX";
const int BioSeqScoreMatrix::blosum62
	[BioSeqScoreMatrix::N_AMINO_TYPES][BioSeqScoreMatrix::N_AMINO_TYPES]
= {
{ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0},
{-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1},
{-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1},
{-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1},
{ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2},
{-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1},
{-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1},
{ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1},
{-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1},
{-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1},
{-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1},
{-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1},
{-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1},
{-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1},
{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2},
{ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0},
{ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0},
{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2},
{-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1},
{ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1},
{-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1},
{-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1},
{ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1}
};

const char* BioSeqScoreMatrix::nucleotide_code = "ACGTMRWSYKVHDBN";
const int BioSeqScoreMatrix::nuc_matrix
	[BioSeqScoreMatrix::N_NUCLEOTIDE_TYPES][BioSeqScoreMatrix::N_NUCLEOTIDE_TYPES]
= {
{ 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
{-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1},
{-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1},
{-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1},
{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0}
};

// BioSeqScoreMatrix
BioSeqScoreMatrix::BioSeqScoreMatrix()
	: m_for_amino_acids(true)
{
	// score matrix used by default
	for (int i = 0; i < N_AMINO_TYPES; i++)
		for (int j = 0; j < N_AMINO_TYPES; j++)
			m_score_matrix[i][j] = blosum62[i][j];
	m_gap_open = -8;
	m_gap_ext = -2;
}

// ~BioSeqScoreMatrix
BioSeqScoreMatrix::~BioSeqScoreMatrix() {}

// for_nucleotides:
void BioSeqScoreMatrix::for_nucleotides()
{
	m_for_amino_acids = false;
	for (int i = 0; i < N_NUCLEOTIDE_TYPES; i++)
		for (int j = 0; j < N_NUCLEOTIDE_TYPES; j++)
			m_score_matrix[i][j] = nuc_matrix[i][j];
	m_gap_open = -2;
	m_gap_ext = -1;
}

// get_score:
int BioSeqScoreMatrix::score(char a1, char a2) const
{
	return m_score_matrix[amino_nucleotide_number(a1)][amino_nucleotide_number(a2)];
}

// amino_nucleotide_number:
size_t BioSeqScoreMatrix::amino_nucleotide_number(char a) const
{
	int i;
	const char* c = (m_for_amino_acids ? amino_acid_code : nucleotide_code);
	for (i = 0; *c != '\0'; i++, c++)
		if (a == *c) return i;
	return i-1; // as 'X' or 'N'
}

// read:
void BioSeqScoreMatrix::read(std::istream& is)
{
	// to be implemented:
	//
	// * assume a certain file format in which a score matrix is given;
	//   format should conform to the order in which amino acids or nucleotides
	//   are listed by const char* amino_acid_codes and const char* nucleotide_codes;
	// * determine whether the matrix is for amino acids or nucleotides
	//   from the format;
	// * for_nucleotides() if it is for nucleotides;
	// * read data into m_score_matrix, m_gap_open, and m_gap_ext
}
