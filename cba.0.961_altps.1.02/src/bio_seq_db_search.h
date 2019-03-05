#ifndef CBA_BIO_SEQ_DB_SEARCH
#define CBA_BIO_SEQ_DB_SEARCH

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <functional>
#include <iostream>
#include <cstdio>
#include "bio_sequence.h"
#include "bio_seq_score_matrix.h"

namespace Cba {

class BioSeqDbSearch {
public:
	BioSeqDbSearch();
	~BioSeqDbSearch();
	void create_tables(const char* target); // for target FASTA file
	void search(const char* query, const char* target,
		int smin, int lmin, const char* output = 0);

private:
	const char* m_query_fasta_file;
	const char* m_target_fasta_file;
	std::string m_segpos_file; // positions of segments in target FASTA file
		// position of segment =
		// position of sequence + position of segment in sequence
	int m_score_min; // minimum acceptable score for alignment
	int m_length_min; // minimum acceptable length of alignment

	typedef std::map<std::string, std::vector<int> > SegPosTable;
	typedef std::vector<int>::const_iterator VICI;
	SegPosTable m_segpos_table;
	std::vector<std::pair<VICI, VICI> > m_segment_matches;
		// target positions of query segments
	std::map<std::string, int> m_segment_scores;
		// BLOSUM62 score for each segment
	BioSeqScoreMatrix m_score_matrix;

	struct MatchingRegion {
		MatchingRegion(int q, int t, int l, int s)
			: q0(q), t0(t), length(l), score(s) {}
		int q0;
		int t0;
		int length;
		int score;
	};
	std::vector<MatchingRegion> m_matching_regions;

	static const char* SEGPOS_FILE_EXT;
	static const int SEGMENT_LENGTH = 3;
	static const int SEQBLOCK_SIZE = 8000000;

	void set_target(const char* target); // target FASTA file
	void update_segpos_table(int seqpos, const BioSequence& seq);
	void write_segpos_table(std::FILE* ofp) const;

	void get_segments_in_query(const BioSequence& seq);
	void get_segment_scores();
	int get_segment_score(const std::string& seg);
	void read_segpos_table(std::FILE* ifp);
	void get_segment_matches(const BioSequence& qseq);
	void scan_db(const BioSequence& qseq, std::ostream& os);
	void compare(const BioSequence& qseq,
		const BioSequence& tseq, int tpos0, std::ostream& os);
	void extend_segment(const BioSequence& qseq, int q0,
		const BioSequence& tseg, int t0);
	bool already_matched(int q, int t) const;
	void print_matching_region(
		const BioSequence& qseq, const BioSequence& tseq,
		const MatchingRegion& region, std::ostream& os) const;
};

}
#endif
