#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <utility>
#include "bio_seq_db_search.h"
#include "bio_sequence.h"
#include "fasta_reader.h"
#include "bio_seq_score_matrix.h"
#include "bio_seq_alignment.h"

using namespace Cba;

const char* BioSeqDbSearch::SEGPOS_FILE_EXT = ".sg";

//
BioSeqDbSearch::BioSeqDbSearch()
	: m_query_fasta_file(0), m_target_fasta_file(0) {}

//
BioSeqDbSearch::~BioSeqDbSearch() {}

// create_tables:
void BioSeqDbSearch::create_tables(const char* target)
{
	if (!target) return;
	set_target(target);

	std::ifstream ifs(m_target_fasta_file);
	FastaReader fas(ifs);
	BioSequence seq;

	std::FILE* ofp = std::fopen(m_segpos_file.c_str(), "wb");
	int seqblock_size = 0;
	int seqpos = ifs.tellg();
	for (int iseq = 0; fas.read(seq); iseq++) {
		update_segpos_table(seqpos, seq);
		seqblock_size += seq.length();
		if (seqblock_size > SEQBLOCK_SIZE) {
			write_segpos_table(ofp);
			m_segpos_table.clear();
			seqblock_size = 0;
		}
		seqpos = ifs.tellg();
	}
	if (seqblock_size) write_segpos_table(ofp);
	std::fclose(ofp);
}

// set_target:
void BioSeqDbSearch::set_target(const char* target)
{
	m_target_fasta_file = target;
	m_segpos_file = m_target_fasta_file;
	m_segpos_file += SEGPOS_FILE_EXT;
}

// update_segpos_table:
void BioSeqDbSearch::update_segpos_table(int seqpos, const BioSequence& seq)
{
	const char* s = seq.str().c_str();
	std::string seg;
	int nsegs = seq.length() - SEGMENT_LENGTH + 1;
	for (int i = 0; i < nsegs; i++) {
		seg.assign(s++, SEGMENT_LENGTH);
		m_segpos_table[seg].push_back(seqpos+i);
	}
}

// write_segpos_table:
void BioSeqDbSearch::write_segpos_table(std::FILE* ofp) const
{
	using std::fwrite;

	int n;
	for (SegPosTable::const_iterator i = m_segpos_table.begin();
		i != m_segpos_table.end(); ++i) {
		// segment
		fwrite(i->first.c_str(), sizeof(char), SEGMENT_LENGTH, ofp);

		// number of segment positions
		n = i->second.size();
		fwrite(&n, sizeof(int), 1, ofp);

		// segment positions
		for (VICI j = i->second.begin(); j != i->second.end(); ++j) {
			n = *j; // segment position
			fwrite(&n, sizeof(int), 1, ofp);
		}
	}
}

// search:
void BioSeqDbSearch::search(const char* query, const char* target,
	int smin, int lmin, const char* output)
{
	if (!(query && target)) return;
	m_query_fasta_file = query;
	set_target(target);
	m_score_min = smin;
	m_length_min = lmin;

	std::ostream* os;
	if (output) os = new std::ofstream(output);
	else os = &std::cout;

	std::ifstream ifs(m_query_fasta_file);
	FastaReader fas(ifs);
	BioSequence seq;

	// for each query sequence
	std::FILE* ifp = std::fopen(m_segpos_file.c_str(), "rb");
	while (fas.read(seq)) {
		get_segments_in_query(seq); // query segments
		std::rewind(ifp);
		read_segpos_table(ifp); // positions of query segments in target
		get_segment_matches(seq);
		get_segment_scores();
		scan_db(seq, *os);
	}
	std::fclose(ifp);
	if (os != &std::cout) delete os;
}

// get_segments_in_query:
void BioSeqDbSearch::get_segments_in_query(const BioSequence& seq)
{
	m_segpos_table.clear();

	const char* s = seq.str().c_str();
	std::string seg;
	int nsegs = seq.length() - SEGMENT_LENGTH + 1;
	for (int i = 0; i < nsegs; i++) {
		seg.assign(s++, SEGMENT_LENGTH);
		m_segpos_table[seg];
			// add a pair of segment (key) and empty position list
	}
}

// read_segpos_table:
void BioSeqDbSearch::read_segpos_table(std::FILE* ifp)
{
	using std::fread;

	char seg[SEGMENT_LENGTH+1];
	seg[SEGMENT_LENGTH] = '\0';
	int nsegs;
	int pos;
	SegPosTable::iterator spti;
	while (fread(seg, sizeof(char), SEGMENT_LENGTH, ifp) == SEGMENT_LENGTH) {
		fread(&nsegs, sizeof(int), 1, ifp); // number of segment positions
		spti = m_segpos_table.find(seg);
		if (spti != m_segpos_table.end()) { // query has this segment 
			for (int i = 0; i < nsegs; i++) {
				fread(&pos, sizeof(int), 1, ifp);
				spti->second.push_back(pos);
			}
		} else { // skip to next segment
			std::fseek(ifp, nsegs*sizeof(int), SEEK_CUR);
		}
	}
}

// get_segment_matches:
void BioSeqDbSearch::get_segment_matches(const BioSequence& qseq)
{
	int nsegs = qseq.length() - SEGMENT_LENGTH + 1;
	m_segment_matches.reserve(nsegs);
	m_segment_matches.clear();
	for (int i = 0; i < nsegs; i++) { // for i-th query segment
		const std::vector<int>& v
			= m_segpos_table[qseq.str().substr(i,SEGMENT_LENGTH)];
		m_segment_matches.push_back(std::make_pair<VICI,VICI>(v.begin(),v.end()));
			// range of identical target segments
	}
}

// get_segment_scores:
void BioSeqDbSearch::get_segment_scores()
{
	m_segment_scores.clear();
	for (SegPosTable::const_iterator i = m_segpos_table.begin();
		i != m_segpos_table.end(); ++i) {
		m_segment_scores[i->first] = get_segment_score(i->first);
	}
}

// get_segment_score:
int BioSeqDbSearch::get_segment_score(const std::string& seg)
{
	int s = 0; // score
	char a; // amino acid
	for (int i = 0; i < SEGMENT_LENGTH; i++) {
		a = seg.at(i);
		s += m_score_matrix.score(a, a);
	}
	return s;
}

// scan_db:
void BioSeqDbSearch::scan_db(const BioSequence& qseq, std::ostream& os)
{
	std::ifstream ifs(m_target_fasta_file);
	FastaReader fas(ifs);
	BioSequence tseq;

	// for each target sequence
	int seqpos = ifs.tellg();
	while (fas.read(tseq)) {
		compare(qseq, tseq, seqpos, os);
		seqpos = ifs.tellg();
	}
}

// compare:
void BioSeqDbSearch::compare(const BioSequence& qseq, const BioSequence& tseq,
	int tpos0, std::ostream& os)
{
	int nqsegs = qseq.length() - SEGMENT_LENGTH + 1; // no. of query segments
	int ntsegs = tseq.length() - SEGMENT_LENGTH + 1; // no. of target segments
	int tpos1 = tpos0 + ntsegs; // last target segment position + 1
	m_matching_regions.clear();

	VICI i;
	for (int qseg = 0; qseg < nqsegs; qseg++) { // for each query segment
		for (i = m_segment_matches[qseg].first;
			i != m_segment_matches[qseg].second && *i < tpos1; ++i) {
			extend_segment(qseq, qseg, tseq, *i-tpos0);
		}
		m_segment_matches[qseg].first = i; // next target segment
	}

	for (std::vector<MatchingRegion>::const_iterator
		j = m_matching_regions.begin(); j != m_matching_regions.end(); ++j)
		if (j->length >= m_length_min && j->score >= m_score_min)
			print_matching_region(qseq, tseq, *j, os);
}

// extend_segment:
void BioSeqDbSearch::extend_segment(
	const BioSequence& qseq, int q0, const BioSequence& tseq, int t0)
{
	if (already_matched(q0, t0)) return;
	int s = m_segment_scores[qseq.str().substr(q0,SEGMENT_LENGTH)];

	// N-terminal <---
	int sn = s;
	int sn_max = s;
	int q1_max = q0;
	int q1, t1;
	for (q1=q0-1, t1=t0-1; q1>=0 && t1>=0; q1--, t1--) {
		sn += m_score_matrix.score(qseq.at(q1), tseq.at(t1));
		if (sn > sn_max) {
			sn_max = sn;
			q1_max = q1;
		}
		if (sn <= 0) break;
	}

	// ---> C-terminal
	int sc = s;
	int sc_max = s;
	int q2_max = q0 + SEGMENT_LENGTH -1;
	int q2, t2;
	for (q2=q0+SEGMENT_LENGTH, t2=t0+SEGMENT_LENGTH;
		q2<qseq.length() && t2<tseq.length(); q2++, t2++) {
		sc += m_score_matrix.score(qseq.at(q2), tseq.at(t2));
		if (sc > sc_max) {
			sc_max = sc;
			q2_max = q2;
		}
		if (sc <= 0) break;
	}

	s = sn_max + sc_max - s;
	m_matching_regions.push_back(MatchingRegion(
		q1_max, t0+(q1_max-q0), q2_max-q1_max+1, s));
}

// already_matched:
bool BioSeqDbSearch::already_matched(int q, int t) const
{
	int p = q - t;
	for (std::vector<MatchingRegion>::const_iterator
		i = m_matching_regions.begin(); i != m_matching_regions.end(); ++i) {
		if ((i->q0 - i->t0) == p  && q < (i->q0 + i->length))
			return true;
	}
	return false;
}

// print_matching_region:
void BioSeqDbSearch::print_matching_region(
	const BioSequence& qseq, const BioSequence& tseq,
	const MatchingRegion& region, std::ostream& os) const
{
	BioSeqAlignment alignment(&qseq, &tseq);
	alignment.set_score(region.score);
	for (int i = 0; i < region.length; i++)
		alignment.add(region.q0+i, region.t0+i);
	alignment.print(os);
}
