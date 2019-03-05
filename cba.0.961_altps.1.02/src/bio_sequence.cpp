#include <iostream>
#include <string>
#include "use_boost.h"
#include "bio_sequence.h"
#include "region.h"

using namespace Cba;

// ~BioSequence:
BioSequence::~BioSequence() {}

// clear:
void BioSequence::clear()
{
	m_id.erase();
	m_annot.erase();
	m_seq.erase();
}

// print:
void BioSequence::print(std::ostream& os, int line_length) const
{
	os << '>' << id() << ' ' << annot() << std::endl;

	// print sequence
	int i = 0; // i-th position in sequence
	int p; // position in line
	while (i < length()) {
		os << at(i++);
		p = i % line_length;
		if (p == 0) os << std::endl;
	}
	if (p) os << std::endl;
}

// find_motif:
int BioSequence::find_motif(const char* motif, Regions& regions)
{
	using std::string;
	regions.clear();

#ifdef CBA_USE_BOOST
	boost::regex q(motif);
	boost::match_results<string::iterator> m;
	string::iterator start0 = m_seq.begin();
	for (string::iterator start = start0;
		boost::regex_search(start, m_seq.end(), m, q); start = m[0].first+1) {
		regions.add((m[0].first-start0), (m[0].second-start0)-1);
	}
#endif

	return regions.size();
}
