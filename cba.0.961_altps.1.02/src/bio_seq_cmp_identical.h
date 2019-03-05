#ifndef CBA_BIO_SEQ_CMP_IDENTICAL_H
#define CBA_BIO_SEQ_CMP_IDENTICAL_H

#include <iostream>
#include <vector>
#include "bio_seq_cmp.h"
#include "bio_seq_equiv_region.h"
#include "bio_sequence.h"

namespace Cba {

class BioSeqCmpIdentical : public BioSeqCmp {
public:
	BioSeqCmpIdentical(size_t min_len_region);
	~BioSeqCmpIdentical();
	void compare(
		const BioSequence* seq1, const BioSequence* seq2, std::ostream& os);

private:
	size_t m_min_len_region;
	const BioSequence* m_seq1;
	const BioSequence* m_seq2;
	std::vector<BioSeqEquivRegion> m_equiv_regions;

	void compare_from(int s1, int s2, int e1, int e2, int l1, int l2);
	void print(std::ostream& os);
};

}
#endif
