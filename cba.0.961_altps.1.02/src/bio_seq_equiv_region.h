#ifndef CBA_BIO_SEQ_EQUIV_REGION_H
#define CBA_BIO_SEQ_EQUIV_REGION_H
#include <functional>

namespace Cba {

struct BioSeqEquivRegion {
	size_t pos1; // start of region in first sequence
	size_t pos2; // start of region in second sequence
	size_t len; // length of region
	BioSeqEquivRegion() {}
	BioSeqEquivRegion(size_t p1, size_t p2, size_t l)
		: pos1(p1), pos2(p2), len(l) {}
};

// function object used in sorting BioSeqEquivRegion objects
struct BioSeqEquivRegionCmp
	: std::binary_function<BioSeqEquivRegion, BioSeqEquivRegion, bool> {
	bool operator() (const BioSeqEquivRegion& e1, const BioSeqEquivRegion& e2) const
	{ return e1.pos1 < e2.pos1; }
};

}
#endif
