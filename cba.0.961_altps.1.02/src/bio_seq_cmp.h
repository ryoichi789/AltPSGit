#ifndef CBA_BIO_SEQ_CMP_H
#define CBA_BIO_SEQ_CMP_H

#include <iostream>
#include "bio_sequence.h"

namespace Cba {

class BioSeqCmp {
public:
	virtual ~BioSeqCmp() = 0;
	virtual void compare(
		const BioSequence* seq1, const BioSequence* seq2, std::ostream& os) = 0;
};

}
#endif
