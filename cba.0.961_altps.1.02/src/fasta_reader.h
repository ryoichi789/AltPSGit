#ifndef CBA_FASTA_READER_H
#define CBA_FASTA_READER_H

#include <iostream>
#include <string>
#include "bio_sequence.h"

namespace Cba {

class FastaReader {
public:
	FastaReader(std::istream& is);
	bool read(BioSequence& seq);

private:
	std::istream& m_is;
	bool read_header(BioSequence& seq);
	void read_seq(BioSequence& seq);
	void read_id_annot(const std::string& header, BioSequence& seq);
};

}

#endif
