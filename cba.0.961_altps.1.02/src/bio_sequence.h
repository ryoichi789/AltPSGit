#ifndef CBA_BIO_SEQUENCE_H
#define CBA_BIO_SEQUENCE_H

#include <string>
#include <iostream>
#include "region.h"

namespace Cba {

class BioSequence {
public:
	virtual ~BioSequence();
	virtual void clear();
	void set_id(const std::string& id);
	void set_annot(const std::string& annot);
	void add(char c);
	void add(const std::string& seq);
	const std::string& id() const;
	const std::string& annot() const;
	int length() const;
	char at(size_t i) const;
	const std::string& str() const;
	void print(std::ostream& os, int line_length = FASTA_LINE_LENGTH) const;
	static const int FASTA_LINE_LENGTH = 70;

	int find_motif(const char* motif, Regions& regions);

private:
	std::string m_id;
	std::string m_annot;
	std::string m_seq;
};

inline void BioSequence::set_id(const std::string& id) { m_id = id; }
inline void BioSequence::set_annot(const std::string& annot) { m_annot = annot; }
inline void BioSequence::add(char c) { m_seq += c; }
inline void BioSequence::add(const std::string& seq) { m_seq += seq; }
inline const std::string& BioSequence::id() const { return m_id; }
inline const std::string& BioSequence::annot() const { return m_annot; }
inline int BioSequence::length() const { return m_seq.length(); }
inline char BioSequence::at(size_t i) const { return m_seq.at(i); }
inline const std::string& BioSequence::str() const { return m_seq; }

}

#endif
