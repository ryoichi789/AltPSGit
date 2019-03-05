#ifndef CBA_DNA_SEQUENCE_H
#define CBA_DNA_SEQUENCE_H

#include <string>
#include <vector>
#include <iostream>
#include <functional>
#include "bio_sequence.h"

namespace Cba {

class DnaSequence : public BioSequence {
public:
	~DnaSequence();
	void clear();
	int find_orfs(); // returns the number of ORFs found;
	void print_orfs(std::ostream& os, int min_len = 1) const;
	bool translate_orf(size_t orf_i, BioSequence& seq, int min_len = 1) const;
		// returns true if orf_i-th ORF is of length >= min_len
	char complement_at(size_t i); // complementary nucleotide at i (0 for invalid i)

private:
	struct Orf {
		int m_reading_frame; // +1,+2,+3, and -1,-2,-3 for complementary chain
		int m_start; // starting position in DNA sequence
		int m_end; // ending position in DNA sequence
		int length() const { return (m_end-m_start+1)/3; }
	};

	// to sort ORFs in order of lengths
	struct OrfCmp : std::binary_function<Orf, Orf, bool> {
		bool operator()(const Orf& a, const Orf& b) const
		{ return a.length() > b.length(); }
	};

	// to print ORF
	class OrfPrint : public std::unary_function<Orf, void> {
	public:
		OrfPrint(std::ostream& os, size_t min_len)
			: m_os(os), m_min_len(min_len) {}
		void operator()(const Orf& orf);
	private:
		std::ostream& m_os;
		int m_min_len;
	};

	std::vector<Orf> m_orfs;
	std::string m_c_seq; // complementary sequence
	std::vector<int> m_stop_codons[6]; // position of stop codons in six frames

	void make_c_seq(); // makes complementary sequence
	void find_stop_codons(); // in all six reading frames
	void find_orfs_in_reading_frame(int frame_i);
	void sort_orfs(); // sort ORFs in order of lengths
	void clear_data();
};

}
#endif
