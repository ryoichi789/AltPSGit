#ifndef CBA_BIT_STRING_H
#define CBA_BIT_STRING_H

#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <string>
#include <iostream>

namespace Cba {

class BitString {
public:
	BitString();
	BitString(size_t l);
	BitString(const char* s);
	~BitString();
	BitString(const BitString& bs);
	BitString& operator=(const BitString& bs);
	BitString& operator=(const char* s);
	bool operator==(const BitString& bs) const;
	bool operator!=(const BitString& bs) const;
	void clear(); // sets all bits to 0
	void init(size_t l); // makes a bit string of length l with each bit 0
	int length() const;
	BitString& set_bit(size_t i, bool bit);
	int get_bit(size_t i) const; // returns negative for error
	long count_ons() const;
	long count_offs() const;
	void bs_and(const BitString& bs); // this <- this & bs (bitwise and)
	void bs_or(const BitString& bs); // this <- this | bs (bitwise or)
	void bs_xor(const BitString& bs); // this <- this xor bs (bitwise exclusive or)
	void flip(size_t i); // flips the i-th bit value (0->1 or 1->0)
	void flip(); // flips all bits
	std::string str() const;

	// binary output and input
	void write(std::ostream& os) const;
	bool read(std::istream& is);
	double get_tanimoto_coeff(const BitString& bs) const;

	static const int BITS_PER_BYTE = 8;
	static const int BITS_PER_UCHAR;

private:
	unsigned char *m_bits;
	size_t m_len; /* length of bit string (number of bits in string) */
	size_t m_size; /* number of chars allocated for bits */

	bool locate_target_bit(size_t i, size_t* byte_pos, size_t* bit_pos) const;
		// locates position of a target bit
		//	as 'bit_pos'-th bit in 'byte_pos'-th element of byte array 'bits'
	static size_t count_ons_in_byte(unsigned char byte);
};

inline int BitString::length() const { return m_len; }
inline bool BitString::operator!=(const BitString& bs) const { return !operator==(bs); }

}
#endif
