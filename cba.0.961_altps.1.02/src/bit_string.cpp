#include <cstring>
#include <iostream>
#include "bit_string.h"

using namespace Cba;

const int BitString::BITS_PER_UCHAR = BitString::BITS_PER_BYTE * sizeof(unsigned char);

//
BitString::BitString() : m_bits(0), m_len(0), m_size(0) {}

//
BitString::BitString(size_t l)
{
	m_len = l;
	if (l == 0) {
		m_bits = 0;
		m_size = 0;
	} else {
		m_size = (l + BITS_PER_UCHAR - 1) / BITS_PER_UCHAR;
		m_bits = new unsigned char[m_size];
		clear();
	}
}

//
BitString::BitString(const char* s)
{
	m_len = std::strlen(s);
	m_size = (m_len + BITS_PER_UCHAR - 1) / BITS_PER_UCHAR;
	m_bits = new unsigned char[m_size];
	for (int i = 0; i < m_len; i++)
		if (s[i] == '1') set_bit(i, 1);
}

//
BitString::~BitString()
{
	delete [] m_bits;
}

//
BitString::BitString(const BitString& bs)
{
	m_len = bs.m_len;
	m_bits = new unsigned char[m_size = bs.m_size];
	for (size_t i = 0; i < m_size; i++)
		m_bits[i] = bs.m_bits[i];
}

//
BitString& BitString::operator=(const BitString& bs)
{
	if (this != &bs) {
		m_len = bs.m_len;
		if (m_size != bs.m_size) {
			delete [] m_bits;
			m_bits = new unsigned char[m_size = bs.m_size];
		}
		for (size_t i = 0; i < m_size; i++)
			m_bits[i] = bs.m_bits[i];
	}
	return *this;
}

//
BitString& BitString::operator=(const char* s)
{
	int l = std::strlen(s);
	init(l);
	for (int i = 0; i < l; i++)
		if (s[i] == '1') set_bit(i, 1);
	return *this;
}

//
bool BitString::operator==(const BitString& bs) const
{
	if (m_len != bs.m_len) return false;
	for (int i = 0; i < m_size; i++)
		if (m_bits[i] != bs.m_bits[i]) return false;
	return true;
}

// clear:
void BitString::clear()
{
	for (size_t i = 0; i < m_size; i++) m_bits[i] = 0;
}

// init:
void BitString::init(size_t l)
{
	if (l == 0) {
		delete [] m_bits;
		m_bits = 0;
		return;
	}

	m_len = l;
	size_t size0 = m_size;
	m_size = (l + BITS_PER_UCHAR - 1) / BITS_PER_UCHAR;
	if (m_size != size0) {
		delete [] m_bits;
		m_bits = new unsigned char[m_size];
	}

	clear();
}

// set_bit:
BitString& BitString::set_bit(size_t i, bool bit)
{
	size_t byte_pos, bit_pos;
	if (locate_target_bit(i, &byte_pos, &bit_pos)) {
		unsigned char byte = (1 << bit_pos);
		if (bit) {
			m_bits[byte_pos] |= byte; // set i-th bit to 1
		} else {
			byte = ~byte;
			m_bits[byte_pos] &= byte; // set i-th bit to 0
		}
	}
	return *this;
}

// locate_target_bit:
bool BitString::locate_target_bit(size_t i, size_t* byte_pos, size_t* bit_pos) const
{
	if (i >= m_len) return false;
	*byte_pos = i / BITS_PER_UCHAR; // target byte within which target bit is located
	*bit_pos = i % BITS_PER_UCHAR; // target bit position in target byte
	return true;
}

// get_bit:
int BitString::get_bit(size_t i) const
{
	size_t byte_pos, bit_pos;
	if (! locate_target_bit(i, &byte_pos, &bit_pos)) return -1;
	unsigned char byte = m_bits[byte_pos];
	return ((byte >> bit_pos) & 1);
}

// count_ons:
long BitString::count_ons() const
{
	size_t nons = 0;
	for (size_t i = 0; i < m_size; i++)
		nons += count_ons_in_byte(m_bits[i]);
	return nons;
}

// count_offs:
long BitString::count_offs() const
{
	return m_len - count_ons();
}

// count_ons_in_byte:
size_t BitString::count_ons_in_byte(unsigned char byte)
{
	size_t nons = 0;
	for (size_t i = 0; i < BITS_PER_UCHAR; i++)
		if ((byte >> i) & 1) nons++;
	return nons;
}

// bs_and:
void BitString::bs_and(const BitString& bs)
{
	if (m_len != bs.m_len) return;
	for (size_t i = 0; i < m_size; i++)
		m_bits[i] = (m_bits[i] & bs.m_bits[i]);
}

// bs_or:
void BitString::bs_or(const BitString& bs)
{
	if (m_len != bs.m_len) return;
	for (size_t i = 0; i < m_size; i++)
		m_bits[i] = (m_bits[i] | bs.m_bits[i]);
}

// bs_xor:
void BitString::bs_xor(const BitString& bs)
{
	if (m_len != bs.m_len) return;
	for (size_t i = 0; i < m_size; i++)
		m_bits[i] = (m_bits[i] ^ bs.m_bits[i]);
}

// flip:
void BitString::flip(size_t i)
{
	int val = get_bit(i);
	if (val == 1) set_bit(i, 0);
	else if (val == 0) set_bit(i, 1);
}

// flip:
void BitString::flip()
{
	for (size_t i = 0; i < m_size; i++)
		m_bits[i] = ~(m_bits[i]);
	// take care of trailing bits (that should always be 0)
	if (size_t i = m_len % BITS_PER_UCHAR)
		m_bits[m_size - 1] &= ~(~0 << i);
}

// str:
std::string BitString::str() const
{
	std::string s;
	int bit_val;
	for (int i = 0; i < m_len; i++) {
		bit_val = get_bit(i);
		if (bit_val == 1) s.append("1");
		else if (bit_val == 0) s.append("0");
		else s.append("-");
	}
	return s;
}

//
void BitString::write(std::ostream& os) const
{
	const char* c;

	// length
	c = (const char*)&m_len;
	os.write(c, sizeof(size_t));

	// bitstring
	os.write((const char*)m_bits, m_size);
}

//
bool BitString::read(std::istream& is)
{
	int l;
	char* c = (char*)&l;
	if (! is.read(c, sizeof(int))) return false;
	init(l);
	if (! is.read((char*)m_bits, m_size)) return false;
	return true;
}

//
double BitString::get_tanimoto_coeff(const BitString& bs) const
{
	if (m_len != bs.m_len) return -1;

	int n1 = 0;
	int n2 = 0;
	int n12 = 0;
	int s1, s2; // bit status
	for (int i = 0; i < m_len; i++) {
		if ((s1 = get_bit(i)) == 1) n1++;
		if ((s2 = bs.get_bit(i)) == 1) n2++;
		if (s1 == 1 && s2 == 1) n12++;
	}
	if (n12 == 0) return 0.0;
	return (double)n12 / (double)(n1 + n2 - n12);
}
