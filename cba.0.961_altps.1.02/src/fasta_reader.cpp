#include <iostream>
#include <string>
#include <cctype>
#include "fasta_reader.h"
#include "bio_sequence.h"
using namespace Cba;
using std::string;

// FastaReader():
FastaReader::FastaReader(std::istream& is) : m_is(is) {}

// read: read first entry from Fasta input stream
bool FastaReader::read(BioSequence& seq)
{
	seq.clear();
	if (read_header(seq)) {
		read_seq(seq);
		return true;
	}
	return false;
}

// read_header: read header to set id and annotation
bool FastaReader::read_header(BioSequence& seq)
{
	string line;
	while (std::getline(m_is, line)) {
		if (line.length() && line[0] == '>') {
			read_id_annot(line, seq);
			return true;
		}
	}
	return false;
}

// read_id_annot: read id and annotation from header line
void FastaReader::read_id_annot(const string& header, BioSequence& seq)
{
	int l = header.length();
	if (l <= 1) return; // no id/annotation

	// id
	int i = header.find_first_of(" \t", 1);
	if (i == string::npos) i = l;
	seq.set_id(header.substr(1, i-1));

	// annotation
	if (i == l) return; // no annotation
	int j = header.find_first_not_of(" \t", i+1);
	if (j != string::npos) seq.set_annot(header.substr(j, l-j));
}

// read_seq:
void FastaReader::read_seq(BioSequence& seq)
{
	string s;
	char c;
	bool hol = true; // true if head of line
	while (m_is.get(c)) {
		if (hol && c == '>') {
			m_is.unget();
			break;
		}
		if (!std::isspace(c)) s += c;
		hol = (c == '\n'); // next is head of line
	}
	seq.add(s);
}
