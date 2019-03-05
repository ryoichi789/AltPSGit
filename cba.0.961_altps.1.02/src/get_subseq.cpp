#include <iostream>
#include <fstream>
#include <cstdlib>
#include <bio_sequence.h>
#include <fasta_reader.h>

using namespace Cba;

int main(int argc, char *argv[])
{
	std::ifstream ifs(argv[1]);
	FastaReader fas(ifs);
	BioSequence seq1, seq2;
	fas.read(seq1);

	int i1 = std::atoi(argv[2]);
	int i2 = std::atoi(argv[3]);
	for (int i = i1; i <= i2; i++)
		seq2.add(seq1.at(i));

	seq2.print(std::cout);
}
