#include <cctype>
#include "genetic_code.h"

using namespace Cba;

// decode:
char GeneticCode::decode(char c1, char c2, char c3)
{
	c1 = purify(c1);
	c2 = purify(c2);
	c3 = purify(c3);

	if (c1 == 'U') {
		///////
		// Uxx
		///////
		if (c2 == 'U') { // UUx
			if (is_y(c3)) return 'F';
			if (is_r(c3)) return 'L';
		} else if (c2 == 'C') { // UCx
			return 'S';
		} else if (c2 == 'A') { // UAx
			if (is_y(c3)) return 'Y';
			if (is_r(c3)) return '*';
		} else if (c2 == 'G') { // UGx
			if (is_y(c3)) return 'C';
			if (c3 == 'A') return '*';
			if (c3 == 'G') return 'W';
		}
	} else if (c1 == 'C') {
		///////
		// Cxx
		///////
		if (c2 == 'U') { // CUx
			return 'L';
		} else if (c2 == 'C') { // CCx
			return 'P';
		} else if (c2 == 'A') { // CAx
			if (is_y(c3)) return 'H';
			if (is_r(c3)) return 'Q';
		} else if (c2 == 'G') { // CGx
			return 'R';
		}
	} else if (c1 == 'A') {
		///////
		// Axx
		///////
		if (c2 == 'U') { // AUx
			if (is_y(c3) || c3 == 'A') return 'I';
			if (c3 == 'G') return 'M';
		} else if (c2 == 'C') { // ACx
			return 'T';
		} else if (c2 == 'A') { // AAx
			if (is_y(c3)) return 'N';
			if (is_r(c3)) return 'K';
		} else if (c2 == 'G') { // AGx
			if (is_y(c3)) return 'S';
			if (is_r(c3)) return 'R';
		}
	} else if (c1 == 'G') {
		///////
		// Gxx
		///////
		if (c2 == 'U') { // GUx
			return 'V';
		} else if (c2 == 'C') { // GCx
			return 'A';
		} else if (c2 == 'A') { // GAx
			if (is_y(c3)) return 'D';
			if (is_r(c3)) return 'E';
		} else if (c2 == 'G') { // GGx
			return 'G';
		}
	}

	return 0;
}

// is_stop_codon:
bool GeneticCode::is_stop_codon(char c1, char c2, char c3)
{
	c1 = purify(c1);
	c2 = purify(c2);
	c3 = purify(c3);
	if (c1 == 'U') {
		if (c2 == 'A' && is_r(c3)) return true;
		if (c2 == 'G' && c3 == 'A') return true;
	}
	return false;
}

// purity:
char GeneticCode::purify(char c)
{
	if ((c = std::toupper(c)) == 'T') return 'U';
	return c;
}
