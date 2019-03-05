#ifndef CBA_GENETIC_CODE
#define CBA_GENETIC_CODE

namespace Cba {

class GeneticCode {
public:
	static char decode(char c1, char c2, char c3 = 'N');
		// returns amino acid symbol for the codon ('*' for stop codon)
	static bool is_stop_codon(char c1, char c2, char c3);

private:
	static char purify(char c); // converts c to upper-case, and T to U
	static bool is_r(char c); // true if purine
	static bool is_y(char c); // true if pyrimidine
};

inline bool GeneticCode::is_r(char c) { return c=='A'||c=='G'||c=='R'; }
inline bool GeneticCode::is_y(char c) { return c=='U'||c=='C'||c=='Y'; }
}
#endif
