#ifndef CBA_SECONDARY_STRUCTURE_H
#define CBA_SECONDARY_STRUCTURE_H

#include <vector>
#include <list>
#include <functional>
#include <iostream>
#include <utility>
#include "protein.h"
#include "position.h"

namespace Cba {

class SecondaryStructure {
public:
	SecondaryStructure(std::vector<Protein*>& proteins);
	void print(std::ostream& os) const;
	void print_seq(std::ostream& os) const; // secondary structure sequence

	int nhelices() const;
	int nstrands() const;
	int nsheets() const;

	struct AlphaHelix {
		AlphaHelix() : from(0), to(0) {}
		AlphaHelix(int i, int j) : from(i), to(j) {}
		int from; // residue number
		int to; // residue number
	};

	struct BetaStrand {
		BetaStrand() : from(0), to(0), sense(true) {}
		BetaStrand(int i, int j, bool p = true) : from(i), to(j), sense(p) {}
		int from; // residue number
		int to; // residue number
		bool sense;
			// true for parallel, false for anti-parallel
			// against first strand of sheet
	};

	struct BetaSheet {
		std::vector<BetaStrand> strands;
		void add(int i0, int i1, bool parallel = true);
	};

	int nhelices_of(const Protein* prot) const;
	int nstrands_of(const Protein* prot) const;
	AlphaHelix get_helix(size_t i, const Protein* prot) const;
	BetaStrand get_strand(size_t i, const Protein* prot) const;
	const BetaSheet& get_sheet(size_t i) const;
	std::pair<Protein*, int> original_residue_number_of(int i) const;

private:
	std::vector<AlphaHelix> m_alpha_helices;
	std::vector<BetaSheet> m_beta_sheets;
	std::vector<Protein*> m_proteins;

	struct Residue {
		Residue(Protein::Residue& residue);
		~Residue();
		Position* n;
		Position* h;
		Position* ca;
		Position* c;
		Position* o;
		Position* my_h; // created if h not available
	};
	std::vector<Residue> m_residues; // positions of residues of all proteins

	struct Bridge {
		int r0;
		int r1;
		bool parallel;
	};

	struct Ladder {
		Ladder() {}
		Ladder(int i0, int i1, int j0, int j1, bool p)
			: r0(i0), r1(i1), s0(j0), s1(j1), parallel(p) {}
		int r0, r1; // start and end of a strand of ladder
		int s0, s1; // start and end of another strand of ladder
		bool parallel;
	};

	struct BetaStrandCmp
		: std::binary_function<BetaStrand, BetaStrand, bool> {
		bool operator()(const BetaStrand& s1, const BetaStrand& s2)
		{ return s1.from < s2.from; }
	};

	static const double CHAIN_BREAK = 2.5;
	static const double MIN_DIST = 0.5;
	static const double HBOND_HIGH = -0.5;
	static const double HBOND_LOW = -9.9;
	static const double Q = 27.888;

	void set_residues();
	void make_hs(int r0, int n);
	bool chain_break(int i) const;
	bool has_hbond(size_t a, size_t d) const;
	void define_alpha_helices();
	void define_beta_sheets();
	void find_turns(std::vector<char>& turns);
	void find_bridges(std::vector<Bridge>& bridges);
	void find_ladders(std::vector<Ladder>& ladders);
	void find_sheet(std::vector<Ladder>& ladders);
	int first_ladder_to_sheet(
		std::vector<Ladder>& ladders, std::deque<bool>& in_sheet) const;
	bool connected_ladders(Ladder& lad1, Ladder& lad2) const;
	void make_sheet(std::vector<Ladder>& ladders);

	void get_secondary_structure_sequence_of(
		const Protein* prot, std::string& s) const;
	void set_secondary_structures();
};

inline int SecondaryStructure::nhelices() const
{ return m_alpha_helices.size(); }
inline int SecondaryStructure::nsheets() const
{ return m_beta_sheets.size(); }
inline const SecondaryStructure::BetaSheet&
SecondaryStructure::get_sheet(size_t i) const
{ return m_beta_sheets.at(i); }

}
#endif
