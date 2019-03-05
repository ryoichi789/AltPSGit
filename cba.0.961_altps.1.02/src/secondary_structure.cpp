#include <vector>
#include <deque>
#include <list>
#include <algorithm>
#include <iostream>
#include <string>
#include "secondary_structure.h"
#include "protein.h"
#include "position.h"
#include "bio_sequence.h"

using namespace Cba;
using std::vector;
using std::list;

//
SecondaryStructure::SecondaryStructure(vector<Protein*>& proteins)
	: m_proteins(proteins)
{
	set_residues();
	define_alpha_helices();
	define_beta_sheets();
	set_secondary_structures();
}

// set_residues:
void SecondaryStructure::set_residues()
{
	for (vector<Protein*>::iterator
		p = m_proteins.begin(); p != m_proteins.end(); ++p)
		for (int i = 0; i < (*p)->nresidues(); i++) {
			m_residues.push_back(Residue((*p)->residue(i)));
		}

	// set amide hydrogens for each protein
	int n = 0;
	for (vector<Protein*>::iterator
		p = m_proteins.begin(); p != m_proteins.end(); ++p) {
		make_hs(n, (*p)->nresidues());
		n += (*p)->nresidues();
	}
}

// make_hs: for n residues from r0
void SecondaryStructure::make_hs(int r0, int n)
{
	// for N-terminal residue
	if (m_residues[r0].h == 0 && m_residues[r0].n != 0) {
		m_residues[r0].my_h = new Position;
		m_residues[r0].h = m_residues[r0].my_h;
		*(m_residues[r0].my_h) = *(m_residues[r0].n);
	}

	Position p;
	for (int i = r0+1; i < r0+n; i++) {
		if (m_residues[i].h || m_residues[i].n == 0) continue;
		m_residues[i].my_h = new Position;
		m_residues[i].h = m_residues[i].my_h;
		if (chain_break(i) ||
			m_residues[i-1].c == 0 || m_residues[i-1].o == 0) {
			*(m_residues[i].my_h) = *(m_residues[i].n);
		} else {
			p = *(m_residues[i-1].c) - *(m_residues[i-1].o);
			*(m_residues[i].my_h) = *(m_residues[i].n) + p.unit();
		}
	}
}

// chain_break: returns true if chain breaks before residue i
bool SecondaryStructure::chain_break(int i) const
{
	if (i <= 0 || i >= m_residues.size()) return true;
	if (m_residues[i-1].c == 0 || m_residues[i].n == 0) return true;
	if (m_residues[i].n->distance_from(*(m_residues[i-1].c)) > CHAIN_BREAK)
		return true;
	return false;
}

// has_hbond:
bool SecondaryStructure::has_hbond(size_t a, size_t d) const
{
	if (m_residues[a].c == 0 || m_residues[a].o == 0 ||
		m_residues[d].n == 0 || m_residues[d].h == 0) return false;

	double d_on = m_residues[a].o->distance_from(*(m_residues[d].n));
	double d_oh = m_residues[a].o->distance_from(*(m_residues[d].h));
	double d_cn = m_residues[a].c->distance_from(*(m_residues[d].n));
	double d_ch = m_residues[a].c->distance_from(*(m_residues[d].h));

	if (d_on < MIN_DIST || d_oh < MIN_DIST ||
		d_cn < MIN_DIST || d_ch < MIN_DIST) return true;
	if ((Q/d_on + Q/d_ch - Q/d_oh - Q/d_cn) < HBOND_HIGH) return true;

	return false;
}

// define_alpha_helices:
void SecondaryStructure::define_alpha_helices()
{
	vector<char> turns(m_residues.size(), ' ');
	find_turns(turns);

	vector<char> alpha_helices(m_residues.size(), ' ');
	for (int i = 1; i < m_residues.size()-4; i++) {
		if((turns[i-1] == '>' || turns[i-1] == 'X') &&
			(turns[i] == '>' || turns[i] == 'X')) {
			for (int j = i; j < i+4; j++) alpha_helices[j] = 'H';
		}
	}

	for (int i = 0; i < m_residues.size()-4; ) {
		if (alpha_helices[i] == 'H') {
			int j = i + 1;
			while (j < m_residues.size() && alpha_helices[j] == 'H')
				if (chain_break(j)) break;
				else j++;
			m_alpha_helices.push_back(AlphaHelix(i,j-1));
			i = j;
		} else {
			i++;
		}
	}
}

// find_turns:
void SecondaryStructure::find_turns(vector<char>& turns)
{
	for (int i = 0; i < m_residues.size()-4; ) {
		if (chain_break(i+1) || chain_break(i+2) ||
			chain_break(i+3) || chain_break(i+4)) {
			i += 4;
		} else {
			if (has_hbond(i,i+4)) {
				if (turns[i] == '<') turns[i] = 'X';
				else turns[i] = '>';
				turns[i+4] = '<';
			}
			i++;
		}
	}
}

// define_beta_sheets:
void SecondaryStructure::define_beta_sheets()
{
	vector<Ladder> ladders;
	find_ladders(ladders);
	find_sheet(ladders);
}

// find_ladders:
void SecondaryStructure::find_ladders(vector<Ladder>& ladders)
{
	vector<Bridge> bridges;
	find_bridges(bridges);

	std::deque<bool> checked(bridges.size(), false);
	int i1, i2, j1, j2;
	for(int i = 0; i < bridges.size(); i++) {
		if (checked[i]) continue;
		i1 = i2 = bridges[i].r0;
		j1 = j2 = bridges[i].r1;
		checked[i] = true;

		for (int j = i+1; j < bridges.size(); j++) {
			if (checked[j]) continue;
			if (bridges[i].parallel && (! bridges[j].parallel)) continue;
			if ((! bridges[i].parallel) && bridges[j].parallel) continue;
			if(i2 == bridges[j].r0) continue;

			if(bridges[i].parallel) { // parallel
				if ((i2+1) == bridges[j].r0 && (j2+1) == bridges[j].r1) {
					i2 = bridges[j].r0;
					j2 = bridges[j].r1;
					checked[j] = true;
				}
			} else { // anti-parallel
				if((i2+1) == bridges[j].r0 && (j1-1) == bridges[j].r1) {
					i2 = bridges[j].r0;
					j1 = bridges[j].r1;
					checked[j] = true;
				}
			}
		}

		if (i1 < i2)
			ladders.push_back(Ladder(i1,i2,j1,j2,bridges[i].parallel));
	}
}

// find_bridges:
void SecondaryStructure::find_bridges(std::vector<Bridge>& bridges)
{
	Bridge b;
	for(int i = 1; i < m_residues.size()-4; i++) {
		if(chain_break(i) || chain_break(i+1)) continue;
		for (int j = i+3; j < m_residues.size()-1; j++) {
			if(chain_break(j) || chain_break(j+1)) continue;

			// parallel bridge
			if ((has_hbond(i-1,j) && has_hbond(j,i+1)) ||
				(has_hbond(j-1,i) && has_hbond(i,j+1))) {
				b.r0 = i;
				b.r1 = j;
				b.parallel = true;
				bridges.push_back(b);
			}

			// antiparallel bridge
			if ((has_hbond(i,j) && has_hbond(j,i)) ||
				(has_hbond(i-1,j+1) && has_hbond(j-1,i+1))) {
				b.r0 = i;
				b.r1 = j;
				b.parallel = false;
				bridges.push_back(b);
			}
		}
	}
}

// find_sheet:
void SecondaryStructure::find_sheet(vector<Ladder>& ladders)
{
	if (ladders.empty()) return;
	std::deque<bool> in_sheet(ladders.size(), false);
	vector<Ladder> ladders_in_sheet;
	int i0, i;
	while ((i0 = first_ladder_to_sheet(ladders, in_sheet)) < ladders.size()) {
		// initially, a new sheet contains only first 'free' ladder
		ladders_in_sheet.clear();
		ladders_in_sheet.push_back(ladders[i0]);
		in_sheet[i0] = true;

		while (i0 < ladders.size()) {
			for (i = 0; i < ladders.size(); i++) {
				if (in_sheet[i]) continue;
				if (connected_ladders(ladders[i0], ladders[i])) {
					ladders_in_sheet.push_back(ladders[i]);
					in_sheet[i] = true;
					break;
				}
			}
			i0 = i;
		}
		make_sheet(ladders_in_sheet);
	}
}

// connect_ladders:
bool SecondaryStructure::connected_ladders(Ladder& lad1, Ladder& lad2) const
{
	// compare first strand of lad1 against first strand of lad2
	if (lad1.r0 <= lad2.r1 && lad2.r0 <= lad1.r1) return true;
	// compare first strand of lad1 against second strand of lad2
	if (lad1.r0 <= lad2.s1 && lad2.s0 <= lad1.r1) return true;
	// compare second strand of lad1 against first strand of lad2
	if (lad1.s0 <= lad2.r1 && lad2.r0 <= lad1.s1) return true;
	// compare second strand of lad1 against second strand of lad2
	if (lad1.s0 <= lad2.s1 && lad2.s0 <= lad1.s1) return true;
	return false;
}

// make_sheet:
void SecondaryStructure::make_sheet(vector<Ladder>& ladders)
{
	if (ladders.empty()) return;
	BetaSheet sheet;

	if (ladders.size() == 1) {
		sheet.add(ladders[0].r0, ladders[0].r1);
		sheet.add(ladders[0].s0, ladders[0].s1, ladders[0].parallel);
		m_beta_sheets.push_back(sheet);
		return;
	}

	// start sheet from strand at edge
	int i0, i1;
	if (ladders[0].r0 <= ladders[1].r1 && ladders[1].r0 <= ladders[0].r1) {
		sheet.add(ladders[0].s0, ladders[0].s1);
		sheet.add(std::min(ladders[0].r0, ladders[1].r0),
			std::max(ladders[0].r1, ladders[1].r1), ladders[0].parallel);
		sheet.add(ladders[1].s0, ladders[1].s1, ladders[1].parallel);
	} else if (ladders[0].r0 <= ladders[1].s1 && ladders[1].s0 <= ladders[0].r1) {
		sheet.add(ladders[0].s0, ladders[0].s1);
		sheet.add(std::min(ladders[0].r0, ladders[1].s0),
			std::max(ladders[0].r1, ladders[1].s1), ladders[0].parallel);
		sheet.add(ladders[1].r0, ladders[1].r1, ladders[1].parallel);
	} else if (ladders[0].s0 <= ladders[1].r1 && ladders[1].r0 <= ladders[0].s1) {
		sheet.add(ladders[0].r0, ladders[0].r1);
		sheet.add(std::min(ladders[0].s0, ladders[1].r0),
			std::max(ladders[0].s1, ladders[1].r1), ladders[0].parallel);
		sheet.add(ladders[1].s0, ladders[1].s1, ladders[1].parallel);
	} else {
		sheet.add(ladders[0].r0, ladders[0].r1);
		sheet.add(std::min(ladders[0].s0, ladders[1].s0),
			std::max(ladders[0].s1, ladders[1].s1), ladders[0].parallel);
		sheet.add(ladders[1].r0, ladders[1].r1, ladders[1].parallel);
	}
	// add strands to sheet	
	for (int i = 2; i < ladders.size(); i++) {
		BetaStrand& strand = sheet.strands.back();
		if (ladders[i].r0 <= strand.to && strand.from <= ladders[i].r1) {
			strand.from = std::min(strand.from, ladders[i].r0);
			strand.to = std::max(strand.to, ladders[i].r1);
			sheet.add(ladders[i].s0, ladders[i].s1, ladders[i].parallel);
		} else if (ladders[i].s0 <= strand.to && strand.from <= ladders[i].s1) {
			strand.from = std::min(strand.from, ladders[i].s0);
			strand.to = std::max(strand.to, ladders[i].s1);
			sheet.add(ladders[i].r0, ladders[i].r1, ladders[i].parallel);
		}
	}
	m_beta_sheets.push_back(sheet);
}

// first_ladder_to_sheet:
int SecondaryStructure::first_ladder_to_sheet(
	vector<Ladder>& ladders, std::deque<bool>& in_sheet) const
{
	// find ladder at edge
	int i0 = ladders.size();
	int nc;
	for (int i = 0; i < ladders.size(); i++) {
		if (in_sheet[i]) continue;
		nc = 0;
		for (int j = 0; j < ladders.size(); j++) {
			if (i == j) continue;
			if (in_sheet[j]) continue;
			if (connected_ladders(ladders[i], ladders[j])) nc++;
		}
		if (nc < 2) {
			i0 = i;
			break;
		}
	}
	if (nc == 2) { // for barrel
		for (int i = 0; i < ladders.size(); i++) {
			if (! in_sheet[i]) {
				i0 = i;
				break;
			}
		}
	}
	return i0;
}

// original_residue_of:
std::pair<Protein*,int>
SecondaryStructure::original_residue_number_of(int i) const
{
	int n = 0;
	for (vector<Protein*>::const_iterator
		p = m_proteins.begin(); p != m_proteins.end(); ++p) {
		i -= n;
		if (i < (n = (*p)->nresidues()))
			return std::make_pair<Protein*,int>(*p,i);
	}
	return std::make_pair<Protein*,int>(0,0);
}

// print:
void SecondaryStructure::print(std::ostream& os) const
{
	int r0, r1;
	std::pair<Protein*,int> resnum;

	// helices
	for (vector<AlphaHelix>::const_iterator
		p = m_alpha_helices.begin(); p != m_alpha_helices.end(); ++p) {
		resnum = original_residue_number_of(p->from);
		r0 = resnum.second;
		resnum = original_residue_number_of(p->to);
		r1 = resnum.second;
		os << "A " << resnum.first->id() << '\t'
			<< r0+1 << '\t' << r1+1 << std::endl;
	}
	os << std::endl;

	// sheets
	for (vector<BetaSheet>::const_iterator
		p = m_beta_sheets.begin(); p != m_beta_sheets.end(); ++p) {
		for (vector<BetaStrand>::const_iterator
			q = p->strands.begin(); q != p->strands.end(); ++q) {
			resnum = original_residue_number_of(q->from);
			r0 = resnum.second;
			resnum = original_residue_number_of(q->to);
			r1 = resnum.second;
			os << "B " << resnum.first->id() << '\t'
				<< r0+1 << '\t' << r1+1 << '\t';
			if (q->sense) os << '1';
			else os << '0';
			os << std::endl;
		}
		os << std::endl;
	}
}

// print_seq:
void SecondaryStructure::print_seq(std::ostream& os) const
{
	BioSequence seq;
	std::string s;
	for (vector<Protein*>::const_iterator
		prot = m_proteins.begin(); prot != m_proteins.end(); ++prot) {
		seq.clear();
		seq.set_id((*prot)->id());
		get_secondary_structure_sequence_of(*prot, s);
		seq.add(s);
		seq.print(os);
	}
}

// get_secondary_structure_sequence_of:
void SecondaryStructure::get_secondary_structure_sequence_of(
	const Protein* prot, std::string& s) const
{
	std::pair<Protein*,int> resnum;
	int r0, r1;

	s.assign(prot->nresidues(), '_');
	for (vector<AlphaHelix>::const_iterator
		p = m_alpha_helices.begin(); p != m_alpha_helices.end(); ++p) {
		resnum = original_residue_number_of(p->from);
		r0 = resnum.second;
		resnum = original_residue_number_of(p->to);
		r1 = resnum.second;
		if (resnum.first == prot)
			for (int i = r0; i <= r1; i++) s[i] = 'a';
	}
	for (vector<BetaSheet>::const_iterator
		p = m_beta_sheets.begin(); p != m_beta_sheets.end(); ++p) {
		for (vector<BetaStrand>::const_iterator
			q = p->strands.begin(); q != p->strands.end(); ++q) {
				resnum = original_residue_number_of(q->from);
				r0 = resnum.second;
				resnum = original_residue_number_of(q->to);
				r1 = resnum.second;
				if (resnum.first == prot)
					for (int i = r0; i <= r1; i++) s[i] = 'b';
		}
	}
}

// nstrands:
int SecondaryStructure::nstrands() const
{
	int n = 0;
	for (vector<BetaSheet>::const_iterator
		p = m_beta_sheets.begin(); p != m_beta_sheets.end(); ++p)
		n += p->strands.size();
	return n;
}

// nhelices_of:
int SecondaryStructure::nhelices_of(const Protein* prot) const
{
	int n = 0;
	std::pair<Protein*,int> resnum;
	for (vector<AlphaHelix>::const_iterator
		p = m_alpha_helices.begin(); p != m_alpha_helices.end(); ++p) {
		resnum = original_residue_number_of(p->from);
		if (resnum.first == prot) n++;
	}
	return n;
}

// nstrands_of:
int SecondaryStructure::nstrands_of(const Protein* prot) const
{
	int n = 0;
	std::pair<Protein*,int> resnum;
	for (vector<BetaSheet>::const_iterator
		p = m_beta_sheets.begin(); p != m_beta_sheets.end(); ++p) {
		for (vector<BetaStrand>::const_iterator
			q = p->strands.begin(); q != p->strands.end(); ++q) {
			resnum = original_residue_number_of(q->from);
			if (resnum.first == prot) n++;
		}
	}
	return n;
}

// get_helix:
SecondaryStructure::AlphaHelix
SecondaryStructure::get_helix(size_t i, const Protein* prot) const
{
	int n = 0;
	int r0, r1;
	std::pair<Protein*,int> resnum;
	for (vector<AlphaHelix>::const_iterator
		p = m_alpha_helices.begin(); p != m_alpha_helices.end(); ++p) {
		resnum = original_residue_number_of(p->from);
		r0 = resnum.second;
		resnum = original_residue_number_of(p->to);
		r1 = resnum.second;
		if (resnum.first == prot) {
			if (i == n) return AlphaHelix(r0, r1);
			n++;
		}
	}
	return AlphaHelix();
}

// get_strand:
SecondaryStructure::BetaStrand
SecondaryStructure::get_strand(size_t i, const Protein* prot) const
{
	int n = 0;
	int r0, r1;
	std::pair<Protein*,int> resnum;
	for (vector<BetaSheet>::const_iterator
		p = m_beta_sheets.begin(); p != m_beta_sheets.end(); ++p) {
		for (vector<BetaStrand>::const_iterator
			q = p->strands.begin(); q != p->strands.end(); ++q) {
			resnum = original_residue_number_of(q->from);
			r0 = resnum.second;
			resnum = original_residue_number_of(q->to);
			r1 = resnum.second;
			if (resnum.first == prot) {
				if (i == n) return BetaStrand(r0, r1);
				n++;
			}
		}
	}
	return BetaStrand();
}

// Residue::Residue:
SecondaryStructure::Residue::Residue(Protein::Residue& residue)
	: n(0), h(0), ca(0), c(0), o(0), my_h(0)
{
	if (residue.n) n = &(residue.n->pos);
	if (residue.h) h = &(residue.h->pos);
	if (residue.ca) ca = &(residue.ca->pos);
	if (residue.c) c = &(residue.c->pos);
	if (residue.o) o = &(residue.o->pos);

	if (residue.name.compare("PRO") == 0) n = 0;
}

// Residue::~Residue:
SecondaryStructure::Residue::~Residue()
{
	if (my_h) delete my_h;
}

// BetaSheet::add:
void SecondaryStructure::BetaSheet::add(int i0, int i1, bool parallel)
{
	if (strands.empty()) {
		strands.push_back(BetaStrand(i0,i1,true));
	} else {
		bool current_sense = strands.back().sense;
		if (parallel) strands.push_back(BetaStrand(i0,i1,current_sense));
		else strands.push_back(BetaStrand(i0,i1,!current_sense));
	}
}

// set_secondary_structures:
void SecondaryStructure::set_secondary_structures()
{
	std::string s;
	for (vector<Protein*>::const_iterator
		prot = m_proteins.begin(); prot != m_proteins.end(); ++prot) {
		get_secondary_structure_sequence_of(*prot, s);
		(*prot)->m_secondary_structures = s;
	}
}
