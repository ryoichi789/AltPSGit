#ifndef CBA_LOCAREA_H
#define CBA_LOCAREA_H

#include <functional>
#include <cstring>
#include <sstream>
#include "app.h"
#include "molecule.h"
#include "protein.h"
#include "superposer.h"


namespace Cba {

struct LocArea{
public:
	LocArea();
	LocArea(int i);
	int c_atom;// serial number of centrul atom
	bool dense;// if number of s_atoms >= 5, dense is true
	std::vector<int> s_atoms;//surface atoms;
	std::vector<int> resn_of_satoms;//residue number (serial) containing surface atom;
	std::vector<int> i_atoms;//inside atoms;
	Position c_pos;//position of c_atom
	Position g_pos;//center of gravity
	Position g_vec;// vector: g_pos to center atom
	std::vector<Position> satom_positions;
	std::vector<std::vector<double> > pp_vec;//property vector
	Superposer spp;
	//-func-
	std::vector<Position> ori_satom_positions(const Molecule& mol);
	std::vector<Position> ori_iatom_positions(const Molecule& mol);
	void calc_patty_vector(const Molecule& mol);//Calculete Propaty-vector at c_atom with PATTY
	int nsatoms() const;
	int niatoms() const;
	double tc(const LocArea& la,int i);
	double tc1(const LocArea& l);
	double tc3(const LocArea& l);
	void output_pdb(const Protein& prot);
	void output_pdb(const Protein& prot, bool n_axis);
	void output_pdb(const Protein& prot, std::string fname, bool n_axis);
	void output_superpose_pdb(const Protein& prot);
	void output_superpose_pdb(const Protein& prot, bool n_axis);
	void output_superpose_pdb(const Protein& prot, std::string fname, bool n_axis);
	LocArea mix(const LocArea& l) const;
	int get_atom_number(int n);
};
inline int LocArea::nsatoms() const{return s_atoms.size();}
inline int LocArea::niatoms() const{return i_atoms.size();}

void set_local_area(const Protein& prot, std::vector<LocArea>& lareas, 
	std::vector<int>& suratoms, std::vector<std::vector<double> >& dmat);
void assign_property_larea(const Protein& prot, std::vector<LocArea>& lareas, 
	const std::vector<int>& suratoms);
std::string chain_id(const Protein& prot);
void print_atom_in_pdb(std::ostream& os, std::string chain, const Protein::Residue& res, const Atom* atom);
void output_pdb(const Protein& prot, std::string fname);
void output_superpose_pdb(const Protein& prot, std::string fname, const Superposer& spp);
std::vector<int> get_serial_number_of_residue(const Protein& prot);
std::pair<int, double> calc_patty_comp(const int i,const Atom* ca, const Atom* ta, const double d0);
std::string i_to_s(int i);
inline std::string i_to_s(int i){std::stringstream s;s << i;return s.str();}
bool patty_match(int a, int b);
double angle(const Position& a, const Position& b);
std::vector<int> get_vector_top(int num, const std::vector<int>& list, const std::vector<double>& score);
std::string patty_atname(int i);
}
#endif
