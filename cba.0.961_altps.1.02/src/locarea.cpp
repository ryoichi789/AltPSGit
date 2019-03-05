#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <ios>
#include <iomanip>
#include "locarea.h"
#include "sdf_entry.h"
#include "sdf_reader.h"
#include "pdb_entry.h"
#include "pdb_reader.h"
#include "molecule.h"

#include <dirent.h>

using namespace Cba;

//get_serial_number_of_residue
std::vector<int> Cba::get_serial_number_of_residue(const Protein& prot)
{
	std::vector<int> number;
	for(int i=0;i<prot.natoms();i++){
		const Atom* a = prot.atom(i);
		number.push_back(prot.residue_containing_atom(a));
	}
	return number;
}

//chain_id
std::string Cba::chain_id(const Protein& prot)
{
	std::string mol_id = prot.id();
	std::string chain;
	if(mol_id.find_last_of('_') != std::string::npos){
		int n = mol_id.find_last_of('_');
		chain = mol_id.substr(n-1,1);
	}else{
		chain = " ";
	}
	return chain;
}

//print_atom_in_pdb
void Cba::print_atom_in_pdb(std::ostream& os, std::string c, const Protein::Residue& r, const Atom* a)
{
	os << "ATOM  "
	   << std::setw(5) << std::right << a->number
	   << " " << a->name
	   << " " << r.name << " " << c
	   << std::setw(5) << std::right << r.number << "   "
	   << std::fixed << std::right << std::setprecision(3) << std::setw(8) << a->pos.x()
	   << std::fixed << std::right << std::setprecision(3) << std::setw(8) << a->pos.y()
	   << std::fixed << std::right << std::setprecision(3) << std::setw(8) << a->pos.z()
	   << "  1.00"
	   << std::fixed << std::right << std::setprecision(2) << std::setw(6) << a->asa
	   << "           "
	   << Atom::symbols[a->element] << "  "
	   << std::endl;

}

//output_pdb
void Cba::output_pdb(const Protein& prot, std::string fname)
{
	std::ofstream os(fname.c_str());
	for(int i=0;i<prot.natoms();i++){
		const Atom* atom = prot.atom(i);
		int resnum = prot.residue_containing_atom(atom);
		const Protein::Residue res = prot.residue(resnum);
		std::string chain = chain_id(prot);
		print_atom_in_pdb(os,chain,res,atom);
	}
}

//output_superpose_pdb
void Cba::output_superpose_pdb(const Protein& prot, std::string fname, const Superposer& spp)
{
	Protein prot2;
	prot.clone(&prot2);
	for(int i=0;i<prot.natoms();i++){
		Atom* atom = prot2.atom(i);
		spp.transform(atom->pos);
	}
	output_pdb(prot2, fname);
}

//calc_patty_comp
std::pair<int, double> Cba::calc_patty_comp(const int i,const Atom* ca, const Atom* ta, const double d0)
{
	std::pair<int, double> pv;
	double dlen = std::fabs(ta->pos.distance_from(ca->pos) - d0 * i);
	if(dlen < d0){
		if(ta->pc_class <= 4){
			pv.first = ta->pc_class - 1;//PATTY :physico-chemical class of Bush & Sheridan (1993)
		}else if(ta->pc_class == 5){
			pv.first = 6;
		}else{
			pv.first = ta->pc_class - 2;
		}
		pv.second = 1.0 - dlen/d0;
	}else{
		pv.first = -1;
		pv.second = 0.0;
	}
	return pv;
}

//set_local_area
void Cba::set_local_area(const Protein& prot, std::vector<LocArea>& lareas, std::vector<int>& suratoms,
	std::vector<std::vector<double> >& dmat)
{
	std::vector<int> resnumbers = get_serial_number_of_residue(prot);
	std::vector<int> inatoms;
	double d;
	// int i : serial number of atom
	for (int i = 0; i < prot.natoms(); i++) {
		LocArea larea(i);
		lareas.push_back(larea);
		if(prot.atom(i)->asa > 0.0 ){
			suratoms.push_back(i);
			lareas[i].s_atoms.push_back(i);
			lareas[i].resn_of_satoms.push_back(resnumbers[i]);
		}else{inatoms.push_back(i);}
	}
	for (int i = 0; i < suratoms.size()-1; i++) {
		int ii = suratoms[i];
		for (int j = i+1; j < suratoms.size(); j++) {
			int jj = suratoms[j];
			d = dmat[ii][jj];
			if(d <= 5.0){
				lareas[ii].s_atoms.push_back(jj);lareas[ii].resn_of_satoms.push_back(resnumbers[jj]);
				lareas[jj].s_atoms.push_back(ii);lareas[jj].resn_of_satoms.push_back(resnumbers[ii]);
			}
		}
	}
	for (int i = 0; i < suratoms.size(); i++) {
		int ii = suratoms[i];
		for (int j = 0; j < inatoms.size(); j++) {
			int jj = inatoms[j];
			d = prot.atom(ii)->pos.distance_from_neighbor(prot.atom(jj)->pos,5.0);
			dmat[ii][jj] = d;dmat[jj][ii] = d;
			if(d >= 0.0){
				lareas[ii].i_atoms.push_back(jj);
			}
		}
		//Calculate center position of gravity at local-area[ii]
		lareas[ii].satom_positions = lareas[ii].ori_satom_positions(prot);
		std::vector<Position> iatom_pos = lareas[ii].ori_iatom_positions(prot);
		std::vector<double> g_xyz(3,0.0);
		for(int i = 0; i <lareas[ii].satom_positions.size(); i++)
		{
			g_xyz[0] += lareas[ii].satom_positions[i].x();
			g_xyz[1] += lareas[ii].satom_positions[i].y();
			g_xyz[2] += lareas[ii].satom_positions[i].z();
		}
		for(int i = 0; i <iatom_pos.size(); i++)
		{
			g_xyz[0] += iatom_pos[i].x();g_xyz[1] += iatom_pos[i].y();g_xyz[2] += iatom_pos[i].z();
		}
		double num = lareas[ii].satom_positions.size() + iatom_pos.size();
		Position p(g_xyz[0]/num,g_xyz[1]/num,g_xyz[2]/num);
		lareas[ii].g_pos = p;
		lareas[ii].c_pos = lareas[ii].satom_positions[0];
		//
		Position nvec(lareas[ii].satom_positions[0].x() - p.x(),
		              lareas[ii].satom_positions[0].y() - p.y(),
					  lareas[ii].satom_positions[0].z() - p.z());
		lareas[ii].g_vec = nvec.unit();
	}
	for (int i = 0; i < suratoms.size(); i++) {
		int ii = suratoms[i];
		if(lareas[ii].nsatoms() >= 3) lareas[ii].dense = true;
	}
}

//assign_property_larea
void Cba::assign_property_larea(const Protein& prot, std::vector<LocArea>& lareas, 
	const std::vector<int>& suratoms)
{
	for(int i=0;i<suratoms.size();i++){
		lareas[suratoms[i]].calc_patty_vector(prot);
	}
}

//angle
double Cba::angle(const Position& a, const Position& b)
{
	double size_a = a.size();
	if(size_a == 0.0) return 0.0;
	double size_b = b.size();
	if(size_b == 0.0) return 0.0;
	return (a.x()*b.x() + a.y()*b.y() + a.z()*b.z()) / (size_a * size_b);
}

//patty_match
bool Cba::patty_match(int a, int b)
{
	if(a == 0 || b == 0){return false;}
	else if(a == b){return true;}
	else if((a == 3 || a == 4) && b == 5){return true;}
	else if((b == 3 || b == 4) && a == 5){return true;}
	else{return false;}
}

//get_vector_top
std::vector<int> Cba::get_vector_top(int num, const std::vector<int>& l, const std::vector<double>& s)
{
	std::vector<int> list;
	std::vector<double> score = s;
	std::sort(score.begin(),score.end(),std::greater<double>());
	for(int i = 0; i < l.size() ; i++){
		if(s[i] >= score[num]) list.push_back(l[i]);
	}
	return list;
}

//patty_atname
std::string Cba::patty_atname(int c)
{
	std::string at;
	if (c == Atom::CATION) at =           "CATION               ";
	else if (c == Atom::ANION) at =       "ANION                ";
	else if (c == Atom::DONOR) at =       "DONOR                ";
	else if (c == Atom::ACCEPTOR) at =    "ACCEPTOR             ";
	else if (c == Atom::POLAR) at =       "POLAR(DONOR/ACCEPTOR)";
	else if (c == Atom::HYDROPHOBIC) at = "HYDROPHOBIC          ";
	else if (c == Atom::NONE) at =        "NONE                 ";
	else at =                             "UNDEFINED            ";
	return at;
}

// LocArea
LocArea::LocArea()
{
	dense = false;
}

LocArea::LocArea(int i){
	c_atom = i;
	dense = false;
}

//ori_satom_positions
std::vector<Position>  LocArea::ori_satom_positions(const Molecule& mol)
{
	std::vector<Position> atom_pos;
	std::vector<int>::iterator it = s_atoms.begin();
	while( it != s_atoms.end() )
	{
		atom_pos.push_back(mol.position_of_atom(*it));
		++it;
	}
	return atom_pos;
}

//ori_iatom_positions
std::vector<Position>  LocArea::ori_iatom_positions(const Molecule& mol)
{
	std::vector<Position> atom_pos;
	std::vector<int>::iterator it = i_atoms.begin();
	while( it != i_atoms.end() )
	{
		atom_pos.push_back(mol.position_of_atom(*it));
		++it;
	}
	return atom_pos;
}

//calc_patty_vector
void LocArea::calc_patty_vector(const Molecule& mol)
{
	std::pair<int, double> pv;
	std::vector<double> vec(6,0.0);
	std::vector<std::vector<double> > pvec(3,vec);
	const Atom* ca = mol.atom(s_atoms[0]); 
	for(int i=0;i<nsatoms();i++){
		const Atom* ta = mol.atom(s_atoms[i]); 
		for(int j = 0;j < 3;j++){
			pv = calc_patty_comp(j, ca, ta, 3.2);
			if(pv.first >= 0){
				if(pv.first == 6){
					pvec[j][2] += pv.second/2.0;
					pvec[j][3] += pv.second/2.0;
				}else{
					pvec[j][pv.first] += pv.second;
				}
			}
		}
	}
	pp_vec = pvec;
}

//output_pdb
void LocArea::output_pdb(const Protein& prot)
{
	output_pdb(prot, false);
}

void LocArea::output_pdb(const Protein& prot, bool n_axis)
{
	std::string outfile("./");
	outfile = outfile + prot.id() + "_la" + i_to_s(c_atom) + ".pdb";
	output_pdb(prot, outfile, n_axis);
}

void LocArea::output_pdb(const Protein& prot, std::string fname, bool n_axis)
{
	std::ofstream os(fname.c_str());

	for(int i=0;i<nsatoms();i++){
		const Atom* atom = prot.atom(s_atoms[i]);
		const Protein::Residue res = prot.residue(resn_of_satoms[i]);
		std::string chain = chain_id(prot);
		print_atom_in_pdb(os,chain,res,atom);
	}
	if(n_axis == true){
		double x,y,z;
		for(int i = 0;i<10;i++){
			x = g_vec.x() * i + g_pos.x();y = g_vec.y() * i + g_pos.y();z = g_vec.z() * i + g_pos.z();
			os << "ATOM  ";
			os << std::setw(5) << std::right << i;
			os << "  H   NON          ";
			os << std::fixed << std::right << std::setprecision(3) << std::setw(8) << x;
			os << std::fixed << std::right << std::setprecision(3) << std::setw(8) << y;
			os << std::fixed << std::right << std::setprecision(3) << std::setw(8) << z;
			os << "  1.00                 H  ";
			os << std::endl;
		}
	}
}

//output_superpose_pdb
void LocArea::output_superpose_pdb(const Protein& prot)
{
	output_superpose_pdb(prot, false);
}

void LocArea::output_superpose_pdb(const Protein& prot, bool n_axis)
{
	std::string outfile("./");
	outfile = outfile + prot.id() + "_la" + i_to_s(c_atom) + ".pdb";
	output_superpose_pdb(prot, outfile, n_axis);
}

void LocArea::output_superpose_pdb(const Protein& prot, std::string fname, bool n_axis)
{
	Protein prot2;
	prot.clone(&prot2);
	for(int i=0;i<nsatoms();i++){
		Atom* atom = prot2.atom(s_atoms[i]);
		spp.transform(atom->pos);
	}

	LocArea loc = *this;
	spp.transform(loc.g_pos);

	Position center = prot2.position_of_atom(loc.s_atoms[0]);
	loc.g_vec = center - loc.g_pos;
	loc.output_pdb(prot2, fname, n_axis);
}

//tc:tanimoto_coefficient
double LocArea::tc(const LocArea& la, int i)
{
	double s = -1.0;
	double aa = 0.0;double bb = 0.0;double ab = 0.0;
	for(int k=0;k<pp_vec[i].size();k++){
		aa += pp_vec[i][k] * pp_vec[i][k];
		bb += la.pp_vec[i][k] * la.pp_vec[i][k];
		ab += pp_vec[i][k] * la.pp_vec[i][k];
	}
	if (aa + bb - ab > 0) s = ab / (aa + bb - ab);
	return s;
}

//tc1:tanimoto_coefficient
double LocArea::tc1(const LocArea& l)
{
	double s = -1.0;
	s = this->tc(l,0);
	return s;
}

//tc3:tanimoto_coefficient
double LocArea::tc3(const LocArea& l)
{
	double s = -1.0;
	s = (this->tc(l,0) + this->tc(l,1) + this->tc(l,2) ) /3.0;
	return s;
}

//mix
LocArea LocArea::mix(const LocArea& l) const
{
	LocArea area;
	area = *this;
	area.dense = true;
	//s_atoms, resn_of_satoms, satom_positions
	std::map<int,int> chk;
	for(int i = 0; i < this->s_atoms.size(); i++){
		chk.insert(std::pair<int,int>(this->s_atoms[i],1));
	}
	for(int i = 0; i < l.s_atoms.size(); i++){
		if(chk[l.s_atoms[i]] != 1){
			area.s_atoms.push_back(l.s_atoms[i]);
			area.resn_of_satoms.push_back(l.resn_of_satoms[i]);
			area.satom_positions.push_back(l.satom_positions[i]);
		}
	}
	//c_pos
	Position p1((this->c_pos.x() + l.c_pos.x())/2.0, (this->c_pos.y() + l.c_pos.y())/2.0, (this->c_pos.z() + l.c_pos.z())/2.0);
	area.c_pos = p1;
	//g_pos
	double m1 = this->s_atoms.size();
	double m2 = l.s_atoms.size();
	Position p2( (m1*this->g_pos.x() + m2*l.g_pos.x())/(m1 + m2),
				 (m1*this->g_pos.y() + m2*l.g_pos.y())/(m1 + m2),
				 (m1*this->g_pos.z() + m2*l.g_pos.z())/(m1 + m2) );
	area.g_pos = p2;
	//g_vec
	Position p3( p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
	area.g_vec = p3.unit();
	return area;
} 

//get_atom_number
int LocArea::get_atom_number(int n)
{
	for(int i = 0; i < s_atoms.size(); i++)
	{
		if(s_atoms[i] == n) return i;
	}
	return -1;
}
