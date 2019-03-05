#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "protein.h"
#include "molecule.h"
#include "atom.h"
#include "bio_sequence.h"

using namespace Cba;
using std::vector;

//
Protein::Protein() : Molecule() {}

//
Protein::~Protein()
{
	clear();
}

// clear:
void Protein::clear()
{
	Molecule::clear();
	m_residues.clear();
	m_secondary_structures.erase();
	m_lbrs.clear();
}

//
Protein::Protein(const Protein& prot)
{
	prot.clone(this);
}

//
Protein& Protein::operator=(const Protein& prot)
{
	if (this != &prot) {
		clear();
		prot.clone(this);
	}
	return *this;
}

// residue_type:
char Protein::residue_type(size_t i) const
{
	const std::string& s = residue(i).name;
	if (s.compare("ALA") == 0) return 'A';
	if (s.compare("ASX") == 0) return 'B';
	if (s.compare("CYS") == 0) return 'C';
	if (s.compare("ASP") == 0) return 'D';
	if (s.compare("GLU") == 0) return 'E';
	if (s.compare("PHE") == 0) return 'F';
	if (s.compare("GLY") == 0) return 'G';
	if (s.compare("HIS") == 0) return 'H';
	if (s.compare("ILE") == 0) return 'I';
	if (s.compare("LYS") == 0) return 'K';
	if (s.compare("LEU") == 0) return 'L';
	if (s.compare("MET") == 0) return 'M';
	if (s.compare("ASN") == 0) return 'N';
	if (s.compare("PRO") == 0) return 'P';
	if (s.compare("GLN") == 0) return 'Q';
	if (s.compare("ARG") == 0) return 'R';
	if (s.compare("SER") == 0) return 'S';
	if (s.compare("THR") == 0) return 'T';
	if (s.compare("VAL") == 0) return 'V';
	if (s.compare("TRP") == 0) return 'W';
	if (s.compare("TYR") == 0) return 'Y';
	if (s.compare("GLX") == 0) return 'Z';
	if (s.compare("PTR") == 0) return 'Y';
	if (s.compare("SEC") == 0) return 'C';
	return 'X';
}

// get_seq:
void Protein::get_seq(BioSequence& seq) const
{
	seq.clear();
	seq.set_id(id());
	seq.set_annot(name());
	for (int i = 0; i < nresidues(); i++)
		seq.add(residue_type(i));
}

// clone:
Protein* Protein::clone(Molecule* protein) const
{
	Protein* p;
	if (protein == this) return 0;
	if (!protein) p = new Protein;
	else if (p = dynamic_cast<Protein*>(protein)) p->clear();
	else return 0; // protein is a non-Protein Molecule object

	if (! Molecule::clone(p)) return 0;
	for (int i = 0; i < nresidues(); i++)
		clone_residue(residue(i), *p);
	return p;
}

// clone_residue:
void Protein::clone_residue(const Residue& r, Protein& p) const
{
	Residue cr;

	int iatom = atom_number_of(r.n);
	if (iatom >= 0) cr.n = p.atom(iatom);
	iatom = atom_number_of(r.h);
	if (iatom >= 0) cr.h = p.atom(iatom);
	iatom = atom_number_of(r.ca);
	if (iatom >= 0) cr.ca = p.atom(iatom);
	iatom = atom_number_of(r.c);
	if (iatom >= 0) cr.c = p.atom(iatom);
	iatom = atom_number_of(r.o);
	if (iatom >= 0) cr.o = p.atom(iatom);
	iatom = atom_number_of(r.oxt);
	if (iatom >= 0) cr.oxt = p.atom(iatom);
	iatom = atom_number_of(r.cb);
	if (iatom >= 0) cr.cb = p.atom(iatom);
	for (vector<Atom*>::const_iterator
		i = r.side_chain_atoms.begin(); i != r.side_chain_atoms.end(); ++i) {
		iatom = atom_number_of(*i);
		if (iatom >= 0) cr.side_chain_atoms.push_back(p.atom(iatom));
	}
	cr.name = r.name;
	cr.number = r.number;
	p.add_residue(cr);
}

// Residue::natoms:
int Protein::Residue::natoms() const
{
	int na = 0;
	if (n) na++;
	if (h) na++;
	if (ca) na++;
	if (c) na++;
	if (o) na++;
	if (oxt) na++;
	if (cb) na++;
	na += side_chain_atoms.size();
	return na;
}

// Residue::clear:
void Protein::Residue::clear()
{
	int na = 0;
	n = h = ca = c = o = oxt = cb = 0;
	side_chain_atoms.clear();
	name.erase();
	number.erase();
	resseq = 0;
}

// remove_hydrogens:
void Protein::remove_hydrogens()
{
	Molecule::remove_hydrogens();
	Atom* atom;
	for (vector<Residue>::iterator
		p = m_residues.begin(); p != m_residues.end(); ++p) {
		p->h = 0;
		// other hydrogens are stored in side_chain_atoms
		for (vector<Atom*>::iterator
			q = p->side_chain_atoms.begin(); q != p->side_chain_atoms.end(); ) {
			atom = *q;
			if (atom->element == Atom::get_element_number("H")) {
				q = p->side_chain_atoms.erase(q);
			} else {
				++q;
			}
		}
	}
}

// assign_atom_types:
void Protein::assign_atom_types()
{
	if (m_atom_types_assigned) return;
	remove_hydrogens();
	for (int i = 0; i < nresidues(); i++) {
		if (residue(i).n) set_type_for_n(residue(i).n);
		if (residue(i).ca) set_type_for_ca(residue(i).ca);
		if (residue(i).c) set_type_for_c(residue(i).c);
		if (residue(i).o) set_type_for_o(residue(i).o);
		if (residue(i).oxt) set_type_for_oxt(residue(i).oxt);
		if (residue(i).cb) set_type_for_cb(residue(i).cb);

		// for side-chain atoms
		if (residue(i).name.compare("CYS") == 0) {
			set_cys_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("ASP") == 0) {
			set_asp_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("GLU") == 0) {
			set_glu_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("PHE") == 0) {
			set_phe_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("HIS") == 0) {
			set_his_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("ILE") == 0) {
			set_ile_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("LYS") == 0) {
			set_lys_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("LEU") == 0) {
			set_leu_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("MET") == 0) {
			set_met_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("ASN") == 0) {
			set_asn_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("PRO") == 0) {
			set_pro_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("GLN") == 0) {
			set_gln_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("ARG") == 0) {
			set_arg_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("SER") == 0) {
			set_ser_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("THR") == 0) {
			set_thr_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("VAL") == 0) {
			set_val_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("TRP") == 0) {
			set_trp_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("TYR") == 0) {
			set_tyr_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("PTR") == 0) {
			set_ptr_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("ACE") == 0) {
			set_ace_atom_property_type(residue(i).side_chain_atoms);
		} else if (residue(i).name.compare("NH2") == 0) {
			set_nh2_atom_property_type(residue(i).side_chain_atoms);
		} else {
			set_unknown_property_type(residue(i).side_chain_atoms);
		}
	}
	m_atom_types_assigned = true;
}

// set_type_for_n:
void Protein::set_type_for_n(Atom* a)
{
	a->properties = "0011100100000110";
	a->pc_class = 3;
}

// set_type_for_ca:
void Protein::set_type_for_ca(Atom* a)
{
	a->properties = "0011000010000000";
	a->pc_class = 6;
}

// set_type_for_c:
void Protein::set_type_for_c(Atom* a)
{
	a->properties = "0100100010000000";
	a->pc_class = 7;
}

// set_type_for_o:
void Protein::set_type_for_o(Atom* a)
{
	a->properties = "0100101000000100";
	a->pc_class = 4;
}

// set_type_for_oxt:
void Protein::set_type_for_oxt(Atom* a)
{
	a->properties = "0011101000000100";
	a->pc_class = 2;
}

// set_type_for_cb:
void Protein::set_type_for_cb(Atom* a)
{
	a->properties = "0010001000000000";
	a->pc_class = 6;
}

//
void Protein::set_cys_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" SG ") == 0) {
			a->properties = "0010001000000000";
			a->pc_class = 6;
		}
	}
}

//
void Protein::set_asp_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0100100010000001";
			a->pc_class = 7;
		} else if (a->name.compare(" OD1") == 0) {
			a->properties = "0100101000000100";
			a->pc_class = 2;
		} else if (a->name.compare(" OD2") == 0) {
			a->properties = "0011101000000100";
			a->pc_class = 2;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_glu_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0011000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD ") == 0) {
			a->properties = "0100100010000001";
			a->pc_class = 7;
		} else if (a->name.compare(" OE1") == 0) {
			a->properties = "0100101000000100";
			a->pc_class = 2;
		} else if (a->name.compare(" OE2") == 0) {
			a->properties = "0011101000000100";
			a->pc_class = 2;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_phe_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0101100010101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD1") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD2") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CE1") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CE2") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CZ ") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_his_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0100100010110000";
			a->pc_class = 7;
		} else if (a->name.compare(" ND1") == 0) {
			a->properties = "0011100100110100";
			a->pc_class = 5;
		} else if (a->name.compare(" CD2") == 0) {
			a->properties = "0101100100110000";
			a->pc_class = 6;
		} else if (a->name.compare(" CE1") == 0) {
			a->properties = "0100100100110000";
			a->pc_class = 6;
		} else if (a->name.compare(" NE2") == 0) {
			a->properties = "0101100100110100";
			a->pc_class = 5;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_ile_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG1") == 0) {
			a->properties = "0010001000000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CG2") == 0) {
			a->properties = "0010000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD1") == 0) {
			a->properties = "0010001000000000";
			a->pc_class = 6;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_lys_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0010000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD ") == 0) {
			a->properties = "0010000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CE ") == 0) {
			a->properties = "0010000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" NZ ") == 0) {
			a->properties = "0010001000000100";
			a->pc_class = 1;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_leu_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0010000010000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD1") == 0) {
			a->properties = "0010001000000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD2") == 0) {
			a->properties = "0010001000000000";
			a->pc_class = 6;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_met_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0010000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" SD ") == 0) {
			a->properties = "0010000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CE ") == 0) {
			a->properties = "0010001000000000";
			a->pc_class = 6;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_asn_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0100100010000000";
			a->pc_class = 7;
		} else if (a->name.compare(" OD1") == 0) {
			a->properties = "0100101000000100";
			a->pc_class = 4;
		} else if (a->name.compare(" ND2") == 0) {
			a->properties = "0011101000000110";
			a->pc_class = 3;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_pro_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0010000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD ") == 0) {
			a->properties = "0010000100000000";
			a->pc_class = 6;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_gln_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0011000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD ") == 0) {
			a->properties = "0100100010000000";
			a->pc_class = 7;
		} else if (a->name.compare(" OE1") == 0) {
			a->properties = "0100101000000100";
			a->pc_class = 4;
		} else if (a->name.compare(" NE2") == 0) {
			a->properties = "0011101000000110";
			a->pc_class = 3;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_arg_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0010000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD ") == 0) {
			a->properties = "0010000100000000";
			a->pc_class = 6;
		} else if (a->name.compare(" NE ") == 0) {
			a->properties = "0011100100000100";
			a->pc_class = 1;
		} else if (a->name.compare(" CZ ") == 0) {
			a->properties = "0100100010000000";
			a->pc_class = 7;
		} else if (a->name.compare(" NH1") == 0) {
			a->properties = "0100101000000100";
			a->pc_class = 1;
		} else if (a->name.compare(" NH2") == 0) {
			a->properties = "0011101000000100";
			a->pc_class = 1;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_ser_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" OG ") == 0) {
			a->properties = "0010001000000100";
			a->pc_class = 5;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_thr_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" OG1") == 0) {
			a->properties = "0010001000000100";
			a->pc_class = 5;
		} else if (a->name.compare(" CG2") == 0) {
			a->properties = "0010001000000000";
			a->pc_class = 6;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_val_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG1") == 0) {
			a->properties = "0010001000000000";
			a->pc_class = 6;
		} else if (a->name.compare(" CG2") == 0) {
			a->properties = "0010001000000000";
			a->pc_class = 6;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_trp_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0101100010110000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD1") == 0) {
			a->properties = "0100100100110000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD2") == 0) {
			a->properties = "0101100010111000";
			a->pc_class = 6;
		} else if (a->name.compare(" NE1") == 0) {
			a->properties = "0011100100110100";
			a->pc_class = 3;
		} else if (a->name.compare(" CE2") == 0) {
			a->properties = "0101100010111000";
			a->pc_class = 7;
		} else if (a->name.compare(" CE3") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CZ2") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CZ3") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CH2") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_tyr_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0101100010101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD1") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD2") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CE1") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CE2") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CZ ") == 0) {
			a->properties = "0101100010101000";
			a->pc_class = 7;
		} else if (a->name.compare(" OH ") == 0) {
			a->properties = "0011101000000100";
			a->pc_class = 5;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_ptr_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CG ") == 0) {
			a->properties = "0101100010101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD1") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CD2") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CE1") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CE2") == 0) {
			a->properties = "0101100100101000";
			a->pc_class = 6;
		} else if (a->name.compare(" CZ ") == 0) {
			a->properties = "0101100010101000";
			a->pc_class = 7;
		} else if (a->name.compare(" OH ") == 0) {
			a->properties = "0011100100000100";
			a->pc_class = 4;
		} else if (a->name.compare(" P  ") == 0) {
			a->properties = "0100100001000000";
			a->pc_class = 7;
		} else if (a->name.compare(" O1P") == 0) {
			a->properties = "0011101000000100";
			a->pc_class = 2;
		} else if (a->name.compare(" O2P") == 0) {
			a->properties = "0011101000000100";
			a->pc_class = 2;
		} else if (a->name.compare(" O3P") == 0) {
			a->properties = "0100101000000100";
			a->pc_class = 2;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_ace_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" CH3") == 0) {
			a->properties = "0100100010000001";
			a->pc_class = 7;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_nh2_atom_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		if (a->name.compare(" N  ") == 0) {
			a->properties = "0011101000000110";
			a->pc_class = 3;
		} else {
			a->properties = "0000000000000000";
			a->pc_class = 0;
		}
	}
}

//
void Protein::set_unknown_property_type(vector<Atom*>& atoms)
{
	Atom* a;
	for (vector<Atom*>::iterator
		p = atoms.begin(); p != atoms.end(); ++p) {
		a = *p;
		a->properties = "0000000000000000";
		a->pc_class = 0;
	}
}

//
void Protein::print(std::ostream& os) const
{
	Molecule::print(os);
	print_residues(os);
}

//
void Protein::print_residues(std::ostream& os) const
{
	int n = 0;
	for (int i = 0; i < nresidues(); i++) {
		const Residue& r = m_residues[i];
		n += r.natoms();
		os << "R";
		os << "\t" << std::setw(4) << std::left << i;
		os << "\t" << std::setw(4) << std::left << r.number;
		os << "\t" << r.name;
		os << "\t" << n-1;
		os << std::endl;
	}
}

//
int Protein::residue_containing_atom(const Atom* a) const
{
	for (int i = 0; i < nresidues(); i++) {
		const Residue& r = residue(i);
		if (r.n == a) return i;
		if (r.h == a) return i;
		if (r.ca == a) return i;
		if (r.c == a) return i;
		if (r.o == a) return i;
		if (r.oxt == a) return i;
		if (r.cb == a) return i;
		for (int j = 0; j < r.side_chain_atoms.size(); j++) {
			if (r.side_chain_atoms[j] == a) return i;
		}
	}
	return -1;
}
