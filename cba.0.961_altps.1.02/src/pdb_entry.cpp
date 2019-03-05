#include <vector>
#include <cstring>
#include <algorithm>
#include <sstream>
#include "pdb_entry.h"
#include "atom.h"
#include "molecule.h"
#include "protein.h"
#include "bond_maker.h"
#include "exception.h"

using namespace Cba;

const int PdbEntry::MAX_LEN_PEPTIDE = 10;
const double PdbEntry::MAX_PEPTIDE_BOND_LENGTH = 1.8;
const double PdbEntry::MAX_COVALENT_BOND_LENGTH = 2.0;
const double PdbEntry::OVERLAPPED_DISTANCE = 0.8;

//
PdbEntry::PdbEntry()
	: m_header(0), m_obslte(0), m_title(0), m_caveat(0), m_compnd(0),
		m_source(0), m_keywds(0), m_expdta(0), m_author(0), m_sprsde(0),
		m_jrnl(0), m_remark1(0), m_remark2(0), m_remark3(0), m_cryst1(0),
		m_origx(0), m_scale(0), m_master(0), m_end(0)
{
}

//
PdbEntry::~PdbEntry()
{
	clear();
}

//
PdbEntry::PdbEntry(const PdbEntry& ent)
{
	ent.clone(this);
}

//
PdbEntry& PdbEntry::operator=(const PdbEntry& ent)
{
	if (this != &ent) {
		clear();
		ent.clone(this);
	}
	return *this;
}

// clear:
void PdbEntry::clear()
{
	for (std::vector<PdbRecord*>::iterator
		i = m_records.begin(); i != m_records.end(); ++i)
		delete *i;
	m_records.clear();
	m_header = 0; m_obslte = 0; m_title = 0; m_caveat = 0; m_compnd = 0;
	m_source = 0; m_keywds = 0; m_expdta = 0; m_author = 0;
	m_revdats.clear();
	m_sprsde = 0; m_jrnl = 0; m_remark1 = 0; m_remark2 = 0; m_remark3 = 0;
	m_remarkns.clear(); m_dbrefs.clear(); m_seqadvs.clear();
	m_seqreses.clear(); m_modreses.clear(); m_hets.clear(); m_hetnams.clear();
	m_hetsyns.clear(); m_formuls.clear(); m_helices.clear(); m_sheets.clear();
	m_turns.clear(); m_ssbonds.clear(); m_links.clear(); m_hydbnds.clear();
	m_sltbrgs.clear(); m_cispeps.clear(); m_sites.clear();
	m_cryst1 = 0; m_origx = 0; m_scale = 0;
	m_mtrices.clear(); m_tvects.clear(); m_models.clear(); m_endmdls.clear();
	m_atoms.clear(); m_sigatms.clear(); m_anisous.clear(); m_siguijs.clear();
	m_ters.clear(); m_conects.clear();
	m_master = 0; m_end = 0;
	m_id.erase();
	m_het_names.clear();
	m_unknowns.clear();
}

// make_entry:
void PdbEntry::make_entry()
{
	assign_id();
	remove_altloc();
	define_residues();
	define_molecules();
	make_het_name_list();
}

// assign_id:
void PdbEntry::assign_id()
{
	if (m_header && m_header->idcode.at(0) != ' ') set_id(m_header->idcode);
	else set_id(std::string("NA"));
}

// define_residues;
void PdbEntry::define_residues()
{
	m_residues.clear();
	if (m_atoms.size() == 0) return;
	for (int i = 1; i < m_atoms.size(); i++)
		if ((m_atoms[i]->resseq != m_atoms[i-1]->resseq) ||
			(m_atoms[i]->icode != m_atoms[i-1]->icode) ||
			(m_atoms[i]->chainid != m_atoms[i-1]->chainid) ||
			(m_atoms[i]->resname.compare(m_atoms[i-1]->resname)))
			m_residues.push_back(i-1);
	m_residues.push_back(m_atoms.size()-1); // last residue
}

// define_molecules:
void PdbEntry::define_molecules()
{
	m_molecules.clear();
	typedef std::vector<PdbRecord*>::const_iterator CI;
	CI p = m_records.begin();
	CI q = m_records.end();
	int ia;
	for (int i = 0; i < m_residues.size(); i++) {
		if (covalently_bonded_to_next(i)) continue;
		p = std::find(p, q, m_atoms[ia = m_residues[i]]);
		if (end_of_molecule(ia, ++p))
			m_molecules.push_back(i);
	}
}

// end_of_molecule:
bool PdbEntry::end_of_molecule(int ia, RCI p) const
{
	if (p == m_records.end()) return true;
	// ia = serial number of the last atom of a residue
	// *p = PdbRecord* right after the PdbAtom* for the last atom
	if (ia == m_atoms.size()-1) return true; // last atom
	if (dynamic_cast<PdbTer*>(*p) || dynamic_cast<PdbEndmdl*>(*p) ||
		dynamic_cast<PdbEnd*>(*p))
		return true;
	if (!dynamic_cast<PdbAtom*>(*p)) {
		if (dynamic_cast<PdbAnisou*>(*p)) {
			RCI q = p;
			if (!dynamic_cast<PdbAtom*>(*++q)) return true;
		} else {
			return true;
		}
	}

	if (m_atoms[ia]->chainid != m_atoms[ia+1]->chainid) return true;

	if (is_het_molecule(ia) || is_het_molecule(ia+1)) return true;

	if (m_atoms[ia]->resname != m_atoms[ia+1]->resname) {
		if (m_atoms[ia]->resname == "HOH" || m_atoms[ia]->resname == "DOD") return true;
		if (m_atoms[ia+1]->resname == "HOH" || m_atoms[ia+1]->resname == "DOD") return true;
	}

	return false;
}

//
bool PdbEntry::is_het_molecule(int ia) const
{
	// each of HETATMs defined in HET records is considered an independent molecule
	// except for PTR, SEC, MSE, ACE, NH2, (see 1d00.pdb for example)
	if (m_atoms[ia]->hetatm &&
		m_atoms[ia]->resname.compare("PTR") && m_atoms[ia]->resname.compare("SEC") &&
		m_atoms[ia]->resname.compare("MSE") && m_atoms[ia]->resname.compare("ACE") &&
		m_atoms[ia]->resname.compare("SEP") && m_atoms[ia]->resname.compare("NH2")) {
		for (int i = 0; i < m_hets.size(); i++) {
			if (m_atoms[ia]->resname.compare(m_hets[i]->hetid) == 0 &&
				m_atoms[ia]->chainid == m_hets[i]->chainid &&
				m_atoms[ia]->resseq == m_hets[i]->seqnum &&
				m_atoms[ia]->icode == m_hets[i]->icode)
				return true;
		}
	}
	return false;
}

// title:
std::string PdbEntry::title() const
{
	if (!m_title) return std::string("NA");
	return m_title->title;
}

// compound:
std::string PdbEntry::compound() const
{
	if (!m_compnd) return std::string("NA");
	return m_compnd->compound;
}

// nresidues:
int PdbEntry::nresidues(int imol) const
{
	if (imol < 0 || imol >= m_molecules.size())
		throw Exception("molecule out of range");
	if (imol == 0) return m_molecules[imol] + 1;
	return m_molecules[imol] - m_molecules[imol-1];
}

// natoms: of a molecule
int PdbEntry::natoms(int imol) const
{
	if (imol < 0 || imol >= m_molecules.size())
		throw Exception("molecule out of range");
	if (imol == 0) return m_residues[m_molecules[imol]] + 1;
	return m_residues[m_molecules[imol]] - m_residues[m_molecules[imol-1]];
}

// natoms: of a residue of a molecule
int PdbEntry::natoms(int imol, int ires) const
{
	if (imol < 0 || imol >= m_molecules.size())
		throw Exception("molecule out of range");
	if (ires < 0 || ires >= nresidues(imol))
		throw Exception("residue out of range");
	if (imol == 0 && ires == 0) return m_residues[0] + 1;
	int jres = residue_begin(imol) + ires;
	return m_residues[jres] - m_residues[jres-1];
}

// residue_begin:
int PdbEntry::residue_begin(int imol) const
{
	int n = 0;
	for (int i = 0; i < imol; i++)
		n += nresidues(i);
	return n;
}

// atom_begin:
int PdbEntry::atom_begin(int imol) const
{
	int i = residue_begin(imol); // first residue
	if (i == 0) return 0;
	return m_residues[i-1] + 1;
}

// chainid:
char PdbEntry::chainid(int imol) const
{
	if (imol < 0 || imol >= m_molecules.size())
		throw Exception("molecule out of range");
	return m_atoms[m_residues[m_molecules[imol]]]->chainid;
}

// resname:
const std::string& PdbEntry::resname(int imol, int ires) const
{
	if (imol < 0 || imol >= m_molecules.size())
		throw Exception("molecule out of range");
	if (ires < 0 || ires >= nresidues(imol))
		throw Exception("residue out of range");
	if (imol == 0 && ires == 0) return m_atoms[0]->resname;
	int jres = residue_begin(imol) + ires;
	return m_atoms[m_residues[jres]]->resname;
}

// resnum:
std::string PdbEntry::resnum(int imol, int ires) const
{
	if (imol < 0 || imol >= m_molecules.size())
		throw Exception("molecule out of range");
	if (ires < 0 || ires >= nresidues(imol))
		throw Exception("residue out of range");

	int a = m_residues[residue_begin(imol)+ires]; //last atom of residue
	std::ostringstream oss;
	oss << m_atoms[a]->resseq << m_atoms[a]->icode;
	return oss.str();
}

// resseq:
int PdbEntry::resseq(int imol, int ires) const
{
	if (imol < 0 || imol >= m_molecules.size())
		throw Exception("molecule out of range");
	if (ires < 0 || ires >= nresidues(imol))
		throw Exception("residue out of range");

	int a = m_residues[residue_begin(imol)+ires]; //last atom of residue
	return m_atoms[a]->resseq;
}

// type_of_molecule:
PdbEntry::MoleculeType PdbEntry::type_of_molecule(int imol) const
{
	if (imol < 0 || imol >= m_molecules.size())
		throw Exception("molecule out of range");
	if (nresidues(imol) == 1 && natoms(imol) == 1) return PdbEntry::ATOM;

	// protein or peptide?
	int n = 0; // number of C-alphas
	for (int i = atom_begin(imol); i <= m_residues[m_molecules[imol]]; i++) {
		if (m_atoms[i]->name.compare(1, 3, "CA ") == 0) n++;
		if (n > MAX_LEN_PEPTIDE) return PdbEntry::PROTEIN;
	}
	if (n >= 2) return PdbEntry::PEPTIDE; // oligopeptides
	if (n == 1) return PdbEntry::SMALL; // amino acid

	// DNA, RNA, or water?
	n = 0; // number of nucleotides
	int nt = 0; // number of Ts
	int nu = 0; // number of Us
	int nw = 0; // number of waters
	for (int i = residue_begin(imol); i <= m_molecules[imol]; i++) {
		const std::string& rn = m_atoms[m_residues[i]]->resname;
		if (rn.compare("  A") == 0 || rn.compare(" +A") == 0) ++n;
		else if (rn.compare("  T") == 0 || rn.compare(" +T") == 0) { ++n; ++nt; }
		else if (rn.compare("  U") == 0 || rn.compare(" +U") == 0) { ++n; ++nu; }
		else if (rn.compare("  G") == 0 || rn.compare(" +G") == 0) ++n;
		else if (rn.compare("  C") == 0 || rn.compare(" +C") == 0) ++n;
		else if (rn.compare("  I") == 0 || rn.compare(" +I") == 0) ++n;
		else if (rn.compare("HOH") == 0 || rn.compare("DOD") == 0) ++nw;
	}
	if (n >= 2) {
		if (nu == 0) return PdbEntry::DNA;
		else return PdbEntry::RNA;
	}
	if (nw == nresidues(imol)) return PdbEntry::WATER;
	/*
	if (nresidues(imol) == 1) {
		std::map<std::string, std::string>::const_iterator p;
		if ((p = m_het_names.find(get_small_molecule_id(imol))) != m_het_names.end())
			return PdbEntry::SMALL;
	}
	return PdbEntry::UNDEFINED;
	*/
	return PdbEntry::SMALL;
}

// nmr_model:
int PdbEntry::nmr_model(int imol) const
{
	if (m_models.empty()) return 0; // not an NMR entry
	if (imol < 0 || imol >= m_molecules.size())
		throw Exception("molecule out of range");

	/*
	std::vector<PdbRecord*>::const_iterator a0
		= std::find(m_records.begin(), m_records.end(), m_atoms[atom_begin(imol)]);
		// start of molecule
	*/

	int n = 0;
	PdbEndmdl* em;
	PdbRecord* first_atom = m_atoms[atom_begin(imol)];
	for (std::vector<PdbRecord*>::const_iterator p = m_records.begin(); p != m_records.end(); ++p) {
		if (*p == first_atom) return n+1;
		if (em = dynamic_cast<PdbEndmdl*>(*p)) n++;
	}

	/*
	std::vector<PdbRecord*>::const_iterator p0 = m_records.begin();
	std::vector<PdbRecord*>::const_iterator p1;
	PdbModel* m;
	for (int i = 0; i < m_models.size(); i++) {
		p0 = std::find(p0, m_records.end(), m_models[i]); // start of model
		if (p0 == m_records.end()) throw Exception("model record not found");

		if (i < m_endmdls.size()) {
			p1 = std::find(p0, m_records.end(), m_endmdls[i]); // end of model
		} else {
			// this can happen, e.g., in 1ANP
			p1 = std::find(p0, m_records.end(), m_end); // end of model
		}
		//if (p1 == m_records.end()) throw Exception("endmdl record not found");

		if (p0 < a0 && a0 < p1) {
			//m = dynamic_cast<PdbModel*>(*p0);
			//return m->serial;
			// some NMR entries starts model numbering from 0, while most starts from 1;
			return (i+1);
		}
		p0 = p1;
	}
	*/

	return 0;
}

// resolution:
float PdbEntry::resolution() const
{
	PdbRemark2* rec;
	for (std::vector<PdbRecord*>::const_iterator
		p = m_records.begin(); p != m_records.end(); ++p) {
		if (rec = dynamic_cast<PdbRemark2*>(*p)) {
			return rec->resolution;
		}
	}
	return 0.0;
}

// create_molecule:
bool PdbEntry::create_molecule(int imol, Molecule& mol) const
{
	mol.clear();

	// id
	std::ostringstream oss;
	oss << id();
	if (chainid(imol) == ' ') oss << '_';
	else oss << chainid(imol);
	oss << '_' << imol;
	mol.set_id(oss.str());

	// name
	mol.set_name(m_atoms[atom_begin(imol)]->resname);

	// type
	int t = type_of_molecule(imol);
	if (t == PdbEntry::PROTEIN) mol.set_type(Molecule::PROTEIN);
	else if (t == PdbEntry::PEPTIDE) mol.set_type(Molecule::PEPTIDE);
	else if (t == PdbEntry::DNA) mol.set_type(Molecule::DNA);
	else if (t == PdbEntry::RNA) mol.set_type(Molecule::RNA);
	else if (t == PdbEntry::SMALL) mol.set_type(Molecule::SMALL);
	else if (t == PdbEntry::WATER) mol.set_type(Molecule::WATER);
	else if (t == PdbEntry::ATOM) mol.set_type(Molecule::ATOM);
	else mol.set_type(Molecule::UNDEFINED);

	Atom a;
	PdbAtom* pa;
	for (int i = atom_begin(imol); i <= m_residues[m_molecules[imol]]; i++) {
		pa = m_atoms[i];
		a.name = pa->name;
		a.number = pa->serial;
		//a.element = get_element_number(a.name);
		a.element = get_element_number(pa);
		a.pos.move_to(pa->x, pa->y, pa->z);
		mol.add_atom(a);
	}

	if (t == PdbEntry::SMALL || t == PdbEntry::PEPTIDE) {
		std::string name;
		name = get_small_molecule_id(imol) + " = "
			+ get_small_molecule_name(get_small_molecule_id(imol));
		mol.set_name(name);
		BondMaker bm;
		bm.make_bonds(mol);
	}

	return true;
}

// create_protein:
bool PdbEntry::create_protein(int imol, Protein& protein) const
{
	if (type_of_molecule(imol) != PROTEIN &&
		type_of_molecule(imol) != PEPTIDE) return false;
	create_molecule(imol, protein);
	if (m_title) protein.set_name(m_title->title);
	else if (m_compnd) protein.set_name(m_compnd->compound);

	// define residues of protein
	Protein::Residue residue;
	Atom* atom;
	int a = 0; // serial number of atom
	for (int i = 0; i < nresidues(imol); i++) {
		residue.clear();
		for (int j = 0; j < natoms(imol, i); j++) {
			atom = protein.atom(a++);
			if (atom->name.compare(" N  ") == 0) residue.n = atom;
			else if (atom->name.compare(" H  ") == 0) residue.h = atom;
			else if (atom->name.compare(" CA ") == 0) residue.ca = atom;
			else if (atom->name.compare(" C  ") == 0) residue.c = atom;
			else if (atom->name.compare(" O  ") == 0) residue.o = atom;
			else if (atom->name.compare(" OXT") == 0) residue.oxt = atom;
			else if (atom->name.compare(" CB ") == 0) residue.cb = atom;
			else if (atom->name.compare("N   ") == 0) residue.n = atom;
			else if (atom->name.compare("H   ") == 0) residue.h = atom;
			else if (atom->name.compare("CA  ") == 0) residue.ca = atom;
			else if (atom->name.compare("C   ") == 0) residue.c = atom;
			else if (atom->name.compare("O   ") == 0) residue.o = atom;
			else if (atom->name.compare("OXT ") == 0) residue.oxt = atom;
			else if (atom->name.compare("CB  ") == 0) residue.cb = atom;
			else residue.side_chain_atoms.push_back(atom);
		}
		residue.name = resname(imol, i);
		residue.number = resnum(imol, i);
		residue.resseq = resseq(imol, i);
		protein.add_residue(residue);
	}

	return true;
}

//
//unsigned char PdbEntry::get_element_number(const std::string& atom_name) const
unsigned char PdbEntry::get_element_number(const PdbAtom* atom) const
{
	if (atom->element == " H") return Atom::get_element_number("H");
	if (atom->element == " C") return Atom::get_element_number("C");
	if (atom->element == " N") return Atom::get_element_number("N");
	if (atom->element == " O") return Atom::get_element_number("O");
	if (atom->element == " P") return Atom::get_element_number("P");
	if (atom->element == " S") return Atom::get_element_number("S");
	if (atom->element == " F") return Atom::get_element_number("F");
	if (atom->element == " I") return Atom::get_element_number("I");
	if (atom->element == " K") return Atom::get_element_number("K");
	if (atom->element == "BR") return Atom::get_element_number("Br");
	if (atom->element == "CL") return Atom::get_element_number("Cl");
	if (atom->element == "NA") return Atom::get_element_number("Na");
	if (atom->element == "CA") return Atom::get_element_number("Ca");
	if (atom->element == "MG") return Atom::get_element_number("Mg");
	if (atom->element == "MN") return Atom::get_element_number("Mn");
	if (atom->element == "FE") return Atom::get_element_number("Fe");
	if (atom->element == "ZN") return Atom::get_element_number("Zn");
	if (atom->element == "HG") return Atom::get_element_number("Hg");

	const std::string& atom_name = atom->name;
	if (atom_name.at(1) == 'H') return Atom::get_element_number("H");
	if (atom_name.at(1) == 'C') return Atom::get_element_number("C");
	if (atom_name.at(1) == 'N') return Atom::get_element_number("N");
	if (atom_name.at(1) == 'O') return Atom::get_element_number("O");
	if (atom_name.at(1) == 'P') return Atom::get_element_number("P");
	if (atom_name.at(1) == 'S') return Atom::get_element_number("S");
	if (atom_name.at(1) == 'I') return Atom::get_element_number("I");
	if (atom_name.compare(0,2,"BR") == 0) return Atom::get_element_number("Br");
	if (atom_name.compare(0,2,"CL") == 0) return Atom::get_element_number("Cl");
	if (atom_name.compare("NA  ") == 0) return Atom::get_element_number("Na");
	if (atom_name.compare("CA  ") == 0) return Atom::get_element_number("Ca");
	if (atom_name.compare("MG  ") == 0) return Atom::get_element_number("Mg");
	if (atom_name.compare("MN  ") == 0) return Atom::get_element_number("Mn");
	if (atom_name.compare("FE  ") == 0) return Atom::get_element_number("Fe");
	if (atom_name.compare("ZN  ") == 0) return Atom::get_element_number("Zn");
	if (atom_name.compare("HG  ") == 0) return Atom::get_element_number("HG");
	return 0;
}

//
void PdbEntry::make_het_name_list()
{
	std::string name;

	m_het_names.clear();
	std::map<std::string, std::string>::iterator p;

	PdbHet* het;
	for (std::vector<PdbHet*>::const_iterator q = m_hets.begin(); q != m_hets.end(); ++q) {
		het = *q;
		if ((p = m_het_names.find(het->hetid)) == m_het_names.end()) { // first occurrence
			if (het->text.empty() || het->text.compare(" ") == 0) m_het_names[het->hetid]; // empty string
			else m_het_names[het->hetid] = het->text + "; ";
		}
	}

	PdbHetnam* hetnam;
	for (std::vector<PdbHetnam*>::const_iterator q = m_hetnams.begin(); q != m_hetnams.end(); ++q) {
		hetnam = *q;
		if ((p = m_het_names.find(hetnam->hetid)) == m_het_names.end()) continue;
		if (hetnam->text.empty() || hetnam->text.compare(" ") == 0) continue;
		m_het_names[hetnam->hetid] += hetnam->text + "; ";
	}

	PdbHetsyn* hetsyn;
	for (std::vector<PdbHetsyn*>::const_iterator q = m_hetsyns.begin(); q != m_hetsyns.end(); ++q) {
		hetsyn = *q;
		if ((p = m_het_names.find(hetsyn->hetid)) == m_het_names.end()) continue;
		if (hetsyn->hetsynonyms.empty() || hetsyn->hetsynonyms.compare(" ") == 0) continue;
		m_het_names[hetsyn->hetid] += hetsyn->hetsynonyms + "; ";
	}

	// dummy
	m_het_names["NA"];

	for (p = m_het_names.begin(); p != m_het_names.end(); ++p) {
		std::string& names = p->second;
		int i = names.rfind("; ");
		if (i < names.length()) names.erase(i);
	}
}

//
const std::string& PdbEntry::get_small_molecule_name(const std::string& id) const
{
	std::map<std::string, std::string>::const_iterator p;
	if ((p = m_het_names.find(id)) == m_het_names.end())
		p = m_het_names.find("NA");
	return p->second;
}

//
bool PdbEntry::peptide_bonded_to_next(int ir) const
{
	if (ir == m_residues.size()-1) return false;
	int ia0, ia1;
	const PdbAtom* c = 0;
	const PdbAtom* n = 0;

	if (ir == 0) ia0 = 0;
	else ia0 = m_residues[ir-1]+1;
	ia1 = m_residues[ir];
	for (int i = ia0; i <= ia1; i++) {
		if (m_atoms[i]->name == " C  ") {
			c = m_atoms[i];
			break;
		}
	}
	if (!c) return false;

	ia0 = ia1 + 1;
	ia1 = m_residues[ir+1];
	for (int i = ia0; i <= ia1; i++) {
		if (m_atoms[i]->name == " N  ") {
			n = m_atoms[i];
			break;
		}
	}
	if (!n) return false;

	double dx = std::fabs(c->x - n->x);
	if (dx > MAX_PEPTIDE_BOND_LENGTH) return false;
	double dy = std::fabs(c->y - n->y);
	if (dx > MAX_PEPTIDE_BOND_LENGTH) return false;
	double dz = std::fabs(c->z - n->z);
	if (dx > MAX_PEPTIDE_BOND_LENGTH) return false;
	double dd = dx*dx + dy*dy + dz*dz;
	if (dd > MAX_PEPTIDE_BOND_LENGTH*MAX_PEPTIDE_BOND_LENGTH) return false;

	return true;
}

//
bool PdbEntry::covalently_bonded_to_next(int ir) const
{
	if (ir == m_residues.size()-1) return false;

	int ia00, ia01;
	if (ir == 0) ia00 = 0;
	else ia00 = m_residues[ir-1]+1;
	ia01 = m_residues[ir];
	if (is_metal(m_atoms[ia00])) return false;

	int ia10, ia11;
	ia10 = ia01 + 1;
	ia11 = m_residues[ir+1];
	if (is_metal(m_atoms[ia10])) return false;

	const PdbAtom* a0 = 0;
	const PdbAtom* a1 = 0;
	double dx, dy, dz, dd;
	double ddmax = MAX_COVALENT_BOND_LENGTH*MAX_COVALENT_BOND_LENGTH;
	double ddmin = OVERLAPPED_DISTANCE*OVERLAPPED_DISTANCE;

	for (int i = ia00; i <= ia01; i++) {
		a0 = m_atoms[i];
		for (int j = ia10; j <= ia11; j++) {
			a1 = m_atoms[j];
			dx = std::fabs(a0->x - a1->x);
			if (dx > MAX_COVALENT_BOND_LENGTH) continue;
			dy = std::fabs(a0->y - a1->y);
			if (dy > MAX_COVALENT_BOND_LENGTH) continue;
			dz = std::fabs(a0->z - a1->z);
			if (dz > MAX_COVALENT_BOND_LENGTH) continue;
			dd = dx*dx + dy*dy + dz*dz;
			if (dd < ddmin) return false; // cf. 1agm, where two conformers of ACR are listed.
		}
	}

	for (int i = ia00; i <= ia01; i++) {
		a0 = m_atoms[i];
		for (int j = ia10; j <= ia11; j++) {
			a1 = m_atoms[j];
			dx = std::fabs(a0->x - a1->x);
			if (dx > MAX_COVALENT_BOND_LENGTH) continue;
			dy = std::fabs(a0->y - a1->y);
			if (dy > MAX_COVALENT_BOND_LENGTH) continue;
			dz = std::fabs(a0->z - a1->z);
			if (dz > MAX_COVALENT_BOND_LENGTH) continue;
			dd = dx*dx + dy*dy + dz*dz;
			if (dd > ddmax) continue;
			return true;
		}
	}

	return false;
}

//
bool PdbEntry::is_metal(const PdbAtom* atom) const
{
	if (atom->element == "NA") return true;
	if (atom->element == "CA") return true;
	if (atom->element == "MG") return true;
	if (atom->element == "MN") return true;
	if (atom->element == "FE") return true;
	if (atom->element == "ZN") return true;
	if (atom->element == "HG") return true;

	if (atom->name == "NA  ") return true;
	if (atom->name == "CA  ") return true;
	if (atom->name == "MG  ") return true;
	if (atom->name == "MN  ") return true;
	if (atom->name == "FE  ") return true;
	if (atom->name == "ZN  ") return true;
	if (atom->name == "HG  ") return true;

	if (atom->resname == "HOH") return true; // though HOH is actually not a metal ...

	return false;
}

//
void PdbEntry::remove_altloc()
{
	if (m_atoms.size() == 0) return;

	PdbAtom* atom0;
	PdbAtom* atom1;
	std::vector<PdbAtom*>::iterator p0 = m_atoms.begin();
	std::vector<PdbAtom*>::iterator p1 = m_atoms.begin();
	for (++p1; p1 != m_atoms.end(); ) {
		atom0 = *p0;
		atom1 = *p1;
		if (atom0->name == atom1->name &&
			atom0->resname == atom1->resname &&
			atom0->chainid == atom1->chainid &&
			atom0->resseq == atom1->resseq &&
			atom0->icode == atom1->icode &&
			atom0->altloc != atom1->altloc) {
			p1 = m_atoms.erase(p1);
		} else {
			++p0;
			++p1;
		}
	}
}

//
int PdbEntry::get_depdate() const
{
	int depdate = 0;
	int k;
	const std::string& sdepdate = m_header->depdate;
	if (sdepdate.length() < 9) return depdate;

	// DD-MMM-YY -> YYYYMMDD

	// year
	if (sdepdate.at(7) == '6' || sdepdate.at(7) == '7' ||
		sdepdate.at(7) == '8' || sdepdate.at(7) == '9') {
		k = std::atoi(sdepdate.substr(7,2).c_str());
		k += 1900;
		depdate += k * 10000;
	} else if (sdepdate.at(7) == '0') {
		k = std::atoi(sdepdate.substr(8,1).c_str());
		k += 2000;
		depdate += k * 10000;
	} else {
		return depdate;
	}

	// month
	if (sdepdate.substr(3,3) == "JAN") depdate += 100;
	else if (sdepdate.substr(3,3) == "FEB") depdate += 200;
	else if (sdepdate.substr(3,3) == "MAR") depdate += 300;
	else if (sdepdate.substr(3,3) == "APR") depdate += 400;
	else if (sdepdate.substr(3,3) == "MAY") depdate += 500;
	else if (sdepdate.substr(3,3) == "JUN") depdate += 600;
	else if (sdepdate.substr(3,3) == "JUL") depdate += 700;
	else if (sdepdate.substr(3,3) == "AUG") depdate += 800;
	else if (sdepdate.substr(3,3) == "SEP") depdate += 900;
	else if (sdepdate.substr(3,3) == "OCT") depdate += 1000;
	else if (sdepdate.substr(3,3) == "NOV") depdate += 1100;
	else if (sdepdate.substr(3,3) == "DEC") depdate += 1200;
	else return depdate;

	// day of month
	if (sdepdate.at(0) == '0') k = 0;
	else if (sdepdate.at(0) == '1') k = 10;
	else if (sdepdate.at(0) == '2') k = 20;
	else if (sdepdate.at(0) == '3') k = 30;
	else return depdate;
	k += std::atoi(sdepdate.substr(1,1).c_str());
	depdate += k;

	return depdate;
}

//
char PdbEntry::get_seq_chain_id(int i) const
{
	if (i < 0 || i >= m_seqreses.size()) return 0;
	return m_seqreses[i]->chainid;
}

//
int PdbEntry::length_of_seqres(int i) const
{
	if (i < 0 || i >= m_seqreses.size()) return 0;
	return m_seqreses[i]->numres;
}

//
bool PdbEntry::get_ith_sequence(BioSequence& seq, int i) const
{
	seq.clear();
	if (i < 0 || i >= m_seqreses.size()) return false;
	return get_protein_sequence(seq, m_seqreses[i]->chainid);
}

//
bool PdbEntry::get_protein_sequence(BioSequence& seq, char chain) const
{
	using std::vector;
	using std::string;

	seq.clear();
	string s = id();
	s.append(1, chain);
	seq.set_id(s);
	seq.set_annot(compound());

	const PdbSeqres* seqres = 0;
	for (vector<PdbSeqres*>::const_iterator p = m_seqreses.begin(); p != m_seqreses.end(); ++p) {
		if ((*p)->chainid == chain) {
			seqres = *p;
			break;
		}
	}
	if (! seqres) return false;

	const vector<string>& resname = seqres->resname;
	for (vector<string>::const_iterator p = resname.begin(); p != resname.end(); ++p) {
		if (*p == "ALA") seq.add('A');
		else if (*p == "ASX") seq.add('B');
		else if (*p == "CYS") seq.add('C');
		else if (*p == "ASP") seq.add('D');
		else if (*p == "GLU") seq.add('E');
		else if (*p == "PHE") seq.add('F');
		else if (*p == "GLY") seq.add('G');
		else if (*p == "HIS") seq.add('H');
		else if (*p == "ILE") seq.add('I');
		else if (*p == "LYS") seq.add('K');
		else if (*p == "LEU") seq.add('L');
		else if (*p == "MET") seq.add('M');
		else if (*p == "ASN") seq.add('N');
		else if (*p == "PRO") seq.add('P');
		else if (*p == "GLN") seq.add('Q');
		else if (*p == "ARG") seq.add('R');
		else if (*p == "SER") seq.add('S');
		else if (*p == "THR") seq.add('T');
		else if (*p == "VAL") seq.add('V');
		else if (*p == "TRP") seq.add('W');
		else if (*p == "TYR") seq.add('Y');
		else if (*p == "GLX") seq.add('Z');
		else if (*p == "SEC") seq.add('C');
		else seq.add('X');
	}

	return true;
}

//
PdbEntry* PdbEntry::clone(PdbEntry* ent) const
{
	if (ent == this) return 0;
	if (ent) ent->clear();
	else ent = new PdbEntry;

	ent->m_id = m_id;
	ent->m_residues = m_residues;
	ent->m_molecules = m_molecules;
	ent->m_het_names = m_het_names;

	if (m_header) {
		ent->m_header = new PdbHeader;
		*(ent->m_header) = *m_header;
	}
	if (m_obslte) {
		ent->m_obslte = new PdbObslte;
		*(ent->m_obslte) = *m_obslte;
	}
	if (m_title) {
		ent->m_title = new PdbTitle;
		*(ent->m_title) = *m_title;
	}
	if (m_caveat) {
		ent->m_caveat = new PdbCaveat;
		*(ent->m_caveat) = *m_caveat;
	}
	if (m_compnd) {
		ent->m_compnd = new PdbCompnd;
		*(ent->m_compnd) = *m_compnd;
	}
	if (m_source) {
		ent->m_source = new PdbSource;
		*(ent->m_source) = *m_source;
	}
	if (m_keywds) {
		ent->m_keywds = new PdbKeywds;
		*(ent->m_keywds) = *m_keywds;
	}
	if (m_expdta) {
		ent->m_expdta = new PdbExpdta;
		*(ent->m_expdta) = *m_expdta;
	}
	if (m_author) {
		ent->m_author = new PdbAuthor;
		*(ent->m_author) = *m_author;
	}
	for (int i = 0; i < m_revdats.size(); i++) {
		if (m_revdats[i]) {
			ent->m_revdats.push_back(new PdbRevdat);
			*(ent->m_revdats.back()) = *(m_revdats[i]);
		} else {
			ent->m_revdats.push_back(0);
		}
	}
	if (m_sprsde) {
		ent->m_sprsde = new PdbSprsde;
		*(ent->m_sprsde) = *m_sprsde;
	}
	if (m_jrnl) {
		ent->m_jrnl = new PdbJrnl;
		*(ent->m_jrnl) = *m_jrnl;
	}
	if (m_remark1) {
		ent->m_remark1 = new PdbRemark1;
		*(ent->m_remark1) = *m_remark1;
	}
	if (m_remark2) {
		ent->m_remark2 = new PdbRemark2;
		*(ent->m_remark2) = *m_remark2;
	}
	if (m_remark3) {
		ent->m_remark3 = new PdbRemark3;
		*(ent->m_remark3) = *m_remark3;
	}
	for (int i = 0; i < m_remarkns.size(); i++) {
		if (m_remarkns[i]) {
			ent->m_remarkns.push_back(new PdbRemarkN);
			*(ent->m_remarkns.back()) = *(m_remarkns[i]);
		} else {
			ent->m_remarkns.push_back(0);
		}
	}
	for (int i = 0; i < m_dbrefs.size(); i++) {
		if (m_dbrefs[i]) {
			ent->m_dbrefs.push_back(new PdbDbref);
			*(ent->m_dbrefs.back()) = *(m_dbrefs[i]);
		} else {
			ent->m_dbrefs.push_back(0);
		}
	}
	for (int i = 0; i < m_seqadvs.size(); i++) {
		if (m_seqadvs[i]) {
			ent->m_seqadvs.push_back(new PdbSeqadv);
			*(ent->m_seqadvs.back()) = *(m_seqadvs[i]);
		} else {
			ent->m_seqadvs.push_back(0);
		}
	}
	for (int i = 0; i < m_seqreses.size(); i++) {
		if (m_seqreses[i]) {
			ent->m_seqreses.push_back(new PdbSeqres);
			*(ent->m_seqreses.back()) = *(m_seqreses[i]);
		} else {
			ent->m_seqreses.push_back(0);
		}
	}
	for (int i = 0; i < m_modreses.size(); i++) {
		if (m_modreses[i]) {
			ent->m_modreses.push_back(new PdbModres);
			*(ent->m_modreses.back()) = *(m_modreses[i]);
		} else {
			ent->m_modreses.push_back(0);
		}
	}
	for (int i = 0; i < m_hets.size(); i++) {
		if (m_hets[i]) {
			ent->m_hets.push_back(new PdbHet);
			*(ent->m_hets.back()) = *(m_hets[i]);
		} else {
			ent->m_hets.push_back(0);
		}
	}
	for (int i = 0; i < m_hetnams.size(); i++) {
		if (m_hetnams[i]) {
			ent->m_hetnams.push_back(new PdbHetnam);
			*(ent->m_hetnams.back()) = *(m_hetnams[i]);
		} else {
			ent->m_hetnams.push_back(0);
		}
	}
	for (int i = 0; i < m_hetsyns.size(); i++) {
		if (m_hetsyns[i]) {
			ent->m_hetsyns.push_back(new PdbHetsyn);
			*(ent->m_hetsyns.back()) = *(m_hetsyns[i]);
		} else {
			ent->m_hetsyns.push_back(0);
		}
	}
	for (int i = 0; i < m_formuls.size(); i++) {
		if (m_formuls[i]) {
			ent->m_formuls.push_back(new PdbFormul);
			*(ent->m_formuls.back()) = *(m_formuls[i]);
		} else {
			ent->m_formuls.push_back(0);
		}
	}
	for (int i = 0; i < m_helices.size(); i++) {
		if (m_helices[i]) {
			ent->m_helices.push_back(new PdbHelix);
			*(ent->m_helices.back()) = *(m_helices[i]);
		} else {
			ent->m_helices.push_back(0);
		}
	}
	for (int i = 0; i < m_sheets.size(); i++) {
		if (m_sheets[i]) {
			ent->m_sheets.push_back(new PdbSheet);
			*(ent->m_sheets.back()) = *(m_sheets[i]);
		} else {
			ent->m_sheets.push_back(0);
		}
	}
	for (int i = 0; i < m_turns.size(); i++) {
		if (m_turns[i]) {
			ent->m_turns.push_back(new PdbTurn);
			*(ent->m_turns.back()) = *(m_turns[i]);
		} else {
			ent->m_turns.push_back(0);
		}
	}
	for (int i = 0; i < m_ssbonds.size(); i++) {
		if (m_ssbonds[i]) {
			ent->m_ssbonds.push_back(new PdbSsbond);
			*(ent->m_ssbonds.back()) = *(m_ssbonds[i]);
		} else {
			ent->m_ssbonds.push_back(0);
		}
	}
	for (int i = 0; i < m_links.size(); i++) {
		if (m_links[i]) {
			ent->m_links.push_back(new PdbLink);
			*(ent->m_links.back()) = *(m_links[i]);
		} else {
			ent->m_links.push_back(0);
		}
	}
	for (int i = 0; i < m_hydbnds.size(); i++) {
		if (m_hydbnds[i]) {
			ent->m_hydbnds.push_back(new PdbHydbnd);
			*(ent->m_hydbnds.back()) = *(m_hydbnds[i]);
		} else {
			ent->m_hydbnds.push_back(0);
		}
	}
	for (int i = 0; i < m_sltbrgs.size(); i++) {
		if (m_sltbrgs[i]) {
			ent->m_sltbrgs.push_back(new PdbSltbrg);
			*(ent->m_sltbrgs.back()) = *(m_sltbrgs[i]);
		} else {
			ent->m_sltbrgs.push_back(0);
		}
	}
	for (int i = 0; i < m_cispeps.size(); i++) {
		if (m_cispeps[i]) {
			ent->m_cispeps.push_back(new PdbCispep);
			*(ent->m_cispeps.back()) = *(m_cispeps[i]);
		} else {
			ent->m_cispeps.push_back(0);
		}
	}
	for (int i = 0; i < m_sites.size(); i++) {
		if (m_sites[i]) {
			ent->m_sites.push_back(new PdbSite);
			*(ent->m_sites.back()) = *(m_sites[i]);
		} else {
			ent->m_sites.push_back(0);
		}
	}
	if (m_cryst1) {
		ent->m_cryst1 = new PdbCryst1;
		*(ent->m_cryst1) = *m_cryst1;
	}
	if (m_origx) {
		ent->m_origx = new PdbOrigx;
		*(ent->m_origx) = *m_origx;
	}
	if (m_scale) {
		ent->m_scale = new PdbScale;
		*(ent->m_scale) = *m_scale;
	}
	for (int i = 0; i < m_mtrices.size(); i++) {
		if (m_mtrices[i]) {
			ent->m_mtrices.push_back(new PdbMtrix);
			*(ent->m_mtrices.back()) = *(m_mtrices[i]);
		} else {
			ent->m_mtrices.push_back(0);
		}
	}
	for (int i = 0; i < m_tvects.size(); i++) {
		if (m_tvects[i]) {
			ent->m_tvects.push_back(new PdbTvect);
			*(ent->m_tvects.back()) = *(m_tvects[i]);
		} else {
			ent->m_tvects.push_back(0);
		}
	}
	for (int i = 0; i < m_models.size(); i++) {
		if (m_models[i]) {
			ent->m_models.push_back(new PdbModel);
			*(ent->m_models.back()) = *(m_models[i]);
		} else {
			ent->m_models.push_back(0);
		}
	}
	for (int i = 0; i < m_endmdls.size(); i++) {
		if (m_endmdls[i]) {
			ent->m_endmdls.push_back(new PdbEndmdl);
			*(ent->m_endmdls.back()) = *(m_endmdls[i]);
		} else {
			ent->m_endmdls.push_back(0);
		}
	}
	for (int i = 0; i < m_atoms.size(); i++) {
		if (m_atoms[i]) {
			ent->m_atoms.push_back(new PdbAtom);
			*(ent->m_atoms.back()) = *(m_atoms[i]);
		} else {
			ent->m_atoms.push_back(0);
		}
	}
	for (int i = 0; i < m_sigatms.size(); i++) {
		if (m_sigatms[i]) {
			ent->m_sigatms.push_back(new PdbSigatm);
			*(ent->m_sigatms.back()) = *(m_sigatms[i]);
		} else {
			ent->m_sigatms.push_back(0);
		}
	}
	for (int i = 0; i < m_anisous.size(); i++) {
		if (m_anisous[i]) {
			ent->m_anisous.push_back(new PdbAnisou);
			*(ent->m_anisous.back()) = *(m_anisous[i]);
		} else {
			ent->m_anisous.push_back(0);
		}
	}
	for (int i = 0; i < m_siguijs.size(); i++) {
		if (m_siguijs[i]) {
			ent->m_siguijs.push_back(new PdbSiguij);
			*(ent->m_siguijs.back()) = *(m_siguijs[i]);
		} else {
			ent->m_siguijs.push_back(0);
		}
	}
	for (int i = 0; i < m_ters.size(); i++) {
		if (m_ters[i]) {
			ent->m_ters.push_back(new PdbTer);
			*(ent->m_ters.back()) = *(m_ters[i]);
		} else {
			ent->m_ters.push_back(0);
		}
	}
	for (int i = 0; i < m_conects.size(); i++) {
		if (m_conects[i]) {
			ent->m_conects.push_back(new PdbConect);
			*(ent->m_conects.back()) = *(m_conects[i]);
		} else {
			ent->m_conects.push_back(0);
		}
	}
	if (m_master) {
		ent->m_master = new PdbMaster;
		*(ent->m_master) = *m_master;
	}
	if (m_end) {
		ent->m_end = new PdbEnd;
		*(ent->m_end) = *m_end;
	}
	for (int i = 0; i < m_unknowns.size(); i++) {
		if (m_unknowns[i]) {
			ent->m_unknowns.push_back(new PdbUnknownRecord);
			*(ent->m_unknowns.back()) = *(m_unknowns[i]);
		} else {
			ent->m_unknowns.push_back(0);
		}
	}

	ent->m_records.resize(m_records.size());
	for (int i = 0; i < m_records.size(); i++) {
		if (m_records[i]) {
			if (m_records[i] == m_header) {
				ent->m_records[i] = ent->m_header;
			} else if (m_records[i] == m_obslte) {
				ent->m_records[i] = ent->m_obslte;
			} else if (m_records[i] == m_title) {
				ent->m_records[i] = ent->m_title;
			} else if (m_records[i] == m_caveat) {
				ent->m_records[i] = ent->m_caveat;
			} else if (m_records[i] == m_compnd) {
				ent->m_records[i] = ent->m_compnd;
			} else if (m_records[i] == m_source) {
				ent->m_records[i] = ent->m_source;
			} else if (m_records[i] == m_keywds) {
				ent->m_records[i] = ent->m_keywds;
			} else if (m_records[i] == m_expdta) {
				ent->m_records[i] = ent->m_expdta;
			} else if (m_records[i] == m_author) {
				ent->m_records[i] = ent->m_author;
			} else if (m_records[i] == m_sprsde) {
				ent->m_records[i] = ent->m_sprsde;
			} else if (m_records[i] == m_jrnl) {
				ent->m_records[i] = ent->m_jrnl;
			} else if (m_records[i] == m_remark1) {
				ent->m_records[i] = ent->m_remark1;
			} else if (m_records[i] == m_remark2) {
				ent->m_records[i] = ent->m_remark2;
			} else if (m_records[i] == m_remark3) {
				ent->m_records[i] = ent->m_remark3;
			} else if (m_records[i] == m_cryst1) {
				ent->m_records[i] = ent->m_cryst1;
			} else if (m_records[i] == m_origx) {
				ent->m_records[i] = ent->m_origx;
			} else if (m_records[i] == m_scale) {
				ent->m_records[i] = ent->m_scale;
			} else if (m_records[i] == m_master) {
				ent->m_records[i] = ent->m_master;
			} else if (m_records[i] == m_end) {
				ent->m_records[i] = ent->m_end;
			} else if (dynamic_cast<PdbRevdat*>(m_records[i])) {
				int j;
				for (j = 0; j < m_revdats.size(); j++) {
					if (m_records[i] == m_revdats[j]) {
						ent->m_records[i] = ent->m_revdats[j];
						break;
					}
				}
				if (j == m_revdats.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbRemarkN*>(m_records[i])) {
				int j;
				for (j = 0; j < m_remarkns.size(); j++) {
					if (m_records[i] == m_remarkns[j]) {
						ent->m_records[i] = ent->m_remarkns[j];
						break;
					}
				}
				if (j == m_remarkns.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbDbref*>(m_records[i])) {
				int j;
				for (j = 0; j < m_dbrefs.size(); j++) {
					if (m_records[i] == m_dbrefs[j]) {
						ent->m_records[i] = ent->m_dbrefs[j];
						break;
					}
				}
				if (j == m_dbrefs.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbSeqadv*>(m_records[i])) {
				int j;
				for (j = 0; j < m_seqadvs.size(); j++) {
					if (m_records[i] == m_seqadvs[j]) {
						ent->m_records[i] = ent->m_seqadvs[j];
						break;
					}
				}
				if (j == m_seqadvs.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbSeqres*>(m_records[i])) {
				int j;
				for (j = 0; j < m_seqreses.size(); j++) {
					if (m_records[i] == m_seqreses[j]) {
						ent->m_records[i] = ent->m_seqreses[j];
						break;
					}
				}
				if (j == m_seqreses.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbModres*>(m_records[i])) {
				int j;
				for (j = 0; j < m_modreses.size(); j++) {
					if (m_records[i] == m_modreses[j]) {
						ent->m_records[i] = ent->m_modreses[j];
						break;
					}
				}
				if (j == m_modreses.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbHet*>(m_records[i])) {
				int j;
				for (j = 0; j < m_hets.size(); j++) {
					if (m_records[i] == m_hets[j]) {
						ent->m_records[i] = ent->m_hets[j];
						break;
					}
				}
				if (j == m_hets.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbHetnam*>(m_records[i])) {
				int j;
				for (j = 0; j < m_hetnams.size(); j++) {
					if (m_records[i] == m_hetnams[j]) {
						ent->m_records[i] = ent->m_hetnams[j];
						break;
					}
				}
				if (j == m_hetnams.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbHetsyn*>(m_records[i])) {
				int j;
				for (j = 0; j < m_hetsyns.size(); j++) {
					if (m_records[i] == m_hetsyns[j]) {
						ent->m_records[i] = ent->m_hetsyns[j];
						break;
					}
				}
				if (j == m_hetsyns.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbFormul*>(m_records[i])) {
				int j;
				for (j = 0; j < m_formuls.size(); j++) {
					if (m_records[i] == m_formuls[j]) {
						ent->m_records[i] = ent->m_formuls[j];
						break;
					}
				}
				if (j == m_formuls.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbHelix*>(m_records[i])) {
				int j;
				for (j = 0; j < m_helices.size(); j++) {
					if (m_records[i] == m_helices[j]) {
						ent->m_records[i] = ent->m_helices[j];
						break;
					}
				}
				if (j == m_helices.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbSheet*>(m_records[i])) {
				int j;
				for (j = 0; j < m_sheets.size(); j++) {
					if (m_records[i] == m_sheets[j]) {
						ent->m_records[i] = ent->m_sheets[j];
						break;
					}
				}
				if (j == m_sheets.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbTurn*>(m_records[i])) {
				int j;
				for (j = 0; j < m_turns.size(); j++) {
					if (m_records[i] == m_turns[j]) {
						ent->m_records[i] = ent->m_turns[j];
						break;
					}
				}
				if (j == m_turns.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbSsbond*>(m_records[i])) {
				int j;
				for (j = 0; j < m_ssbonds.size(); j++) {
					if (m_records[i] == m_ssbonds[j]) {
						ent->m_records[i] = ent->m_ssbonds[j];
						break;
					}
				}
				if (j == m_ssbonds.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbLink*>(m_records[i])) {
				int j;
				for (j = 0; j < m_links.size(); j++) {
					if (m_records[i] == m_links[j]) {
						ent->m_records[i] = ent->m_links[j];
						break;
					}
				}
				if (j == m_links.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbHydbnd*>(m_records[i])) {
				int j;
				for (j = 0; j < m_hydbnds.size(); j++) {
					if (m_records[i] == m_hydbnds[j]) {
						ent->m_records[i] = ent->m_hydbnds[j];
						break;
					}
				}
				if (j == m_hydbnds.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbSltbrg*>(m_records[i])) {
				int j;
				for (j = 0; j < m_sltbrgs.size(); j++) {
					if (m_records[i] == m_sltbrgs[j]) {
						ent->m_records[i] = ent->m_sltbrgs[j];
						break;
					}
				}
				if (j == m_sltbrgs.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbCispep*>(m_records[i])) {
				int j;
				for (j = 0; j < m_cispeps.size(); j++) {
					if (m_records[i] == m_cispeps[j]) {
						ent->m_records[i] = ent->m_cispeps[j];
						break;
					}
				}
				if (j == m_cispeps.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbSite*>(m_records[i])) {
				int j;
				for (j = 0; j < m_sites.size(); j++) {
					if (m_records[i] == m_sites[j]) {
						ent->m_records[i] = ent->m_sites[j];
						break;
					}
				}
				if (j == m_sites.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbMtrix*>(m_records[i])) {
				int j;
				for (j = 0; j < m_mtrices.size(); j++) {
					if (m_records[i] == m_mtrices[j]) {
						ent->m_records[i] = ent->m_mtrices[j];
						break;
					}
				}
				if (j == m_mtrices.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbTvect*>(m_records[i])) {
				int j;
				for (j = 0; j < m_tvects.size(); j++) {
					if (m_records[i] == m_tvects[j]) {
						ent->m_records[i] = ent->m_tvects[j];
						break;
					}
				}
				if (j == m_tvects.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbModel*>(m_records[i])) {
				int j;
				for (j = 0; j < m_models.size(); j++) {
					if (m_records[i] == m_models[j]) {
						ent->m_records[i] = ent->m_models[j];
						break;
					}
				}
				if (j == m_models.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbEndmdl*>(m_records[i])) {
				int j;
				for (j = 0; j < m_endmdls.size(); j++) {
					if (m_records[i] == m_endmdls[j]) {
						ent->m_records[i] = ent->m_endmdls[j];
						break;
					}
				}
				if (j == m_endmdls.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbAtom*>(m_records[i])) {
				int j;
				for (j = 0; j < m_atoms.size(); j++) {
					if (m_records[i] == m_atoms[j]) {
						ent->m_records[i] = ent->m_atoms[j];
						break;
					}
				}
				if (j == m_atoms.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbSigatm*>(m_records[i])) {
				int j;
				for (j = 0; j < m_sigatms.size(); j++) {
					if (m_records[i] == m_sigatms[j]) {
						ent->m_records[i] = ent->m_sigatms[j];
						break;
					}
				}
				if (j == m_sigatms.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbAnisou*>(m_records[i])) {
				int j;
				for (j = 0; j < m_anisous.size(); j++) {
					if (m_records[i] == m_anisous[j]) {
						ent->m_records[i] = ent->m_anisous[j];
						break;
					}
				}
				if (j == m_anisous.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbSiguij*>(m_records[i])) {
				int j;
				for (j = 0; j < m_siguijs.size(); j++) {
					if (m_records[i] == m_siguijs[j]) {
						ent->m_records[i] = ent->m_siguijs[j];
						break;
					}
				}
				if (j == m_siguijs.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbTer*>(m_records[i])) {
				int j;
				for (j = 0; j < m_ters.size(); j++) {
					if (m_records[i] == m_ters[j]) {
						ent->m_records[i] = ent->m_ters[j];
						break;
					}
				}
				if (j == m_ters.size()) ent->m_records[i] = 0;
			} else if (dynamic_cast<PdbConect*>(m_records[i])) {
				int j;
				for (j = 0; j < m_conects.size(); j++) {
					if (m_records[i] == m_conects[j]) {
						ent->m_records[i] = ent->m_conects[j];
						break;
					}
				}
				if (j == m_conects.size()) ent->m_records[i] = 0;
			} else {
				int j;
				for (j = 0; j < m_unknowns.size(); j++) {
					if (m_records[i] == m_unknowns[j]) {
						ent->m_records[i] = ent->m_unknowns[j];
						break;
					}
				}
				if (j == m_unknowns.size()) ent->m_records[i] = 0;
			}
		} else {
			ent->m_records[i] = 0;
		}
	}

	return ent;
}
