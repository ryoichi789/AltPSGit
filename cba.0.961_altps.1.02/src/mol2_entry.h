#ifndef CBA_MOL2_ENTRY_H
#define CBA_MOL2_ENTRY_H

#include <vector>
#include <string>
#include "molecule.h"
#include "atom.h"

namespace Cba {

class Mol2Reader; // forward declaration

class Mol2Entry {
public:
	Mol2Entry();
	~Mol2Entry();

	const std::string& name() const;
	const std::string& comment() const;
	const std::string& type() const;
	int natoms() const;
	int nbonds() const;
	bool create_molecule(Molecule& mol) const; // returns false in failure
	void clear();

private:
	struct MoleculeRecord { // @<TRIPOS>MOLECULE record
		std::string name;
		int natoms;
		int nbonds;
		std::string type;
		std::string charge_type;
		std::string status_bit;
		std::string comment;

		MoleculeRecord() : natoms(0), nbonds(0) {}
		void clear();
	};

	struct AtomRecord { // @<TRIPOS>ATOM record
		int id;
		std::string name;
		double x;
		double y;
		double z;
		std::string type;
		int subst_id;
		std::string subst_name;
		double charge;
		std::string status_bit;
	};

	struct BondRecord { // @<TRIPOS>BOND record
		int id;
		int atom1;
		int atom2;
		std::string type;
		std::string status_bit;
	};

	MoleculeRecord m_mol;
	std::vector<AtomRecord> m_atoms;
	std::vector<BondRecord> m_bonds;

	void set_atom_element(Atom& a, const std::string& atom_type) const;
	int get_atom_serial_number(int atom1) const;
	int get_bond_type(const std::string& bond_type) const;

	friend class Mol2Reader;
};

inline const std::string& Mol2Entry::name() const { return m_mol.name; }
inline const std::string& Mol2Entry::comment() const { return m_mol.comment; }
inline const std::string& Mol2Entry::type() const { return m_mol.type; }
inline int Mol2Entry::natoms() const { return m_atoms.size(); }
inline int Mol2Entry::nbonds() const { return m_bonds.size(); }

}
#endif
