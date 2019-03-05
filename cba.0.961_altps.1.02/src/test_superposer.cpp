#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include "test_superposer.h"
#include "pdb_reader.h"
#include "pdb_entry.h"
#include "sdf_reader.h"
#include "sdf_entry.h"
#include "molecule.h"
#include "protein.h"
#include "superposer.h"

using namespace Cba;

TestSuperposer the_app;

// TestSuperposer:
TestSuperposer::TestSuperposer() {}

// set_param_default:
void TestSuperposer::set_param_default()
{
	m_prog = "test_superposer";
	m_molfile = 0;
	m_sdf = false;
}

// set_param:
void TestSuperposer::set_param(char *argv[])
{
	try {
		parse_argv(argv);
	} catch(BadArgv) {
		usage();
		std::exit(2);
	}
	std::for_each(m_params.begin(), m_params.end(), SetParam(this));
	if (!m_molfile) {
		usage();
		std::exit(2);
	}
}

// run:
int TestSuperposer::run()
{
	std::ifstream ifs(m_molfile);
	Molecule mol1, mol2;
	Superposer spp;

	if (m_sdf) {
		SdfReader reader(ifs);
		SdfEntry ent;
		while (reader.read(ent)) {
			if (ent.create_molecule(mol1)) {
				mol1.clone(&mol2);
				mol2.shift_by(1.0, 1.0, 1.0);
				superpose_molecules(mol1, mol2, spp);
			}
		}
	} else {
		PdbEntry ent;
		PdbReader::read(ent, ifs);
		for (int i = 0; i < ent.nmolecules(); i++) {
			if (ent.create_molecule(i, mol1)) {
				mol1.clone(&mol2);
				mol2.shift_by(1.0, 1.0, 1.0);
				superpose_molecules(mol1, mol2, spp);
			}
		}
	}

	return 0;
}

// usage:
void TestSuperposer::usage()
{
	using std::cerr;
	using std::endl;
	cerr << "usage: " << m_prog << " [-key value ...]" << endl;
	cerr << "\tpossible keys and values are as follows:" << endl;
	cerr << "\tf: molecule file in PDB (default) or SDF format [string] (unique; mandatory)" << endl;
	cerr << "\tsdf: set true if molecule file is in SDF format [t/f] (unique; f by default)" << endl;
	cerr << "\tdebug: run in debug mode if set true [t/f] (f by default)" << endl;
}

// SetParam::operator():
void TestSuperposer::SetParam::operator()(const Param& param)
{
	using std::strcmp;
	if (strcmp(param.key, "f") == 0) {
		m_obj->m_molfile = param.value;
	} else if (strcmp(param.key, "sdf") == 0) {
		if (param.value[0] == 't') m_obj->m_sdf = true;
		else m_obj->m_sdf = false;
	} else if (strcmp(param.key, "debug") == 0) {
		if (param.value[0] == 't') m_obj->m_debug = true;
		else m_obj->m_debug = false;
	} else {
		std::cerr << "unrecognized parameter " << param.key << std::endl;
		m_obj->usage();
		std::exit(2);
	}
}

// superpose_molecules:
void TestSuperposer::superpose_molecules(
	Molecule& mol1, Molecule& mol2, Superposer& spp) const
{
	using std::cout;
	using std::endl;

	std::vector<Position*> ps1;
	std::vector<Position*> ps2;
	for (int i = 0; i < mol1.natoms(); i++)
		ps1.push_back(&(mol1.atom(i)->pos));
	for (int i = 0; i < mol2.natoms(); i++)
		ps2.push_back(&(mol2.atom(i)->pos));
	double rmsd = spp.superpose(ps1, ps2);

	cout << mol1.id() << ' ' << mol2.id() << ' '
		<< std::fixed << std::setprecision(3) << std::setw(8)
		<< rmsd << endl;
	spp.print(cout);
	cout << endl;
}
