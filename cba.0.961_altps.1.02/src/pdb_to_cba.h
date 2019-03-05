/*
 * pdb_to_cba:
 *	converts protein structure data in PDB format
 *	into CBA-formatted data with physico-chemical properties of atoms added
 */
#ifndef CBA_PDB_TO_CBA_H
#define CBA_PDB_TO_CBA_H

#include <list>
#include "app.h"
#include "pdb_entry.h"
#include "protein.h"

namespace Cba {

class PdbToCba : public App {
public:
	PdbToCba();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

protected:
	void usage();

private:
	char* m_pdb_file; // input
	char* m_cba_file; // output
	bool m_asa; // calculate asa if true
	bool m_protein_only; // process only proteins if true
	std::list<char> m_chains; // target chains; if empty, all chains
	PdbEntry m_pdb;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(PdbToCba* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		PdbToCba* m_obj;
	};
	friend class SetParam;

	bool is_target(int imol) const;
	void make_cba_file_name();
};

}
#endif
