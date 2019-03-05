#ifndef CBA_MOL_FIN_CMP_H
#define CBA_MOL_FIN_CMP_H

#include <functional>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include "app.h"
#include "mol_fin.h"

namespace Cba {

class MolFinCmp : public App {
public:
	MolFinCmp();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();

	const char* m_fquery;
	const char* m_ftargetlist;
	std::vector<std::string> m_ftargets;
	const char* m_fout;
	double m_tc_cut;
	double m_pv_cut;
	int m_top_n;
	int m_nqueries;
	int m_ntargets;

	struct Hit {
		std::string idt;
		int fileno;
		double tc;
		double pv;
	};
	std::list<Hit> m_hits;
	int m_hits_capacity;
	static const int HITS_CAPACITY;
	struct CmpHit : std::binary_function<Hit,Hit,bool> {
		bool operator()(const Hit& h1, const Hit& h2) const
		{ return h1.tc > h2.tc; }
	};

	void set_targets();
	void print_search_info(std::ostream& os);
	void get_molid(const char* file, int i, const MolFin& mf, std::string& id) const;
	void print_hits(std::ostream& os, const std::string& idq, int nq, int nt);
	void scan_targets(int nq, const std::string& idq, MolFin& mfq, std::ostream& os);

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(MolFinCmp* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		MolFinCmp* m_obj;
	};
	friend class SetParam;
};

}
#endif
