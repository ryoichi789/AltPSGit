#ifndef CBA_FIND_ORFS_H
#define CBA_FIND_ORFS_H

#include <functional>
#include <iostream>
#include "app.h"
#include "dna_sequence.h"

namespace Cba {

class FindOrfs : public App {
public:
	FindOrfs();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	char* m_seq_file;
	char* m_out_file;
	bool m_print_fasta;
	int m_min_len;
	static const int MIN_LEN = 20;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(FindOrfs* obj) : m_obj(obj) {}
		void operator() (const Param& param);
	private:
		FindOrfs* m_obj;
	};
	friend class SetParam;

	void print_orfs(DnaSequence& dna, std::ostream& os) const;
};

}
#endif
