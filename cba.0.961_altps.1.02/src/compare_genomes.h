#ifndef CBA_COMPARE_GENOMES_H
#define CBA_COMPARE_GENOMES_H

#include <functional>
#include <iostream>
#include <fstream>
#include "app.h"

namespace Cba {

class CompareGenomes : public App {
public:
	CompareGenomes();
	void set_param_default();
	void set_param(char *argv[]);
	int run();
private:
	void usage();
	static const int MIN_LEN_EQUIV_SUBSEQ = 15;

	const char* m_dna_file1;
	const char* m_dna_file2;
	const char* m_output_file;
	int m_min_len_equiv_subseq;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(CompareGenomes* obj) : m_obj(obj) {}
		void operator() (const Param& param);
	private:
		CompareGenomes* m_obj;
	};
	friend class SetParam;
};

}

#endif
