#ifndef CBA_ALIGN_SEQ_H
#define CBA_ALIGN_SEQ_H

#include <functional>
#include "app.h"

namespace Cba {

class AlignSeq : public App {
public:
	AlignSeq();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	bool m_dna; // DNA sequence comparison?
	const char* m_fasta1;
	const char* m_fasta2;
	const char* m_output_file;
	const char* m_score_file; // file containing substitution matrix
	int m_max_nalignments;
	int m_min_score;
	static const int MAXNALIGN = 10;
	static const int MINSCORE = 80;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(AlignSeq* obj) : m_obj(obj) {}
		void operator() (const Param& param);
	private:
		AlignSeq* m_obj;
	};
	friend class SetParam;
};

}
#endif
