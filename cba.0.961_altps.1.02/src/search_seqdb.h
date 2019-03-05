#ifndef CBA_SEARCH_SEQDB_H
#define CBA_SEARCH_SEQDB_H

#include <functional>
#include "app.h"

namespace Cba {

class SearchSeqdb : public App {
public:
	SearchSeqdb();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	const char* m_query_fasta_file;
	const char* m_target_fasta_file;
	const char* m_output_file;
	int m_score_min;
	int m_length_min;
	bool m_make_tables; // make segment tables if true

	static const int SCORE_MIN = 80;
	static const int LENGTH_MIN = 50;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(SearchSeqdb* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		SearchSeqdb* m_obj;
	};
	friend class SetParam;
};

}
#endif
