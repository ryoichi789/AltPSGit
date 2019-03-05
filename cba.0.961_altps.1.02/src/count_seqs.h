#ifndef CBA_COUNT_SEQS_H
#define CBA_COUNT_SEQS_H

#include <functional>
#include <string>
#include <vector>
#include "app.h"

namespace Cba {

class CountSeqs : public App {
public:
	CountSeqs();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	std::string m_progname; // file name without directory path
	const char* m_seqfile;
	std::vector<std::string> m_ids; // ID codes of sequences to get

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(CountSeqs* obj) : m_obj(obj) {}
		void operator()(const Param& param);
	private:
		CountSeqs* m_obj;
	};
	friend class SetParam;
};

}
#endif
