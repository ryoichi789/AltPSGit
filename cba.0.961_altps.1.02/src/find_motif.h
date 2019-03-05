#ifndef CBA_FIND_MOTIF_H
#define CBA_FIND_MOTIF_H

#include <functional>
#include <iostream>
#include "app.h"
#include "bio_sequence.h"
#include "region.h"

namespace Cba {

class FindMotif : public App {
public:
	FindMotif();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	
	const char* m_motif; // motif in Perl regular expression
	const char* m_fasta_file; // sequence file in Fasta format
	const char* m_out_file; // file for output

	void print_match(
		const BioSequence& seq, const Region& region, std::ostream& os) const;
	static const size_t L_FLANKING_REGION = 5;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(FindMotif* obj) : m_obj(obj) {}
		void operator() (const Param& param);
	private:
		FindMotif* m_obj;
	};
	friend class SetParam;
};

}
#endif
