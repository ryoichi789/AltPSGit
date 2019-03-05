#ifndef CBA_PREDICT_TM_H
#define CBA_PREDICT_TM_H

#include <functional>
#include <iostream>
#include <vector>
#include "app.h"
#include "bio_sequence.h"
#include "trans_membrane_prediction.h"

namespace Cba {

class PredictTm : public App {
public:
	PredictTm();
	void set_param_default();
	void set_param(char *argv[]);
	int run();

private:
	void usage();
	const char* m_fasta_file; // sequence file in Fasta format
	const char* m_out_file; // file for output

	void print(const BioSequence& seq,
		const std::vector<TransMembranePrediction::Region>& rs,
		std::ostream& os) const;

	class SetParam : public std::unary_function<Param, void> {
	public:
		SetParam(PredictTm* obj) : m_obj(obj) {}
		void operator() (const Param& param);
	private:
		PredictTm* m_obj;
	};
	friend class SetParam;
};

}
#endif
