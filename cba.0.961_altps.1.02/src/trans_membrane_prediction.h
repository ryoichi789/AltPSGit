#ifndef CBA_TRANS_MEMBRANE_PREDICTION_H
#define CBA_TRANS_MEMBRANE_PREDICTION_H

#include <vector>
#include <functional>
#include "bio_sequence.h"

namespace Cba {

class TransMembranePrediction {
public:
	struct Region {
		int s0; // start of core hydrophobic region
		int e0; // end of core hydrophobic region
		int s; // start of extended hydrophobic region
		int e; // end of extended hydrophobic region
		double h; // hydrophobicity of region
	};

	void predict(const BioSequence& seq, std::vector<Region>& regions);
	void sort(std::vector<Region>& regions) const; // by position in sequence

private:
	std::vector<double> m_hydropathy; // for each residue
	int m_seq_length;

	// constants used in prediction
	static const int LSEG = 17;
	static const double C1; // = -9.02;
	static const double C2; // = 14.27;
	// if C1*h + C2 < 0.0,
	// segment is predicted to be trans-membrane
	// (h = mean hydropathy of sequence segment of length LSEG)
	static const double MIN_H; // = -C2/C1-1.0;

	void assign_hydropathy(const BioSequence& seq);
	bool find_tm(Region& r) const;
	void extend_tm(Region& r) const;
	bool is_tm(double h) const;
	void mask_region(int s, int e);

	struct Sort : std::binary_function<Region, Region, bool> {
		bool operator()(const Region& r1, const Region& r2) const
		{ return r1.s0 < r2.s0; }
	};
};

inline bool TransMembranePrediction::is_tm(double h) const
{ return (C1*h + C2) < 0.0; }

}
#endif
