#include <vector>
#include <algorithm>
#include "trans_membrane_prediction.h"
#include "hydropathy.h"

using namespace Cba;

const double TransMembranePrediction::C1 = -9.02;
const double TransMembranePrediction::C2 = 14.27;
const double TransMembranePrediction::MIN_H = -C2/C1-1.0;

// predict: according to the algorithm by:
// Petr Klein, Minoru Kanehisa and Charles Delisi (1985).
//   Biochimica et Biophysica Acta, 815: 468-476.
void TransMembranePrediction::predict(const BioSequence& seq, 
	std::vector<Region>& regions)
{
	regions.clear();
	if (seq.length() < LSEG) return;
	assign_hydropathy(seq);
	Region r;
	while (find_tm(r)) {
		extend_tm(r);
		regions.push_back(r);
		mask_region(r.s, r.e);
	}
}

//
void TransMembranePrediction::mask_region(int s, int e)
{
	s -= LSEG - 1;
	s = (s >= 0 ? s : 0);
	e += LSEG - 1;
	e = (e <= m_seq_length-LSEG ? e : m_seq_length-LSEG);
	for (int i = s; i <= e; i++)
		m_hydropathy[i] = MIN_H;
}

// assign_hydropathy:
void TransMembranePrediction::assign_hydropathy(const BioSequence& seq)
{
	m_seq_length = seq.length();

	Hydropathy hindex;
	m_hydropathy.clear();
	m_hydropathy.reserve(m_seq_length);
	for (int i = 0; i < m_seq_length; i++)
		m_hydropathy.push_back(hindex.value_for(seq.at(i)));

	double h;
	for (int i = 0; i < m_seq_length-LSEG+1; i++) {
		h = 0;
		for (int j = i; j < i+LSEG; j++)
			h += m_hydropathy[j];
		m_hydropathy[i] = h/LSEG;
	}
}

// find_tm:
bool TransMembranePrediction::find_tm(Region& r) const
{
	// for each sequence segment of length LSEG
	double h;
	double hmax = MIN_H;
	int s0 = -1;
	for (int i = 0; i < m_seq_length-LSEG+1; i++) {
		if (m_hydropathy[i] > hmax) {
			hmax = m_hydropathy[i];
			s0 = i; 
		}
	}

	if (is_tm(hmax)) {
		r.s0 = s0;
		r.e0 = s0+LSEG-1;
		r.h = hmax;
		return true;
	}
	return false;
}

// extend_tm:
void TransMembranePrediction::extend_tm(Region& r) const
{
	// N-terminal <---
	r.s = r.s0;
	for (int i = r.s0-1; i >= 0; i--)
		if (is_tm(m_hydropathy[i])) r.s--;
		else break;

	// ---> C-terminal
	r.e = r.e0;
	for (int i = r.s0+1; i < (m_hydropathy.size()-LSEG+1); i++)
		if (is_tm(m_hydropathy[i])) r.e++;
		else break;
}

// sort:
void TransMembranePrediction::sort(std::vector<Region>& regions) const
{
	std::sort(regions.begin(), regions.end(), Sort());
}
