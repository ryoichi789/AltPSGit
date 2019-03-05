#include <vector>
#include <cmath>
#include <algorithm>
#include "asa_calculator.h"
#include "atom.h"

using namespace Cba;
using std::vector;

const double AsaCalculator::SOLVENT_RADIUS = 1.4;
const double AsaCalculator::ATOM_SLICING_STEP = 0.2;
const double AsaCalculator::PI = 4.0 * std::atan(1.0);

// calc:
void AsaCalculator::calc(vector<Atom*>& atoms, double rsolv, double dz)
{
	set_params(atoms, rsolv, dz);

	// for each atom
	vector<double> zs;
	double ri2, taccarc, z, r, r2;
	for (int iatom = 0; iatom < atoms.size(); iatom++) {
		ignore_atoms_for(iatom, atoms); // ignore distant atoms

		// z-coordinates of circles for atom spheres
		get_zs(iatom, zs, atoms);
		ri2 = m_rs[iatom] * m_rs[iatom];
		taccarc = 0.0; // total accessible arc
		// for each circle of iatom
		for (vector<double>::const_iterator
			q = zs.begin(); q != zs.end(); ++q) {
			z = *q;
			r2 = ri2 -
				(atoms[iatom]->pos.z() - z) * (atoms[iatom]->pos.z() - z);
			r = (r2 < 0.0) ? 0.0 : std::sqrt(r2); // radius of circle for iatom
			taccarc += get_accarc_in_circle(iatom, z, r, atoms);
		}
		
		// asa = sum of r * sin(dt) * accarc * (dzstep / sin(dt))
		//	= sum of r * accarc * dzstep
		//	= r * taccarc * dzstep */
		atoms[iatom]->asa = m_rs[iatom] * m_zstep * taccarc;
	}
}

// set_params:
void AsaCalculator::set_params(const vector<Atom*>& atoms, double rsolv, double dz)
{
	if (rsolv > 0.0) m_rsolv = rsolv;
	else m_rsolv = AsaCalculator::SOLVENT_RADIUS;

	if (dz > 0.0) m_zstep = dz;
	else m_zstep = AsaCalculator::ATOM_SLICING_STEP;

	m_igna.clear();
	m_rs.resize(atoms.size());
	for (int i = 0; i < atoms.size(); i++) {
		m_rs[i] = atoms[i]->vdwr() + m_rsolv;
		m_igna.push_back(false);
	}

	m_rmax = 0.0;
	for (int i = 0; i < atoms.size(); i++)
		if (m_rs[i] > m_rmax) m_rmax = m_rs[i];
}

// ignore_atoms_for:
void AsaCalculator::ignore_atoms_for(size_t iatom, const vector<Atom*>& atoms)
{
	double dx, dy, dz, dd;
	const Position& pi = atoms[iatom]->pos;
	for (int jatom = 0; jatom < atoms.size(); jatom++) {
		const Position& pj = atoms[jatom]->pos;
		if (jatom == iatom) { m_igna[jatom] = true; continue; }

		dx = std::fabs(pi.x() - pj.x());
		if (dx > 2.0 * m_rmax) { m_igna[jatom] = true; continue; }

		dy = std::fabs(pi.y() - pj.y());
		if (dy > 2.0 * m_rmax) { m_igna[jatom] = true; continue; }

		dz = std::fabs(pi.z() - pj.z());
		if (dz > 2.0 * m_rmax) { m_igna[jatom] = true; continue; }

		dd = dx*dx + dy*dy + dz*dz;
		if (dd > 4.0 * m_rmax * m_rmax) { m_igna[jatom] = true; continue; }

		m_igna[jatom] = false;
	}
}

// get_zs:
void AsaCalculator::get_zs(size_t iatom, vector<double>& zs,
	const vector<Atom*>& atoms)
{
	//double a = atoms[iatom]->pos.z() - m_rs[iatom] + m_zstep;
	double a = atoms[iatom]->pos.z() - m_rs[iatom] + 0.5*m_zstep;
	double b = atoms[iatom]->pos.z() + m_rs[iatom];
	zs.clear();
	do {
		zs.push_back(a);
		a += m_zstep;
	} while(a < b);
}
		
// get_accarc_in_circle:
double AsaCalculator::get_accarc_in_circle(
	int iatom, double z, double r, vector<Atom*>& atoms)
{
	vector<Arc> barcs; // arcs covered (buried) by other atoms
	double rr, r2;
	for (int jatom = 0; jatom < atoms.size(); jatom++) {
		if (m_igna[jatom]) continue;
		const Position& pj = atoms[jatom]->pos;
		rr = m_rs[jatom]*m_rs[jatom] - (pj.z()-z)*(pj.z()-z);
		if (rr <= 0.0) continue;
		r2 = std::sqrt(rr); // radius of sphere for jatom
		// circle-circle overlap
		barcs.push_back(cc_overlap(
			atoms[iatom]->pos.x(), atoms[iatom]->pos.y(), r,
			pj.x(), pj.y(), r2));
	}
			
	// combine overlapping arcs
	vector<Arc> barcs2; // non-overlapping arcs covered by other atoms
	combine_arcs(barcs, barcs2);

	// accessible arcs
	double accarc;
	if (barcs2.empty()) {
		accarc = 2.0 * PI;
	} else {
		accarc = 0.0;
		accarc += barcs2[0].ti;
		for (int i = 1; i < barcs2.size(); i++)
			accarc += (barcs2[i].ti - barcs2[i-1].te);
		accarc += (2.0*PI - barcs2[barcs2.size()-1].te);
	}
	if(accarc > 2.0*PI) accarc = 2.0 * PI;

	return accarc;
}

// find a circle-circle overlap;
// returns an arc in circle 1 covered by circle 2
AsaCalculator::Arc
AsaCalculator::cc_overlap(
	double x1, double y1, double r1, double x2, double y2, double r2)
{
	Arc barc;
	double dd12 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
	double d12 = std::sqrt(dd12);
	if (((r1+r2) <= d12) || ((d12+r2) <= r1)) return barc;
		// circle 1 and circle 2 do not intersect, or circle 1 includes circle 2
	if((d12+r1) <= r2) { // circle 1 is included within circle 2
		barc.ti = 0.0; barc.te = 2.0 * PI;
		return barc;
	}

	double dx = x2 - x1;
	double dy = y2 - y1;
	double mt = std::acos(dx/d12);
	if (dy < 0) mt = 2.0*PI - mt;
	double a = r1*r1 + dd12 - r2*r2;
	double b = 2.0*r1*d12;
	double c = a/b;
	if(c > 1.0) c = 1.0;
	double s = std::acos(c);
	double t1 = mt + s;
	if (t1 > 2.0*PI) t1 -= 2.0*PI;
	double t2 = mt - s;
	if (t2 < 0.0) t2 += 2.0*PI;
	if (t1 > t2) { double d = t1; t1 = t2; t2 = d; }
	double t3 = (t1+t2) / 2.0;
	double t4 = t3 - PI;
	if (t4 < 0.0) t4 += 2.0*PI;
	double x3 = r1*std::cos(t3) + x1;
	double y3 = r1*std::sin(t3) + y1;
	double x4 = r1*std::cos(t4) + x1;
	double y4 = r1*std::sin(t4) + y1;
	double d3 = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2);
	double d4 = (x4-x2)*(x4-x2) + (y4-y2)*(y4-y2);
	if (d3 <= r2*r2) {
		barc.ti = t1; barc.te = t2;
	} else if (d4 <= r2*r2) {
		barc.ti = t2; barc.te = t1;
	} else {
		if(d3 < d4) { barc.ti = t1; barc.te = t2; }
		else { barc.ti = t2; barc.te = t1; }
	}

	return barc;
}

// combine overlapping arcs
void AsaCalculator::combine_arcs(const vector<Arc>& arcs0, vector<Arc>& arcs)
{
	// divide an arc into two if it crosses over the angle 0
	arcs.clear();
	for (vector<Arc>::const_iterator p = arcs0.begin(); p != arcs0.end(); ++p) {
		const Arc& arc0 = *p;
		if (arc0.ti < 0.0) continue;
		if (arc0.ti < arc0.te) {
			arcs.push_back(arc0);
		} else {
			Arc a;
			a.ti = 0.0; a.te = arc0.te;
			arcs.push_back(a);
			a.ti = arc0.ti; a.te = 2.0 * PI;
			arcs.push_back(a);
		}
	}

	// remove an arc if it is included within another arc
	vector<Arc> arcs_tmp;
	for (int i = 0; i < arcs.size(); i++) {
		bool included = false;
		for (int j = 0; j < arcs.size(); j++) {
			if (i == j) continue;
			if ((arcs[j].ti <= arcs[i].ti) && (arcs[i].te <= arcs[j].te)) {
				if ((arcs[j].ti != arcs[i].ti) || (arcs[i].te != arcs[j].te)) {
					included = true;
					break;
				}
			}
		}
		if (! included) arcs_tmp.push_back(arcs[i]);
	}
	std::sort(arcs_tmp.begin(), arcs_tmp.end(), CompArcs()); // sort arcs

	// combine overlapping arcs
	arcs.clear();
	for (int i = 0; i < arcs_tmp.size(); i++) {
		bool overlapping = false;
		int n = arcs.size();
		for (int j = 0; j < n; j++) {
			if(arcs_tmp[i].ti <= arcs[j].te) {
				arcs[j].te = arcs_tmp[i].te;
				overlapping = true;
				break;
			}
		}
		if (! overlapping) arcs.push_back(arcs_tmp[i]);
	}
}	

//
bool AsaCalculator::CompArcs::operator()(const Arc& a1, const Arc& a2) const
{
	if (a1.ti < a2.ti) return true;
	if ((a1.ti == a2.ti) && (a1.te < a2.te)) return true;
	return false;
}
