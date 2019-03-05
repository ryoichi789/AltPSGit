#ifndef CBA_POSITION_H
#define CBA_POSITION_H

#include <cmath>
#include "matrix.h"

namespace Cba {

// Position represents position of a point in three-dimensional space
class Position {
public:
	Position();
	Position(double x, double y, double z);
	Position(const Position& p);
	Position& operator=(const Position& p);
	~Position() {}
	void clear();

	double x() const;
	double y() const;
	double z() const;

	Position& move_to(double x, double y, double z);
	Position& move_to(const Position& p);
	Position& shift_by(double x, double y, double z);
	Position& shift_by(const Position& p);
	Position& rotate_by(const Matrix& rot); // rot must be a 3x3 rotation matrix
	Position& operator+=(const Position& p);
	Position& operator-=(const Position& p);
	Position& operator*=(double a);
	Position& operator/=(double a);

	double distance_from(double x = 0, double y = 0, double z = 0) const;
	double distance_from(const Position& p) const;
	double square_distance_from(double x = 0, double y = 0, double z = 0) const;
	double square_distance_from(const Position& p) const;
	bool is_within_range(const Position& pos, double range) const;

	double distance_from_neighbor(const Position& p, double range) const;
		// if p is within range, returns distance; else returns -1
	double square_distance_from_neighbor(const Position& p, double range) const;
		// if p is within range, returns square distance; else returns -1

	double size() const; // returns length from origin to this point
	Position unit() const; // returns 'unit vector' from 0 toward this point

private:
	double m_x;
	double m_y;
	double m_z;

	friend Position operator+(const Position& p1, const Position& p2);
	friend Position operator-(const Position& p1, const Position& p2);
};

inline double Position::x() const { return m_x; }
inline double Position::y() const { return m_y; }
inline double Position::z() const { return m_z; }
inline Position& Position::move_to(const Position& p)
{ return move_to(p.m_x, p.m_y, p.m_z); }   
inline Position& Position::shift_by(const Position& p)
{ return shift_by(p.m_x, p.m_y, p.m_z); }
inline Position& Position::operator+=(const Position& p) { return shift_by(p); }
inline Position& Position::operator-=(const Position& p)
{ return shift_by(-p.m_x, -p.m_y, -p.m_z); }
inline double Position::distance_from(double x, double y, double z) const
{ return std::sqrt(square_distance_from(x, y, z)); }
inline double Position::distance_from(const Position& p) const
{ return std::sqrt(square_distance_from(p)); }
inline double Position::square_distance_from(const Position& p) const
{ return square_distance_from(p.m_x, p.m_y, p.m_z); }
inline double Position::size() const { return distance_from(); }

// helper functions

Position operator+(const Position& p1, const Position& p2);
Position operator-(const Position& p1, const Position& p2);

typedef Position Point;
}
#endif
