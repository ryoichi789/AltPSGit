#include <cmath>
#include "position.h"
#include "exception.h"

using namespace Cba;

//
Position::Position() : m_x(0), m_y(0), m_z(0) {}

//
Position::Position(double x, double y, double z) : m_x(x), m_y(y), m_z(z) {}

//
Position::Position(const Position& p)
{
	m_x = p.m_x;
	m_y = p.m_y;
	m_z = p.m_z;
}

// operator=:
Position& Position::operator=(const Position& p)
{
	m_x = p.m_x;
	m_y = p.m_y;
	m_z = p.m_z;
	return *this;
}

// clear:
void Position::clear() { m_x = m_y = m_z = 0.0; }

// move_to:
Position& Position::move_to(double x, double y, double z)
{
	m_x = x;
	m_y = y;
	m_z = z;
	return *this;
}

// shift_by:
Position& Position::shift_by(double x, double y, double z)
{
	m_x += x;
	m_y += y;
	m_z += z;
	return *this;
}

// rotate_by:
Position& Position::rotate_by(const Matrix& rot)
{
	if (rot.nrows() != 3) return *this;
	if (rot.ncolumns() != 3) return *this;

	double x = rot.elem(0,0)*m_x + rot.elem(0,1)*m_y + rot.elem(0,2)*m_z;
	double y = rot.elem(1,0)*m_x + rot.elem(1,1)*m_y + rot.elem(1,2)*m_z;
	double z = rot.elem(2,0)*m_x + rot.elem(2,1)*m_y + rot.elem(2,2)*m_z;
	m_x = x;
	m_y = y;
	m_z = z;
	return *this;
}

// operator*=:
Position& Position::operator*=(double a)
{
	m_x *= a;
	m_y *= a;
	m_z *= a;
	return *this;
}

// operator/=:
Position& Position::operator/=(double a)
{
	if (a == 0.0) throw Exception("dividing by 0");
	m_x /= a;
	m_y /= a;
	m_z /= a;
	return *this;
}

// square_distance_from:
double Position::square_distance_from(double x, double y, double z) const
{
	double dx = m_x - x;
	double dy = m_y - y;
	double dz = m_z - z;
	return (dx*dx + dy*dy + dz*dz);
}

// is_within_range:
bool Position::is_within_range(const Position& pos, double range) const
{
	double dx = std::fabs(x() - pos.x());
	if (dx > range) return false;

	double dy = std::fabs(y() - pos.y());
	if (dy > range) return false;

	double dz = std::fabs(z() - pos.z());
	if (dz > range) return false;

	double dd = dx*dx + dy*dy + dz*dz;
	if (dd > range*range) return false;

	return true;
}

//
double Position::distance_from_neighbor(const Position& p, double range) const
{
	double d = square_distance_from_neighbor(p, range);
	if (d < 0.0) return -1.0;
	return std::sqrt(d);
}

//
double Position::square_distance_from_neighbor(const Position& p, double range) const
{
	double dx = std::fabs(x() - p.x());
	if (dx > range) return -1.0;

	double dy = std::fabs(y() - p.y());
	if (dy > range) return -1.0;

	double dz = std::fabs(z() - p.z());
	if (dz > range) return -1.0;

	double dd = dx*dx + dy*dy + dz*dz;
	if (dd > range*range) return -1.0;

	return dd;
}

// unit:
Position Position::unit() const
{
	Position p(*this);
	double s = size();
	if (s > 0) p /= s;
	return p;
}

// operator+:
Position Cba::operator+(const Position& p1, const Position& p2)
{
	Position p(p1);
	p.shift_by(p2);
	return p;
}

// operator-:
Position Cba::operator-(const Position& p1, const Position& p2)
{
	Position p(p1);
	p.shift_by(-p2.x(), -p2.y(), -p2.z());
	return p;
}
