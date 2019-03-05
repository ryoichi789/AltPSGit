#ifndef CBA_REGION_H
#define CBA_REGION_H

#include <vector>
#include <utility>

namespace Cba {

typedef std::pair<int, int> Region;

class Regions {
public:
	void clear();
	void add(size_t first, size_t second);
		// region from 'first' to 'second' (NOTE: includes 'second')
	const Region& operator[](size_t i) const;
	int size() const;

private:
	std::vector<Region> m_regions;
};

inline void Regions::clear() { m_regions.clear(); }
inline void Regions::add(size_t first, size_t second)
	{ m_regions.push_back(Region(first, second)); }
inline int Regions::size() const { return m_regions.size(); }
inline const Region& Regions::operator[](size_t i) const { return m_regions[i]; }

}
#endif
