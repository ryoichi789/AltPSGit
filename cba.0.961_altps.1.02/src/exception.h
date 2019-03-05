#ifndef CBA_EXCEPTION_H
#define CBA_EXCEPTION_H

#include <iostream>
#include <string>

namespace Cba {

class Exception {
public:
	Exception() {}
	Exception(const char* m);
	void set_message(const char* m);
	void print_message() const;
	const char* get_message() const;

private:
	std::string m_message;
};

inline Exception::Exception(const char* m) { m_message = m; }
inline void Exception::set_message(const char* m) { m_message = m; }
inline void Exception::print_message() const { std::cerr << m_message << std::endl; }
inline const char* Exception::get_message() const { return m_message.c_str(); }

}

#endif
