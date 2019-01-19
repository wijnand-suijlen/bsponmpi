#ifndef BSPONMPI_EXCEPTION_H
#define BSPONMPI_EXCEPTION_H

#include <stdexcept>
#include <sstream>
#include <string>

namespace bsplib {

struct exception : std::exception{ 
public:
    explicit exception(const char * name) 
      : m_stream()
      , m_buf()
    { m_stream << name ; }

    exception( const exception & e )
      : m_stream( e.m_stream.str() )
      , m_buf( )
    {}

    template <typename T>
    exception & operator<<( const T & val )
    { m_stream << val; return *this; }

    virtual const char * what() const throw()
    { m_buf = m_stream.str();
      return m_buf.c_str(); 
    }

    std::string str() const { return m_stream.str(); }

private:
    std::ostringstream m_stream;
    mutable std::string m_buf;
};



}

#endif
