#ifndef BSPONMPI_CONFIG_HPP
#define BSPONMPI_CONFIG_HPP

#include <string>
#include <sstream>
#include <cstdlib>

namespace bsplib {

    template <class T> bool read_env( const char * str, T & result ) {
        const char * env = std::getenv( str );
        if (env) {
            std::istringstream s( env );
            s >> result;
            if (!s)
                return true;
        }
        return false;
    }


}


#endif
