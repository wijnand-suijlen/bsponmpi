#ifndef BSPONMPI_CONFIG_HPP
#define BSPONMPI_CONFIG_HPP

#include <string>
#include <sstream>
#include <cstdlib>

namespace bsplib {

    template <class T> bool read_env( const char * str, T & result ) {
#if HAS_DUPENV_S_WIN32
        char * env = NULL;
        if (_dupenv_s( &env, NULL , str) )
            return false;
#else
        const char * env = std::getenv( str );
#endif
        if (env) {
            std::istringstream s( env );
#if HAS_DUPENV_S_WIN32
            free(env);
#endif
            s >> result;
            if (!s)
                return true;
        }
        return false;
    }


}


#endif
