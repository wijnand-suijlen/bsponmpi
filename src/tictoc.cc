#include "tictoc.h"
#include <cstdio>
#include <iostream>
#include <iomanip>

namespace bsplib {
    std::vector< TicToc::Record > TicToc::s_timers;
    unsigned TicToc::s_pow = 1;
    unsigned TicToc::s_pos = 0;

    void TicToc::printStats(std::ostream & out)
    {
        const char * names[] = 
          { "DUMMY", "SYNC", "PUT", "GET", "HPPUT", "HPGET", "BSMP", 
            "MPI_META_A2A", "MPI_SMALL_A2A", "MPI_LARGE_RECV",
            "MPI_LARGE_SEND", "MPI_PUT", "MPI_UNBUF",
            "UNKNOWN"
          };

        out << "      usec        MB   MB/sec     Count\n";
        for (size_t i = 1; i < s_timers.size(); ++i ) {

            if (s_timers[i].count == 0 ) continue;
           
            double us = usec( s_timers[i].time );
            double mb = s_timers[i].nbytes * 1e-6;
            double mbs = s_timers[i].nbytes / us;
            double count = s_timers[i].count;

            out << std::setprecision(3)
                << std::setw(10) << us 
                << std::setw(10) << mb 
                << std::setw(10) << mbs
                << std::setw(10) << count;

            size_t pos = i;
            size_t pow = 1;
            do {
                unsigned j = (pos / pow) % N_CATEGORIES;
                pow *= N_CATEGORIES;

                out << ' ' << std::setw(16) << names[j];

            } while (pos > pow);
            out << '\n';
        }
    }


}
