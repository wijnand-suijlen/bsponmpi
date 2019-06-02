#include "tictoc.h"
#include "mpisizet.h"
#include <cstdio>
#include <iostream>
#include <iomanip>

namespace bsplib {
    std::vector< TicToc::Record > TicToc::s_timers;
    unsigned TicToc::s_pow = 1;
    unsigned TicToc::s_pos = 0;



    void TicToc::print_stats(MPI_Comm comm, std::ostream & out)
    {
        const char * names[] = 
          { "DUMMY", "TOTAL", "SYNC", "PUT", "GET", "HPPUT", "HPGET", 
            "BSMP", "MPI_META_A2A", "MPI_SMALL_A2A", "MPI_LARGE_RECV",
            "MPI_LARGE_SEND", "MPI_PUT", "MPI_UNBUF",
            "IMBALANCE",
            "UNKNOWN"
          };
        const char * op[] = { "Min", "Max", "Avg" };
        enum { Min, Max, Avg };
        const int s=0, mb=1, mbs=2, count=3;

        int P; MPI_Comm_size( comm, &P);
        std::vector< double > allrecords(4*P), record(4);

        out << "        seconds PID        MB PID    MB/sec PID     Count PID\n";
        out << std::setprecision(3);
        for (size_t i = 1; i < s_timers.size(); ++i ) {

            std::size_t local_count = s_timers[i].count, max_count = 0;
            MPI_Allreduce( &local_count, &max_count, 1, MY_MPI_SIZE_T, 
                    MPI_MAX, comm );
            if (max_count == 0 ) continue;

            record[s]     = sec( s_timers[i].time );
            record[mb]    = s_timers[i].nbytes * 1e-6;
            record[mbs]   = record[mb] / record[s];
            record[count] = 1.0 * s_timers[i].count;

            MPI_Allgather( record.data(), 4, MPI_DOUBLE,
                    allrecords.data(), 4, MPI_DOUBLE,
                    comm );

            for (int j = 0; j < 3; ++j ) {
                out << op[j] << "  ";

                for (int k = 0; k < 4; ++k ) {
                    double x = allrecords[k];
                    int pid  = 0;
                    for (int p = 1; p < P; ++p) {
                        double y = allrecords[p*4+k];
                        switch( j )  {
                            case Min: if ( y < x ) { x = y; pid = p; }
                                      break;
                            case Max: if ( y > x ) { x = y; pid = p; }
                                      break;
                            case Avg: x += y;
                                      break;
                        }
                    }
                    if (j == Avg) { 
                        out << std::setw(10) << (x/P) << "    ";
                    } else {
                        out << std::setw(10) << x << std::setw(4) << pid;
                    }

                }

                size_t pos = i;
                size_t pow = 1;
                do {
                    unsigned j = (pos / pow) % N_CATEGORIES;
                    pow *= N_CATEGORIES;

                    out << ' ' << std::setw(16) << names[j];

                } while (pos > pow);
                out << '\n';
            }
            out << '\n';    
        }
    }


}
