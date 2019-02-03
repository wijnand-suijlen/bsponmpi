#ifndef BSPONMPI_TICTOC_HPP
#define BSPONMPI_TICTOC_HPP

#include <mpi.h>

#if HAS_CLOCK_GETTIME
  #include <time.h>
#endif

#include <iosfwd>
#include <cassert>
#include <vector>


namespace bsplib {


class TicToc
{
public:
#if HAS_CLOCK_GETTIME
    typedef struct timespec Time;
#else
    typedef double Time;
#endif

    enum Category { DUMMY, TOTAL, SYNC, PUT, GET, HPPUT, HPGET, 
        BSMP, MPI_META_A2A, MPI_SMALL_A2A, MPI_LARGE_RECV,
        MPI_LARGE_SEND, MPI_PUT, MPI_UNBUF,
        IMBALANCE,
        N_CATEGORIES };

    TicToc( Category c, size_t bytes = 0u)
       : m_record( push(c) )
       , m_start( getTime() )
       , m_nbytes( bytes )
     {
     }

    void add_bytes( size_t bytes ) 
    { m_nbytes += bytes; }

    ~TicToc()
    {
        Record & r = s_timers[m_record];
        r.time = timeAdd( r.time, 
                timeDiff( getTime(), m_start ) );
        r.nbytes += m_nbytes;
        r.count += 1;
        pop();
    }


    static void print_stats(MPI_Comm comm, std::ostream & out);

private:
#if HAS_CLOCK_GETTIME
    static Time getTime() 
    { Time t; clock_gettime(CLOCK_MONOTONIC, &t); return t; }

    static Time timeAdd( Time a, Time b )
    { a.tv_sec += b.tv_sec; 
      a.tv_nsec += b.tv_nsec;
      if ( a.tv_nsec >= 1000000000l ) {
          a.tv_sec += 1;
          a.tv_nsec -= 1000000000l;
      } else if ( a.tv_nsec < 0) {
          a.tv_sec -= 1;
          a.tv_nsec += 1000000000l;
      }
      return a;
    }

    static Time timeDiff( Time a, Time b )
    { a.tv_sec -= b.tv_sec; 
      a.tv_nsec -= b.tv_nsec;
      if ( a.tv_nsec >= 1000000000l ) {
          a.tv_sec += 1;
          a.tv_nsec -= 1000000000l;
      } else if ( a.tv_nsec < 0) {
          a.tv_sec -= 1;
          a.tv_nsec += 1000000000l;
      }
      return a;
    }

    static Time zero() 
    { Time t; t.tv_sec = 0; t.tv_nsec = 0; return t;}

    static double sec( Time x ) { 
        return x.tv_sec + 1e-9 * x.tv_nsec;
    }

#else
    static Time getTime()
    { return MPI_Wtime(); }

    static Time timeAdd( Time a, Time b) { return a + b; }
    static Time timeDiff( Time a, Time b ) { return a - b; }
    static Time zero() { return 0.0; }
    static double sec( Time x ) { return x; }
#endif
    struct Record {
        Time time;
        size_t nbytes;
        size_t count;

        Record() : time(zero()), nbytes(0), count(0) {}
    };

   
    static unsigned push(Category c ) {
        s_pos += s_pow * c;
        s_pow *= N_CATEGORIES;
        assert( s_pow < 1000000 );

        if (s_timers.size() <= s_pos)
            s_timers.resize( s_pos + 1);

        return s_pos;
    }

    static void pop() {
        s_pow /= N_CATEGORIES;
        assert( s_pow > 0 );
        s_pos %= s_pow;
    }

    static std::vector< Record > s_timers;
    static unsigned s_pow;
    static unsigned s_pos;

    unsigned m_record;
    Time m_start;
    size_t m_nbytes;
};



}


#endif
