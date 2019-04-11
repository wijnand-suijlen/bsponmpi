#include "lincost.h"

#include <algorithm>
#include <cassert>
#include <limits>

#include "util.h"

namespace bsplib {

LinCost::LinCost( double alpha, double beta )
     : m_alpha( alpha )
     , m_beta( beta )
     , m_nprocs(0)
     , m_sizes()
     , m_counts()
{}

void LinCost::reset(int nprocs)
{
    m_nprocs = nprocs;
    m_sizes.clear(); m_counts.clear();
    m_sizes.reserve( nprocs + 1);
    m_counts.reserve( nprocs + 1);
    m_sizes.push_back(0); // we insert 0 as an extra option
}

void LinCost::send( std::size_t size )
{
    m_sizes.push_back( size );
}

std::size_t LinCost::get_bruck_vol()
{
    if ( m_nprocs == 1 ) return m_sizes.back();

    // Make sure the list m_sizes consists of unique numbers, sorted in
    // increasing order, while m_counts contains the multiplicity
    std::sort( m_sizes.begin(), m_sizes.end() );
    unsigned k = 0;
    for ( unsigned i = 0, j=0; i < m_sizes.size(); i=j ) {

        for (j=i+1 ; j < m_sizes.size(); ++j ) 
            if ( m_sizes[i] != m_sizes[j] )
                break;

        m_counts.push_back( j-i );
        m_sizes[k] = m_sizes[i];
        k++;
    }
    m_sizes.resize(k);
    assert( m_counts.size() == k );
    assert( m_sizes[0] == 0 );
    assert( m_counts[0] >= 1);
    m_counts[0] -= 1; // remove one occurrence of zero


    // now we can find the minimum cost
    double min_cost = std::numeric_limits<double>::max();
    std::size_t min_m = 0;
    std::size_t p2p_volume = 0, relation = 0;
    int logP = 1 + int_log2( m_nprocs-1 );
    for ( unsigned i = m_sizes.size(); i > 1; --i ) {
        std::size_t m = m_sizes[i-1];

        double cost = logP * ( m_alpha + m_beta * m_nprocs * m ) + 
            relation * m_alpha + p2p_volume * m_beta;

        if (cost < min_cost) {
            min_cost = cost;
            min_m = m;
        }

        relation += m_counts[i-1];
        std::size_t next_m = m_sizes[i-2];
        p2p_volume += relation * (m - next_m);
    }
    reset(m_nprocs);
    // now test for m = zero
    if ( relation * m_alpha + p2p_volume * m_beta < min_cost ) {
        return 0;
    }
    return min_m;
}

}
