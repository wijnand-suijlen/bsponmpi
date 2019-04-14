#ifndef BSPONMPI_LINCOST_H
#define BSPONMPI_LINCOST_H

#include <vector>
#include <cstddef>

namespace bsplib {

/** Does communication optimisation on the basis of the assumption that sending
 * a point-to-point message costs T(n) = alpha + beta n, where alpha and beta
 * are parameters of the network. In particular, this class computes the amount
 * of data that should be sent through Bruck's algorithm while the rest goes
 * through point-to-point messages.
 *
 * Assuming a linear cost model, we search M that minimizes total communication cost
 *
 * \f[ (\alpha + \beta M P) \log P + \sum_{i=1}{P} \mathbf{1}\left[ X_i > M \right](\alpha + \beta(X_i - M)) \f]
 *
 * The left term comes from Bruck's algorithm, while the right term is the cost of 
 * sending the remaining messages through point-to-point communications.
 */
class LinCost {
public:
    LinCost( double alpha, double beta );

    void reset( int nprocs );

    void send( std::size_t size )
    { m_sizes.push_back( size ); }

    std::size_t get_bruck_vol() ;

private:
    double m_alpha, m_beta;
    int m_nprocs;
    std::vector< std::size_t > m_sizes, m_counts;
};


}

#endif
