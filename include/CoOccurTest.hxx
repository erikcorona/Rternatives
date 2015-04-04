/* 
 * File:   CoOccurTest.hxx
 * Author: Erik Corona
 *
 */

#ifndef COOCCURTEST_HXX
#define	COOCCURTEST_HXX

#include <algorithm>
#include <boost/math/special_functions/gamma.hpp>
#include <limits>
#include <map>
#include <math.h>
#include <vector>

using namespace boost::math;

namespace fastR{
class CoOccurTest
{
public:
    /**
     * number of times category 1 occurs is a+c
     * number of times category 2 occurs is a+b
     * number of times category 1 doesn't occur is b+d
     * number of times category 2 doesn't occur is c+d
     * total number of times anything occurs is a+b+c+d
     * 
     *     | c1 |!c1 |
     *-----|----|----|--------
     *  c2 | a  | b  | a+b
     * ----|----|----|--------
     * !c2 | c  | d  | c+d
     *-----|----|----|--------
     *     |a+c |b+d | a+b+c+d
     * 
     * @param a number of times category 1 and category b co-occur
     * @param b number of times category 2 co-occurs with anything that isn't category 1
     * @param c number of times category 1 co-occurs with anything that isn't category 2
     * @param d number of times neither category 1 nor category 2 occur
     */
    CoOccurTest(const long a, const long b, const long c, const long d) :
            aa{a > d ? d : a},
            bb{b},
            cc{c},
            dd{a > d ? a : d},
            k{lgamma((double)a + b + 1) + lgamma((double)c + d + 1) + lgamma((double)a + c + 1) + lgamma((double)b + d + 1) - lgamma((double)a + b + c + d + 1)},
            a_d_switched{a > d}
    {
    }

    /**
     * Computes a two-tailed test
     * @return sum of all p-values that are equal than or lower than what is
     * observed among all possible configurations with fixed margins
     */
    long double two_tailed() const
    {
        long double sig{0}, ap;
        long a = aa, b = bb, c = cc, d = dd;
        auto pval = configuration_pvalue(a,b,c,d);
        
        do {
            ap = configuration_pvalue(a--, b++, c++, d--);
            if (ap <= pval) sig += ap;
        } while (ap / (sig - ap) > 0.00000001 && a >= 0 && d >= 0);

        a = aa+1, b = bb-1, c = cc-1, d = dd+1;
        
        do ap = configuration_pvalue(a++, b--, c--, d++);
        while (ap > pval && b >= 0 && c >= 0);
        
        sig += ap;
        if (b >= 0 && c >= 0)
        {
            do {
                ap = configuration_pvalue(a++, b--, c--, d++);
                if (ap < pval) sig += ap;
            } while (ap / (sig - ap) > 0.00000001 && b >= 0 && c >= 0);
        }
        
        return sig;
    }

    /**
     * Used to test for an odds ratio that is less than what is observed
     * @return sum of p-values representing an equal or smaller odds ratio
     * than observed
     */
    long double less() const
    {
        long double sig{0}, ap;

        long a = aa, b = bb, c = cc, d = dd;

        do {
            ap = configuration_pvalue(a--, b++, c++, d--);
            sig += ap;
        } while (ap / (sig - ap) > 0.00000001 && a >= 0 && d >= 0);
        return sig;
    }

    /**
     * Used to test for an odds ratio greater than what is observed
     * @return sum of p-values representing an equal than or greater odds ratio
     * than observed
     */
    long double greater() const
    {
        long double sig{0}, ap;

        long a = aa, b = bb, c = cc, d = dd;

        do {
            ap = configuration_pvalue(a++, b--, c--, d++);
            sig += ap;
           } while (ap / (sig - ap) > 0.00000001 && b >= 0 && c >= 0);
        return sig;
    }

    /**
     * Computes the odds ratio
     * @return odds ratio
     */
    double oddsRatio() const
    {
        return ((double)aa/bb)/((double)cc/dd);
    }

    std::string toString(std::string message = "") const
    {
        long tmpa = a_d_switched ? dd : aa;
        long tmpd = a_d_switched ? aa : dd;
        std::string ret = message;
        ret.append(message.size() > 0 ? ":\t" : "")
                .append(std::to_string(tmpa)).append(", ")
                .append(std::to_string(bb)).append(", ")
                .append(std::to_string(cc)).append(", ")
                .append(std::to_string(tmpd)).append(", ")
        .append("Configuration p-value: ").append(std::to_string(configuration_pvalue(aa,bb,cc,dd)))
                .append("\tLess Significance: ").append(std::to_string(less()))
                .append("\tGreater Significance: ").append(std::to_string(greater()))
                .append("\tTwo-Tailed Significance: ").append(std::to_string(two_tailed()))
                .append("\tOR: ").append(std::to_string(oddsRatio()));
        ret.append("\n");
        return ret;
    }

private:

    /**
    * Computes the probability of this exact table, not cumulative so it can't be used for significance
    **/
    long double configuration_pvalue(const long a, const long b, const long c, const long d) const
    {
        return std::exp(k - lgamma((double)a + 1) - lgamma((double)b + 1) - lgamma((double)c + 1) - lgamma((double)d + 1));
    }

    const long aa, bb, cc, dd; // the counts in a 2x2 contingency table
    const long double k; // this is a constant used to compute p-values for each table configuration
    
    // a and d can be switched without affecting significance results and this
    // was done to avoid having to write special code when a is greater than d
    const bool a_d_switched; 
};
}
#endif	/* COOCCURRTEST_HXX */
