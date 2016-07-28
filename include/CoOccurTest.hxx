/* 
 * File:   CoOccurTest.hxx
 * Author: Erik Corona
 *
 */

#ifndef COOCCURTEST_HXX
#define	COOCCURTEST_HXX

#include <boost/math/special_functions/gamma.hpp>
#include "Hypothesis.hxx"

namespace fastR
{
constexpr long double percentChangeCutoff = 0.00000001;


class CoOccurTest
{
public:
    
    friend std::ostream& operator<< (std::ostream& out, CoOccurTest& test);
    
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
    CoOccurTest(const long double a, const long double b, const long double c, const long double d) :
            aa{a > d ? d : a},
            bb{b},
            cc{c},
            dd{a > d ? a : d},
            k{lgamma(a + b + 1) + lgamma(c + d + 1) + lgamma(a + c + 1) + lgamma(b + d + 1) - lgamma(a + b + c + d + 1)},
            a_d_switched{a > d}
    {
    }

    /**
     * Computes a p-value for seeing an observaton as extreme or more extreme
     * than this one.
     * @param pvalue_type specify whether the p-value should be computed with
     *                    respect to lower, greater, or two tailed side of the
     *                    distribution 
     * @return significance p-value
     */
    long double significance(alternative pvalue_type) const
    {
        if (pvalue_type == less)
            return less_impl();
        else if (pvalue_type == greater)
            return greater_impl();
        return two_tailed_impl();
    }

    /**
     * Computes the odds ratio
     * @return odds ratio
     */
    long double oddsRatio() const
    {
        return (aa/bb)/(cc/dd);
    }

private:

    /**
     * Computes a two-tailed test
     * @return sum of all p-values that are equal than or lower than what is
     * observed among all possible configurations with fixed margins
     */
    long double two_tailed_impl() const
    {
        long double sig{0}, ap;
        long a = aa, b = bb, c = cc, d = dd;
        auto pval = configuration_pvalue(a,b,c,d);
        
        do {
            ap = configuration_pvalue(a--, b++, c++, d--);
            if (ap <= pval) sig += ap;
        } while (ap / (sig - ap) > percentChangeCutoff && a >= 0 && d >= 0);

        a = aa+1, b = bb-1, c = cc-1, d = dd+1;
        
        do ap = configuration_pvalue(a++, b--, c--, d++);
        while (ap > pval && b >= 0 && c >= 0);
        
        sig += ap;
        if (b >= 0 && c >= 0)
        {
            do {
                ap = configuration_pvalue(a++, b--, c--, d++);
                if (ap < pval) sig += ap;
            } while (ap / (sig - ap) > percentChangeCutoff && b >= 0 && c >= 0);
        }
        
        return sig;
    }

    /**
     * Used to test for an odds ratio that is less than what is observed
     * @return sum of p-values representing an equal or smaller odds ratio
     * than observed
     */
    long double less_impl() const
    {
        long double sig{0}, ap;

        long a = aa, b = bb, c = cc, d = dd;

        do {
            ap = configuration_pvalue(a--, b++, c++, d--);
            sig += ap;
        } while (ap / (sig - ap) > percentChangeCutoff && a >= 0 && d >= 0);
        return sig;
    }

    /**
     * Used to test for an odds ratio greater than what is observed
     * @return sum of p-values representing an equal than or greater odds ratio
     * than observed
     */
    long double greater_impl() const
    {
        long double sig{0}, ap;

        long double a = aa, b = bb, c = cc, d = dd;

        do {
            ap = configuration_pvalue(a++, b--, c--, d++);
            sig += ap;
           } while (ap / (sig - ap) > percentChangeCutoff && b >= 0 && c >= 0);
        return sig;
    }
    
    std::string toString(std::string message = "") const
    {
        long double tmpa = a_d_switched ? dd : aa;
        long double tmpd = a_d_switched ? aa : dd;
        std::string ret = message;
        ret.append(message.size() > 0 ? ":\t" : "")
                .append(std::to_string(tmpa)).append(", ")
                .append(std::to_string(bb)).append(", ")
                .append(std::to_string(cc)).append(", ")
                .append(std::to_string(tmpd)).append(", ")
        .append("Configuration p-value: ").append(std::to_string(configuration_pvalue(aa,bb,cc,dd)))
                .append("\tLess Significance: ").append(std::to_string(significance(alternative::less)))
                .append("\tGreater Significance: ").append(std::to_string(significance(alternative::greater)))
                .append("\tTwo-Tailed Significance: ").append(std::to_string(significance(alternative::two_tailed)))
                .append("\tOR: ").append(std::to_string(oddsRatio()));
        ret.append("\n");
        return ret;
    }
    /**
    * Computes the probability of this exact table, not cumulative so it can't be used for significance
    **/
    long double configuration_pvalue(const long a, const long b, const long c, const long d) const
    {
        return std::exp(k - lgamma(a + 1) - lgamma(b + 1) - lgamma(c + 1) - lgamma(d + 1));
    }

    const long double aa, bb, cc, dd; // the counts in a 2x2 contingency table
    const long double k; // this is a constant used to compute p-values for each table configuration
    
    // a and d can be switched without affecting significance results and this
    // was done to avoid having to write special code when a is greater than d
    const bool a_d_switched; 
};

std::ostream& operator<< (std::ostream& out, CoOccurTest& test)
{
    out << test.toString();
    return out;
}

} // namespace fastR closing bracket
#endif	/* COOCCURRTEST_HXX */
