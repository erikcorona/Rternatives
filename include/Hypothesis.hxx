//
// Created by Erik Corona on 7/28/16.
//

#ifndef RTERNATIVES_HYPOTHESIS_HXX
#define RTERNATIVES_HYPOTHESIS_HXX

namespace fastR
{
    /**
     * Specifies whether seeking probability of seeing lower or equal co-occurrence
     * counts
     */
    enum alternative { less, greater, two_tailed };
    enum TiesMethod {average,first,random,max,min};
}
#endif //RTERNATIVES_HYPOTHESIS_HXX
