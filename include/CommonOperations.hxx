//
// Created by Erik Corona on 11/24/16.
//

#ifndef RTERNATIVES_COMMONOPERATIONS_HXX
#define RTERNATIVES_COMMONOPERATIONS_HXX

#include <vector>
#include <iostream>
#include "Hypothesis.hxx"
#include <random>

namespace fastR
{





#define sort2_with_index \
    for (h = sincs[t]; t < 16; h = sincs[++t]) { \
    for (i = lo + h; i <= hi; i++) {     \
        itmp = indx[i];             \
        j = i;                             \
        while (j >= lo + h && less(indx[j - h], itmp)) {         \
        indx[j] = indx[j - h]; j -= h; }             \
        indx[j] = itmp;                         \
    }                                 \
    }

    constexpr int sincs[17] = {
            1073790977, 268460033, 67121153, 16783361, 4197377, 1050113,
            262913, 65921, 16577, 4193, 1073, 281, 77, 23, 8, 1, 0
    };

    /* Needs indx set to  0:(n-1)  initially.
   Also used by do_options and  ../gnuwin32/extra.c
   Called with rho != R_NilValue only from do_rank, when NAs are not involved.
 */
    void orderVector1(std::vector<int> &indx, unsigned long n, std::vector<double> &key, bool decreasing) {
        int i, j, h, t, lo = 0, hi = n - 1;
        int itmp;

        if (n < 2)
            return;

        std::vector<double> x = key;

        /* Shell sort isn't stable, so add test on index */

        for (t = 0; sincs[t] > hi - lo + 1; t++);

        if (decreasing) {
#define less(a, b) (x[a] < x[b] || (x[a] == x[b] && a > b))
            sort2_with_index
#undef less
        } else {
#define less(a, b) (x[a] > x[b] || (x[a] == x[b] && a > b))
            sort2_with_index
#undef less
        }
    }


    std::random_device rd;
    std::mt19937 gen(rd());
    /**
     * sample takes a sample of the specified size from the elements of x using either with or without replacement.
     * @param x  A vector of one or more elements from which to choose
     * @param size a non-negative integer giving the number of items to choose
     * @param replace Should sampling be with replacement?, currently it's with replacement @todo add replacement option
     */
    std::vector<double> sample (std::vector<double> &x, unsigned long size)
    {
        std::uniform_int_distribution<unsigned long> dis(0, x.size()-1);
        std::vector<double> res;
        res.resize(size);

        for(int i = 0; i < size; i++)
        {
            unsigned long index = dis(gen);
            res[i] = x[index];
        }

        return res;
    }


    int equal(int i, int j, std::vector<double> &x)
    {
        return x[i] == x[j];
    }

    /* FUNCTION: rank(x, length, ties.method) */
    auto do_rank(std::vector<double> &x, unsigned long n, const fastR::TiesMethod ties_kind) {
        std::vector<double> rank(n);

        if (n > 0) {
            std::vector<int> in(n);
            for (int i = 0; i < n; i++)
                in[i] = i;
            orderVector1(in, n, x, false);
            int j;
            for (int i = 0; i < n; i = j + 1) {
                j = i;

//                std::cout << in[j] << " " << in[j+1] << " " << fastR::equal(in[j], in[j + 1], x) << std::endl;
                while ((j < n - 1) && fastR::equal(in[j], in[j + 1], x))
                    j++;
                switch (ties_kind) {
                    case average:
                        for (int k = i; k <= j; k++)
                            rank[in[k]] = (i + j + 2) / 2.;
                        break;
                    case max:
                        for (int k = i; k <= j; k++)
                            rank[in[k]] = j + 1;
                        break;
                    case min:
                        for (int k = i; k <= j; k++)
                            rank[in[k]] = i + 1;
                        break;
                    default:
                        std::cerr << "ties_kind Not Supported" << std::endl;
                }
            }
        }
        return rank;
    }

    // tiesMethod = c("average", "first","random")
    auto rank(std::vector<double> &rank, fastR::TiesMethod tiesMethod)
    {
        switch (tiesMethod)
        {
            case fastR::average:
            case fastR::min    :
            case fastR::max    :
                return do_rank(rank, rank.size(), tiesMethod);
            case fastR::first  :
                std::cerr << "Not Yet Implemented" << std::endl;
//                return sort.list(sort.list(rank));
            case fastR::random :
                std::cerr << "Not Yet Implemented" << std::endl;
//                return sort.list(order(rank,runif(rank.size())));
            default:
                return std::vector<double>();
        }

    }
}

#endif //RTERNATIVES_COMMONOPERATIONS_HXX
