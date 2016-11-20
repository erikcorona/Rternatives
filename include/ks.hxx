//
// Created by Erik Corona on 11/18/16.
//

#ifndef KS_KS_HXX
#define KS_KS_HXX

#include "Hypothesis.hxx"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <unordered_set>


constexpr double M_1_SQRT_2PI{0.398942280401432677939946059934};

bool hasDuplicates(std::vector<double> &v) {
    std::unordered_set<double> map;
    for (double &num : v)
    {
        if(map.count(num))
            return true;
        map.insert(num);
    }
    return false;
}

/**
 * This implementation matches that found in R, including starting with 1 as the first position for the order function
 * @param array
 * @return
 */
std::vector<int> order(std::vector<double> &array)
{
    std::vector<std::pair<double,int>> nums(array.size());
    int count = 1;
    for(int i = 0; i < array.size();i++)
        nums[i] = std::make_pair(array[i], count++);

    std::sort(nums.begin(), nums.end(), [&](auto &a, auto &b){ return a.first < b.first; });

    std::vector<int> pos(array.size());

    for(int i = 0; i < array.size();i++)
        pos[i] = nums[i].second;
    return pos;
}

std::vector<double> cumsum(std::vector<double> &x)
{
    std::vector<double> sums(x.size());
    sums[0] = x[0];
    for (int i = 1 ; i < x.size() ; i++)
        sums[i] = x[i] + sums[i-1];
    return sums;
}

std::vector<double> diff(std::vector<double> &x)
{
    std::vector<double> diffs(x.size()-1);

    for(int i = 1; i < x.size(); i++)
        diffs[i-1] = x[i] - x[i-1];
    return diffs;
}

double maxAbsValue(std::vector<double> &x)
{
    double maxValue = std::abs(x[0]);

    for(int i = 1; i < x.size(); i++)
        if(std::abs(x[i]) > maxValue)
            maxValue = std::abs(x[i]);
    return maxValue;
}

double maxValue(std::vector<double> &x)
{
    double maxValue = x[0];

    for(int i = 1; i < x.size(); i++)
        if(x[i] > maxValue)
            maxValue = x[i];
    return maxValue;
}

double minValue(std::vector<double> &x)
{
    double minValue = x[0];

    for(int i = 1; i < x.size(); i++)
        if(x[i] < minValue)
            minValue = x[i];
    return minValue;
}


/* Two-sided two-sample */
static double psmirnov2x(double x, int m, int n)
{
    double md, nd, q, w;
    int i, j;

    if(m > n)
    {
        i = n;
        n = m;
        m = i;
    }
    md = (double) m;
    nd = (double) n;
    /*
       q has 0.5/mn added to ensure that rounding error doesn't
       turn an equality into an inequality, eg abs(1/2-4/5)>3/10
    */
    q = (0.5 + floor(x * md * nd - 1e-7)) / (md * nd);
    std::vector<double> u ((unsigned long) (n + 1));

    for(j = 0; j <= n; j++)
        u[j] = ((j / nd) > q) ? 0 : 1;

    for(i = 1; i <= m; i++) {
        w = (double)(i) / ((double)(i + n));
        if((i / md) > q)
            u[0] = 0;
        else
            u[0] = w * u[0];
        for(j = 1; j <= n; j++) {
            if(fabs(i / md - j / nd) > q)
                u[j] = 0;
            else
                u[j] = w * u[j] + u[j - 1];
        }
    }
    return u[n];
}


void pkstwo(int n, double *x, double tol)
{
/* x[1:n] is input and output
 *
 * Compute
 *   \sum_{k=-\infty}^\infty (-1)^k e^{-2 k^2 x^2}
 *   = 1 + 2 \sum_{k=1}^\infty (-1)^k e^{-2 k^2 x^2}
 *   = \frac{\sqrt{2\pi}}{x} \sum_{k=1}^\infty \exp(-(2k-1)^2\pi^2/(8x^2))
 *
 * See e.g. J. Durbin (1973), Distribution Theory for Tests Based on the
 * Sample Distribution Function.  SIAM.
 *
 * The 'standard' series expansion obviously cannot be used close to 0;
 * we use the alternative series for x < 1, and a rather crude estimate
 * of the series remainder term in this case, in particular using that
 * ue^(-lu^2) \le e^(-lu^2 + u) \le e^(-(l-1)u^2 - u^2+u) \le e^(-(l-1))
 * provided that u and l are >= 1.
 *
 * (But note that for reasonable tolerances, one could simply take 0 as
 * the value for x < 0.2, and use the standard expansion otherwise.)
 *
 */
    double new2, old, s, w, z;
    int i, k, k_max;

    k_max = (int) sqrt(2 - log(tol));

    for(i = 0; i < n; i++) {
        if(x[i] < 1) {
            z = - (M_PI_2 * M_PI_4) / (x[i] * x[i]);
            w = log(x[i]);
            s = 0;
            for(k = 1; k < k_max; k += 2) {
                s += exp(k * k * z - w);
            }
            x[i] = s / M_1_SQRT_2PI;
        }
        else {
            z = -2 * x[i] * x[i];
            s = -1;
            k = 1;
            old = 0;
            new2 = 1;
            while(fabs(old - new2) > tol) {
                old = new2;
                new2 += 2 * s * exp(z * k * k);
                s *= -1;
                k++;
            }
            x[i] = new2;
        }
    }
}

double pKS2(double statistic, double stol)
{
    double tol = stol;
    double ans = statistic;
    pkstwo(1, &ans, tol);
    return ans;
}

double pkstwo (double x, double tol = 1e-06)
{
    return x > 0 ? pKS2(x,tol) : 0;
}

struct KSVal
{
    double statistic, pValue;
    fastR::alternative method;
};

//                                                          c("two_sided", "less", "greater")
// exact has to be set explicitly, not automated good rule is exact = x.size() * y.size() < 10000
// two.sided = "two-sided", less = "the CDF of x lies below that of y", greater = "the CDF of x lies above that of y"
template<typename ITER>
KSVal kstest2sample2 (ITER xbegin, ITER xend, ITER ybegin, ITER yend, const fastR::alternative method, const bool exact)
{
    double nx = std::distance(xbegin, xend);
    double ny = std::distance(ybegin, yend);

    std::vector<double> w(xbegin, xend);
    w.reserve(nx + ny);
    w.insert(w.end(), ybegin, yend);

    std::vector<double> z(w.size());

    auto ranks = order(w);

    for(int i = 0; i < ranks.size(); i++)
        z[i] = ranks[i] <= nx ? 1 / nx : -1 / ny;
    z = cumsum(z);

    bool TIES = false;
    if (hasDuplicates(w))
    {
        if (exact) {
            std::cerr <<  "cannot compute exact p-value with ties" << std::endl;
            exact = false;
        }
        else
            std::cerr << "p-value will be approximate in the presence of ties" << std::endl;

        std::sort(w.begin(), w.end());

        auto aDiff = diff(w);

        std::vector<double> newZ(0);
        newZ.reserve(aDiff.size());
        for(int i = 0; i < aDiff.size(); i++)
            if(aDiff[i] != 0)
                newZ.push_back(z[i]);
        newZ.push_back(z[nx+ny-1]);

        z = std::move(newZ);
        TIES = true;
    }

    double statistic;
    switch(method)
    {
        case fastR::two_tailed : statistic =  maxAbsValue(z); break;
        case fastR::greater    : statistic =  maxValue   (z); break;
        case fastR::less       : statistic = -minValue   (z); break;
    }

    double PVAL = std::numeric_limits<double>::quiet_NaN();
    if (exact && method == fastR::two_tailed && !TIES)
        PVAL = 1 - psmirnov2x(statistic, (int) nx, (int)ny);

    double n = (nx * ny)/(nx + ny);
    if (std::isnan(PVAL))
        PVAL = method == fastR::two_tailed ? 1 - pkstwo(sqrt(n) * statistic) : exp(-2 * n * statistic*statistic);

    PVAL = std::min(std::max(PVAL,0.0),1.0);

    KSVal RVAL = {statistic, PVAL, method};

    return RVAL;
}




#endif //KS_KS_HXX
