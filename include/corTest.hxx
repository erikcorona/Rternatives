//
// Created by Erik Corona on 7/23/16.
//

#ifndef RTERNATIVES_CORTEST_HXX
#define RTERNATIVES_CORTEST_HXX

#include <iostream>
#include <vector>
#include <unordered_map>
#include "utils.hxx"
#include "Hypothesis.hxx"

namespace fastR {
    class CorTest {
    public:
        CorTest(std::vector<double>& aa, std::vector<double>& bb) : n{aa.size()} {
            a = aa;
            b = bb;
        }

        double n_0() {
            return n * (n - 1) / 2.0;
        }

        int n_2()
        {
            return count_ties(b);
        }


        std::vector<int> tiesVector(std::vector<double>& values)
        {
            std::unordered_map<double, int> freq;
            for(int i = 0; i < n; i++)
                freq[values[i]]++;

            std::vector<int> ties;
            for(auto& cnts : freq)
                if(cnts.second > 1)
                    ties.push_back(cnts.second);
            return ties;
        }

        int count_ties(std::vector<double>& values)
        {
            std::unordered_map<double, int> freq;
            for(int i = 0; i < n; i++)
                freq[values[i]]++;

            int sum{0};
            for(auto& cnts : freq)
                sum += cnts.second * (cnts.second-1)/2;
            return sum;
        }

        int n_1()
        {
            return count_ties(a);
        }

        int concordant_minus_discordant()
        {
            int diff{0};
            for(int i = 0; i < n; i++)
                for(int k = i+1; k < n; k++)
                {
                    diff += ((a[i] < a[k] && b[i] < b[k]) || (a[i] > a[k] && b[i] > b[k]));
                    diff -= ((a[i] < a[k] && b[i] > b[k]) || (a[i] > a[k] && b[i] < b[k]));
                }

            return diff;
        }

        double concordant()
        {
            int concord{0};
            for(int i = 0; i < n; i++)
                for(int k = i+1; k < n; k++)
                    if ((a[i] < a[k] && b[i] < b[k]) || (a[i] > a[k] && b[i] > b[k]))
                        concord++;
            return concord;
        }

        double discordant()
        {
            int discord{0};
            for(int i = 0; i < n; i++)
                for(int k = i+1; k < n; k++)
                    if ((a[i] < a[k] && b[i] > b[k]) || (a[i] > a[k] && b[i] < b[k]))
                        discord++;
            return discord;
        }

        double tau_a()
        {
            return (concordant() - discordant())/(n_0());
        }

        double tau_b()
        {
            return concordant_minus_discordant()/
                    std::sqrt((n_0()-n_1())*(n_0()-n_2()));
        }

        double computeV2(std::vector<int>& tiesVector_a, std::vector<int>& tiesVector_b)
        {
            double v1{0};
            for(int v : tiesVector_a)
            {
                double inner = v*(v-1)*(v-2);

                for(int k : tiesVector_b)
                    v1 += inner*k*(k-1)*(k-2)/(9*n*(n-1)*(n-2));
            }
            return v1;
        }

        double computeV1(std::vector<int>& tiesVector_a, std::vector<int>& tiesVector_b)
        {
            double v1{0};
            for(int v : tiesVector_a)
            {
                double inner = v*(v-1);

                for(int k : tiesVector_b)
                    v1 += inner*k*(k-1)/(2*n*(n-1));
            }
            return v1;
        }

        double z_b()
        {
            std::vector<int> tiesA = tiesVector(a);
            std::vector<int> tiesB = tiesVector(b);
            double v1 = computeV1(tiesA,tiesB);
            double v2 = computeV2(tiesA,tiesB);

            double vt{0};
            for(int val : tiesA)
                vt += val*(val-1)*(2*val+5);

            double vu{0};
            for(int val : tiesVector(b))
                vu += val*(val-1)*(2*val+5);

            double v0 = n*(n-1)*(2*n+5);
            double v = (v0-vt-vu)/(18+v1+v2);
            return concordant_minus_discordant()/
                    (std::sqrt(v));
        }

        /**
         * Probability of tau_b
         * @param alt specify whether it's two-sided, one sided lower, or one sided greater
         * @return p-value associated with tau_b
         */
        double p_b(fastR::alternative alt)
        {
            double z = z_b();
            double lesser = 0.5 * (1 + erf(z / std::sqrt(2)));
            if(alt == greater)
                return 1 - lesser;
            else if(alt == less)
                return lesser;
            else
                return std::min(lesser, 1 - lesser)*2;
        }

        CorTest randomDraw()
        {
            std::vector<double> randA = numUtils::resample(a);
            std::vector<double> randB = numUtils::resample(b);
            CorTest perm(randA, randB);
            return perm;
        }

        double pValue(fastR::alternative alt)
        {
            if(alt != fastR::greater) {
                std::cerr << "Not Yet Implemented!";
                return -1;
            }

            double b = tau_b();
            affirm = 0;
            count  = 0;
            while(affirm < 500 && count < 1000000)
            {
                if (randomDraw().tau_b() > b)
                    affirm++;
                count++;
            }
            return (affirm+1.0)/(count+1.0);
        }

        CorTest(CorTest && m) : n{m.n}
        {
            a = std::move(m.a);
            b = std::move(m.b);
        }

        void showPvalueStats()
        {
            std::cout << affirm << "/" << count << std::endl;
        }

    private:
        std::vector<double> a, b;
        const unsigned long n;
        int affirm;
        int count;
    };
}
#endif //RTERNATIVES_CORTEST_HXX
