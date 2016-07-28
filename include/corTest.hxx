//
// Created by Erik Corona on 7/23/16.
//

#ifndef RTERNATIVES_CORTEST_HXX
#define RTERNATIVES_CORTEST_HXX

#include <vector>
#include "utils.hxx"

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

        double n_2()
        {
            int ties{0};
            for(int i = 0; i < n; i++)
                for(int k = i+1; k < n; k++)
                    if (b[i] == b[k])
                        ties++;
            return ties;
        }

        double n_1()
        {
            int ties{0};
            for(int i = 0; i < n; i++)
                for(int k = i+1; k < n; k++)
                    if (a[i] == a[k])
                        ties++;
            return ties;
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
            return (concordant() - discordant())/
                    std::sqrt((n_0()-n_1())*(n_0()-n_2()));
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
            int affirm{0};
            int count{0};
            while(affirm < 500 && count < 1000000)
            {
                if (randomDraw().tau_b() > b)
                    affirm++;
                count++;
            }
            std::cout << affirm << " / " << count << std::endl;
            return (affirm+1.0)/(count+1.0);
        }

        CorTest(CorTest && m) : n{m.n}
        {
            a = std::move(m.a);
            b = std::move(m.b);
        }
    private:
        std::vector<double> a, b;
        const unsigned long n;
    };
}
#endif //RTERNATIVES_CORTEST_HXX
