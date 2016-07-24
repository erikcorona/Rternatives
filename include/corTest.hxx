//
// Created by Erik Corona on 7/23/16.
//

#ifndef RTERNATIVES_CORTEST_HXX
#define RTERNATIVES_CORTEST_HXX

#include <vector>

namespace fastR{

    class CorTest
    {
    public:
        CorTest(std::vector<double>& aa, std::vector<double>& bb)
        {
            a = aa;
            b = bb;
            n = a.size();
        }

        double n_0()
        {
            return n*(n-1)/2;
        }

        double tau_a()
        {
            int discordant{0}, concordant{0};
            for(int i = 0; i < n; i++)
                for(int k = i+1; k < n; k++)
                {
                    if (a[i] < a[k])
                    {
                        if (b[i] < b[k])
                            concordant++;
                        else if (b[i] > b[k])
                            discordant++;
                        else
                            std::cout << "problem1: " << a[i] << ", " << a[k] << ", " << b[i] << ", " << b[k] << std::endl;
                    }
                    else if (a[i] > a[k])
                    {
                        if (b[i] > b[k])
                            concordant++;
                        else if (b[i] < b[k])
                            discordant++;
                        else
                            std::cout << "problem2: " << a[i] << ", " << a[k] << ", " << b[i] << ", " << b[k] << std::endl;
                    }
                    else
                        std::cout << "problem3: " << a[i] << ", " << a[k] << ", " << b[i] << ", " << b[k] << std::endl;
                }

            return (concordant - discordant)/(n_0());
        }

    private:
        std::vector<double> a, b;
        double n;
    };
}
#endif //RTERNATIVES_CORTEST_HXX
