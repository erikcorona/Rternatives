//
// Created by Erik Corona on 7/26/16.
//

#ifndef RTERNATIVES_UTILS_HXX
#define RTERNATIVES_UTILS_HXX

#include <vector>
#include <random>

namespace numUtils
{
    std::random_device r;

    std::vector<double> resample(std::vector<double> &values)
    {
        std::mt19937 eng{r()};

        std::vector<double> ret(values.size());
        std::uniform_int_distribution<long> dist(0, values.size() - 1);
        for (int i = 0; i < values.size(); i++)
            ret[i] = values[dist(eng)];
        return ret;
    }
}
#endif //RTERNATIVES_UTILS_HXX
