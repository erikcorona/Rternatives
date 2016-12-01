//
// Created by Erik Corona on 11/24/16.
//

#include <corTest.hxx>
#include "CommonOperations.hxx"
#include "gtest/gtest.h"


TEST(MESSAGES, test_mean)
{
    std::vector<double> aaa = {0.92062, 0.92062, 0.92062, 0.92062, 0.98921, 0.98921, 0.98921};
    std::vector<double> bbb = {0, 0, 0, 0, 0, 0, 1};

    fastR::CorTest cor(aaa,bbb);

    EXPECT_NEAR(cor.tau_b(),0.4714045,0.000001);
    EXPECT_NEAR(cor.z_b(),1.154701,0.00000000000000001); //@todo fix this. Do not use z_b()

    auto rranks = fastR::rank(aaa, fastR::TiesMethod::average);

    for(auto &b : rranks)
        std::cout << b << ", ";
    std::cout << std::endl;

    std::cout.flush();
}
