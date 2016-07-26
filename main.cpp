#include <iostream>
#include <CoOccurTest.hxx>
#include <CorTest.hxx>

using namespace std;

int main()
{
    fastR::CoOccurTest co(30,30,10,20);

    std::vector<double> a = {  95,3,6,77,23,1,335,13,37, 18, 29};
    std::vector<double> b = {952,13,2,37,33,11,332,82,7, 8, 9};

    std::cout << a.size() << std::endl;
    fastR::CorTest cor(a,b);
    cout << "co: " << co.oddsRatio() << endl;
    cout << "cor: " << cor.tau_a() << endl;

    std::vector<double> aa = {23,23,27,27,39,41,45,49,50,53,53,54,56,57,58,58,60,61};
    std::vector<double> bb = {9.5,27.9,7.8,17.8,31.4,25.9,27.4,25.2,31.1,34.7,42.0,29.1,32.5,30.3,33.0,33.8,41.1,34.5};

    fastR::CorTest cor2(aa,bb);
    cout << "cor2a: " << cor2.tau_a() << endl;
    cout << "cor2b: " << cor2.tau_b() << endl;
    cout << "cor2b-p: " << cor2.pValue(fastR::greater) << endl;
    return 0;
}