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
    return 0;
}