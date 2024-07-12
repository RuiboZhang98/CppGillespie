#include <iomanip>
#include <iostream>
#include <map>
#include <random>

#include <Eigen/Dense> // path: G:\\code\\cfiles\\eigen-3.4.0
#include <Eigen/Core>

// function prototype
template <typename Derived>
int increment_type(const Eigen::ArrayBase<Derived> &weights, std::mt19937_64 &rng, int &decrement);
 
int main()
{
    Eigen::ArrayXXd weights {   
        {1, 2, 0},             
        {0, 2, 1},              
        {0, 0, 0}
    };
    
    //long long unsigned int seed = 1;
    std::random_device rd;
    std::mt19937_64 mt(rd());
    int decre_type = -1;

    int incre_type = increment_type(weights, mt, decre_type);
    if(decre_type == -1){
        std::cout << "Birth Event: type " << incre_type << std::endl;
    } else {
        std::cout << "Mutation Event: From type " << decre_type << " to type "<< incre_type << std::endl;
    }

    return 0;
}

template <typename Derived>
int increment_type(const Eigen::ArrayBase<Derived> &weights, std::mt19937_64 &rng, int &decrement)
{
    // The function gives the type (a integer) that increases in a population change.
    int ntype = weights.rows();
    std::discrete_distribution<> d(weights.reshaped().begin(),weights.reshaped().end());

    int random_entry = d(rng);
    int c = random_entry / ntype;
    int r = random_entry - c * ntype;

    std::cout << " rd = " << random_entry << std::endl;
    std::cout << " row: " << r << std::endl;
    std::cout << " col: " << c << std::endl;

    if(r != c) decrement = r;

    return c;
}