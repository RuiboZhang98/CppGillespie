#include "TumorGenerator.h"

using std::cout,std::endl;

TumorGenerator::TumorGenerator(const Eigen::ArrayXd& param_initial_population,
                    const Eigen::ArrayXXd& param_transition_rates,
                    const double& param_tmax, const double& param_dt):
                    // collect parameters
                    initial_population(param_initial_population),
                    transition_rates(param_transition_rates),
                    tmax(param_tmax),
                    dt(param_dt)
{
    ntype = param_initial_population.size();
}

void TumorGenerator::print_status(){
    cout << " time = " << time << endl;
    cout << "current population = " << population.transpose() << endl;
}

void TumorGenerator::print_parameters(){
    // print tumor parameters
    cout << "initial population = " << initial_population.transpose() << endl; 
    cout << "transition rate matrix: \n" << transition_rates << endl;
    cout << "maximum time = " << tmax << endl;
    cout << "grid size of time (delta t) = " << dt << endl;
}

void TumorGenerator::print_thread(){
    cout << "Running on thread" << std::this_thread::get_id() << endl;
}

void TumorGenerator::initiation(){
        time = 0;
        population = initial_population;
}

//double TumorGenerator::lifespan(const Eigen::Ref<const Eigen::ArrayXXd>& weights){
double TumorGenerator::lifespan(const Eigen::ArrayXXd& weights){
    // The function computes the waiting time for next population change.
    double lambda = weights.sum();
    std::exponential_distribution<> exp{ lambda };
    return exp(rnd_generator);
} 

//template <typename Derived> //std::array<int,2>, don't need REF
std::tuple<int,int> TumorGenerator::increment_type(const Eigen::ArrayXXd &weights)
{
    // The function gives the type (a integer) that increases its population during 
    // population change.
    int ntype = weights.rows();
    std::discrete_distribution<> d(weights.reshaped().begin(),weights.reshaped().end());

    int random_entry = d(rnd_generator);
    int c = random_entry / ntype;
    int r = random_entry - c * ntype;

    int decrement = -1;
    if(r != c) decrement = r;

    return {decrement, c};
}

void TumorGenerator::evolve_step(){

    Eigen::ArrayXXd weights = transition_rates.colwise() * population;

    auto [decre_type, incre_type] = increment_type(weights); //edit 8/14

    population(incre_type)++;
    if (decre_type != -1) population(decre_type)--;
    time += lifespan(weights);
}

// Get a single realization of a tumor 
Eigen::ArrayXXd TumorGenerator::single_tumor(unsigned seed){
    
    // set the seed of the random number generator
    rnd_generator.seed(seed);
    const int datalen = (int)(tmax / dt) + 1;  // length of recorded data
    const Eigen::VectorXd record_time = Eigen::VectorXd::LinSpaced(datalen, 0, tmax); // record population at these populations
    Eigen::ArrayXXd time_population = Eigen::ArrayXXd::Constant(datalen, ntype + 1, -1.0); // initiate the data table

    time_population.block(0, 0, datalen, 1) = record_time; // flush the first column with the record times

    // initialize a tumor
    initiation();
    Eigen::ArrayXd population_old = population;

    // data_index indicates the current row of the result
    int data_index = 0;

    while(data_index < datalen){
        for ( ; (data_index < datalen) && (record_time(data_index) <= time) ; data_index++){
                time_population.block(data_index, 1, 1, ntype) = population_old.transpose();
        }
        population_old = population;
        evolve_step();
    }

    return time_population;
}