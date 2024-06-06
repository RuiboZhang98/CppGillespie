#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <vector>

class Gene { // The Gene class
    public:
    void print_gene_info(){
        std::cout << "order:" << order << "name:" << name << std::endl;
    }

    protected:
    int order;
    std::string name;
};

class Ocogene: public Gene {
    public: 
    Ocogene(const std::string& filePathAndName);
    void print_gene_info();

    private:
    std::string filePathAndName;
    float mutation_rate0;
    float growth_adv0;   
};

class TumorSuppressor: public Gene {
    public: 
    TumorSuppressor(const std::string& filePathAndName);
    void print_gene_info();

    private:
    std::string filePathAndName;
    float mutation_rate0;
    float mutation_rate1;
    float growth_adv0;
    float growth_adv1;
};

using std::fstream;
using std::ifstream;
using std::string;

using std::cout;
using std::endl;


Ocogene::Ocogene(const string& filePathAndName){
    ifstream file(filePathAndName);
    string str;

    // order
    std::getline(file, str);
    order = std::stoi(str);

    // name
    std::getline(file, str);
    name = str;

    // mutation rate 0 
    std::getline(file, str);
    mutation_rate0 = std::stof(str);

    // fitness advantage
    std::getline(file, str);
    growth_adv0 = std::stof(str);

}

TumorSuppressor::TumorSuppressor(const string& filePathAndName){
    ifstream file(filePathAndName);
    string str;

    // order
    std::getline(file, str);
    order = std::stoi(str);

    // name
    std::getline(file, str);
    name = str;

    // mutation rate 0 
    std::getline(file, str);
    mutation_rate0 = std::stof(str);

    // mutation rate 1
    std::getline(file, str);
    mutation_rate1 = std::stof(str);

    // fitness advantage 0
    std::getline(file, str);
    growth_adv0 = std::stof(str);

    // fitness advantage 1
    std::getline(file, str);
    growth_adv1 = std::stof(str);

}

void Ocogene::print_gene_info(){
    cout << "order: " << order << " ocogene:" << name << endl
         << "mutation rate 0:" << mutation_rate0 << endl
         << "fitness advantage:" << growth_adv0 << endl;
    cout << endl;
}

void TumorSuppressor::print_gene_info(){
    cout << "order: " << order << " TSG:" << name << endl
         << "mutation rate 0:" << mutation_rate0 << endl
         << "mutation rate 1:" << mutation_rate1 << endl
         << "fitness advantage 0 :" << growth_adv0 << endl
         << "fitness advantage 1 :" << growth_adv1 << endl;
    cout << endl;
}