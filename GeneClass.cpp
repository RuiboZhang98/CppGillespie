#include <iostream>
#include <fstream>
#include <string>

class Gene { // The Gene class
    public:
        int order;
        std::string name;
        void print_gene_info(){
            std::cout << "name:";
        }
};

class Ocogene: public Gene {
    public: 
    Ocogene(const std::string& filePathAndName);

    std::string filePathAndName;

    float mutation_rate0;
    float growth_adv0;   

    void print_gene_info();
};

using std::fstream;
using std::ifstream;
using std::string;

using std::cout;		// for debug mode
using std::endl;		// for debug mode


Ocogene::Ocogene(const std::string& filePathAndName){
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

//print_gene_info::TumorSuppressor(){
//    cout << "mutation rate 0:" << mutation_rate0
//         << "mutation rate 1:" << mutation_rate1
//}

void Ocogene::print_gene_info(){
    std::cout << "mutation rate 0:" << mutation_rate0 << std::endl
         << "fitness advantage:" << growth_adv0;
}

int main(){
    Ocogene Kras("Kras.txt");
    Kras.print_gene_info();
    return 0;
}