// F. Borando - francesco.borando@studenti.unimi.it
// M. Garbellini - matteo.garbellini@studenti.unimi.it
// Università degli Studi di Milano
// Dept. of Physics
// January 2019
//
// GRAVITATIONAL N-BODY SIMULATION
// tool to analyze the dispersion of a distribution of energy data

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

double get_dispersion(const double* v_energy, const int& m_num_of_bodies);

int
main(int argc, char* argv[])
{
    if(!argv[1])
    {
        std::cerr << "Please insert filename as argv[1]!" << std::endl;
        std::exit (EXIT_FAILURE);
    }
    
    int index;
    int num_of_bodies=0;
    std::ifstream input;
    
    //find num_of_bodies
    input.open(argv[1]);
    if(!input)
    {
        std::cerr << "Error while opening the input file!" << std::endl;
        std::exit (EXIT_FAILURE);
    }
    
    while (input.good())
    {
        input >> index >> index;
        num_of_bodies++;
    }
    
    input.close();
    
    double energy[num_of_bodies];
    
    //find dispersion
    input.open(argv[1]);
    if(!input)
    {
        std::cerr << "Error while opening the input file!" << std::endl;
        std::exit (EXIT_FAILURE);
    }
    
    while (input.good())
    {
        input >> index >> energy[index];
    }
    
    input.close();
    
    std::cout << "Dispersion: " << get_dispersion(energy, num_of_bodies) << std::endl;
    
    return 0;
}

double get_dispersion(const double* v_energy, const int& m_num_of_bodies)
{
    double initial_energy = v_energy[0];
    long double dispersion=0;
    for (int i=0;i<m_num_of_bodies;i++)
    {
        dispersion += std::abs((initial_energy - v_energy[i]));
    }
    return (double) (dispersion/m_num_of_bodies);
}
