// F. Borando - francesco.borando@studenti.unimi.it
// M. Garbellini - matteo.garbellini@studenti.unimi.it
// Università degli Studi di Milano
// Dept. of Physics
// January 2019
//
// GRAVITATIONAL N-BODY SIMULATION
//transforms the result of a computation of the n-body problem of type
//x_position  y_position  z_position kin_energy pot_energy
//into new files useful for energetic analysis:
//• clean energy data for plotting
//• a "table" of percentage errors

#include <iostream>
#include <fstream>
#include <vector>

int
main(int argc, char* argv[])
{
    if(!argv[1])
    {
        std::cerr << "Please insert filename as argv[1]!" << std::endl;
        std::exit (EXIT_FAILURE);
    }
    
    std::ifstream input;
    std::ofstream output_clean_energy;
    std::ofstream output_errors;
    
    input.open(argv[1]);
    if(!input)
    {
        std::cerr << "Error while opening the input file!" << std::endl;
        std::exit (EXIT_FAILURE);
    }
    
    int N;
    double x, y, z, kinetic_energy, potential_energy;
    long double energy_total, initial_energy_total;
    double energy_error;
    int cycle_counter=0;
    
    std::string str = "energy_clean_";
    str += argv[1];
    output_clean_energy.open(str);
    if(!output_clean_energy)
    {
        std::cerr << "Error while opening the output_clean_energy file!" << std::endl;
        std::exit (EXIT_FAILURE);
    }
    
    std::string str_2 = "energy_errors_";
    str_2 += argv[1];
    output_errors.open(str_2);
    if(!output_errors)
    {
        std::cerr << "Error while opening the output_errors file!" << std::endl;
        std::exit (EXIT_FAILURE);
    }
    
    input >> N;
    
    while (input.good())
    {
        energy_total=0;
        for (int i=0;i<N;i++)
        {
            input >> x >> y >> z >> kinetic_energy >> potential_energy;
            energy_total += kinetic_energy;
            energy_total += potential_energy;
        }
        if(cycle_counter==0)
        {
            initial_energy_total=energy_total;
        }
        
        //error calculated via simple proportion scheme
        energy_error = 100*(initial_energy_total - energy_total) / initial_energy_total;
    
        //handling annoying whitespaces in output
        if(cycle_counter==0)
        {
            output_clean_energy << cycle_counter << '\t' << energy_total;
            output_errors << cycle_counter << '\t' << energy_error;
        }
        else
        {
            output_clean_energy << std::endl << cycle_counter << '\t' << energy_total;
            output_errors << std::endl << cycle_counter << '\t' << energy_error;
        }
        cycle_counter++;
    }
    
    input.close();
    output_clean_energy.close();
    output_errors.close();
    
    return 0;
}
