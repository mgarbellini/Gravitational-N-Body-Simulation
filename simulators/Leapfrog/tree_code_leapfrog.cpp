// F. Borando - francesco.borando@studenti.unimi.it
// M. Garbellini - matteo.garbellini@studenti.unimi.it
// Università degli Studi di Milano
// Dept. of Physics
// January 2019
//
// GRAVITATIONAL N-BODY SIMULATION
// Using tree method with leapfrog algorithm

#include <iostream>
#include <fstream>
#include <vector>
#include "algebra.hpp"
#include "body.hpp"
#include "node.hpp"
#include "tree.hpp"


template <typename DT>
double get_space_width(std::vector< body<DT> > m_all_bodies);

int
main(int argc, char* argv[])
{
    if(!argv[1])
    {
        std::cerr << "Please filename.txt as argv[1]!" << std::endl;
        std::exit (EXIT_FAILURE);
    }

    if(!argv[2])
    {
        std::cerr << "Please timestep as argv[2]!" << std::endl;
        std::exit (EXIT_FAILURE);
    }

    if(!argv[3])
    {
        std::cerr << "Please insert iteration_num as argv[3]!" << std::endl;
        std::exit (EXIT_FAILURE);
    }

    int N=0;
    std::vector< body<double> > all_bodies_vector;
    double t=0;
    const double timestep=std::atof(argv[2]);
    const int iteration_num=std::atof(argv[3]);
    double space_width;
    double mass_to_push, pos_x_coord_to_push, pos_y_coord_to_push, pos_z_coord_to_push, vel_x_coord_to_push, vel_y_coord_to_push, vel_z_coord_to_push;

    //reading of initial values and creation of the bodies
    std::ifstream input;
    std::ofstream output;
    std::ofstream output_velocity;
    input.open(argv[1]);

    if(!input)
    {
        std::cerr << "Error while opening the input file" << std::endl;
        std::exit (EXIT_FAILURE);
    }

    std::string output_name = "computational_res_leapfrog_";
    output_name += argv[1];
    output_name += "_t_";
    output_name += argv[2];
    output_name += ".txt";
    output.open(output_name);

    if(!output)
    {
        std::cerr << "Error while opening the output file" << std::endl;
        std::exit (EXIT_FAILURE);
    }

    std::string output_name_2 = "computational_res_leapfrog_";
    output_name_2 += argv[1];
    output_name_2 += "_t_";
    output_name_2 += argv[2];
    output_name_2 += "_velocity.txt";
    output_velocity.open(output_name_2);


    if(!output_velocity)
    {
        std::cerr << "Error while opening the output_velocity file" << std::endl;
        std::exit (EXIT_FAILURE);
    }

    while (input.good())
    {
        input >> mass_to_push
        >> pos_x_coord_to_push >> pos_y_coord_to_push >> pos_z_coord_to_push
        >> vel_x_coord_to_push >> vel_y_coord_to_push >> vel_z_coord_to_push;

        body<double> body_to_push(mass_to_push, pos_x_coord_to_push, pos_y_coord_to_push, pos_z_coord_to_push, vel_x_coord_to_push, vel_y_coord_to_push, vel_z_coord_to_push);

        all_bodies_vector.push_back(body_to_push);
        N++;
    }

    input.close(); //closing initial-conditions input

    output << N << '\t' << timestep;
    int log_count=0;

    while ( t < iteration_num*timestep )
    {
        space_width = get_space_width(all_bodies_vector); // variable space width (recalculated at every iteration)
        oct_tree<double> m_tree(all_bodies_vector, space_width); // tree reconstructed from scratch at every iteration

        m_tree.evolve(timestep, all_bodies_vector);

        // save positions in output file
        if(log_count!=iteration_num-1)
        {
            for(std::vector< body<double> >::const_iterator it = all_bodies_vector.begin(); it != all_bodies_vector.end(); ++it)
            {
                output << std::endl << it->get_position_x_coord() << '\t' << it->get_position_y_coord() << '\t' << it->get_position_z_coord();
            }
        }

        // last: velocity_backup for checkpoint, save positions AND velocities
        else
        {
            for(std::vector< body<double> >::const_iterator it = all_bodies_vector.begin(); it != (all_bodies_vector.end() - 1); ++it)
            {
                output_velocity << it->get_mass() << '\t' << it->get_position_x_coord() << '\t' << it->get_position_y_coord() << '\t' << it->get_position_z_coord() << '\t' << it->get_velocity_x_coord() << '\t' << it->get_velocity_y_coord() << '\t' << it->get_velocity_z_coord() << std::endl;
            }

            //avoiding undesired std::endl at the and of the file
            std::vector< body<double> >::const_iterator it = all_bodies_vector.end() - 1;
            output_velocity << it->get_mass() << '\t' << it->get_position_x_coord() << '\t' << it->get_position_y_coord() << '\t' << it->get_position_z_coord() << '\t' << it->get_velocity_x_coord() << '\t' << it->get_velocity_y_coord() << '\t' << it->get_velocity_z_coord();
        }

        t+=timestep;

        std::cout << "Iter: " << log_count << std::endl;
        log_count++;
    }

    output.close(); //closing computational-results output

    return 0;
}

template <typename DT>
double get_space_width(std::vector< body<DT> > m_all_bodies)
{
    DT max_x=m_all_bodies[0].get_position_x_coord(),
    max_y=m_all_bodies[0].get_position_y_coord(),
    max_z=m_all_bodies[0].get_position_z_coord(),
    min_x=max_x,
    min_y=max_y,
    min_z=max_z;

    for(typename std::vector< body<DT> >::const_iterator it = m_all_bodies.begin(); it != m_all_bodies.end(); it++)
    {
        if(it->get_position_x_coord() > max_x)
        {
            max_x=it->get_position_x_coord();
        }
        if(it->get_position_y_coord() > max_y)
        {
            max_y=it->get_position_y_coord();
        }
        if(it->get_position_z_coord() > max_z)
        {
            max_z=it->get_position_z_coord();
        }
        if(it->get_position_x_coord() < min_x)
        {
            min_x=it->get_position_x_coord();
        }
        if(it->get_position_y_coord() < min_y)
        {
            min_y=it->get_position_y_coord();
        }
        if(it->get_position_z_coord() < min_z)
        {
            min_z=it->get_position_z_coord();
        }
    }

    //bodies everywhere
    if(max_x>=0 && min_x<=0 && max_y>=0 && min_y<=0 && max_z>=0 && min_z<=0)
    {
        return get_max ( get_max(max_x,-min_x), get_max(max_y,-min_y), get_max(max_z, -min_z))*2.1;
    }
    //all in 4
    if (max_x>=0 && min_x>=0 && max_y>=0 && min_y>=0 && max_z>=0 && min_z>=0)
    {
        return get_max(max_x, max_y, max_z)*2.1;
    }
    //all in 5
    if (max_x<=0 && min_x<=0 && max_y>=0 && min_y>=0 && max_z>=0 && min_z>=0)
    {
        return get_max(-min_x, max_y, max_z)*2.1;
    }
    //all in 1
    if (max_x<=0 && min_x<=0 && max_y<=0 && min_y<=0 && max_z>=0 && min_z>=0)
    {
        return get_max(-min_x, -min_y, max_z)*2.1;
    }
    //all in 0
    if (max_x>=0 && min_x>=0 && max_y<=0 && min_y<=0 && max_z>=0 && min_z>=0)
    {
        return get_max(max_x, -min_y, -min_z)*2.1;
    }
    //all in 7
    if (max_x>=0 && min_x>=0 && max_y>=0 && min_y>=0 && max_z<=0 && min_z<=0)
    {
        return get_max(max_x, max_y, -min_z)*2.1;
    }
    //all in 6
    if (max_x<=0 && min_x<=0 && max_y>=0 && min_y>=0 && max_z<=0 && min_z<=0)
    {
        return get_max(-min_x, max_y, -min_z)*2.1;
    }
    //all in 2
    if (max_x<=0 && min_x<=0 && max_y<=0 && min_y<=0 && max_z<=0 && min_z<=0)
    {
        return get_max(-min_x, -min_y, -min_z)*2.1;
    }
    //all in 3
    if (max_x>=0 && min_x>=0 && max_y<=0 && min_y<=0 && max_z<=0 && min_z<=0)
    {
        return get_max(max_x, -min_y, -min_z)*2.1;
    }
    //all on the RIGHT side - 0 || 4 || 3 || 7
    if (max_x>=0 && min_x>=0)
    {
        return get_max( get_max(max_y, -min_y) , get_max(max_z, -min_z), max_x)*2.1;
    }
    //all on the LEFT side - 1 || 5 || 6 || 2
    if (max_x<=0 && min_x<=0)
    {
        return get_max( get_max(max_y, -min_y), get_max(max_z, -min_z), -min_x)*2.1;
    }
    //all on the UPPER side - 0 || 4 || 1 || 5
    if (max_z>=0 && min_z>=0)
    {
        return get_max( get_max(max_x, -min_x), get_max(max_y, -min_y), max_z)*2.1;
    }
    //all on the LOWER side - 3 || 7 || 2 || 6
    if (max_z<=0 && min_z<=0)
    {
        return get_max( get_max(max_x, -min_x), get_max(max_y, -min_y), -min_z)*2.1;
    }
    //all on the FRONT side - 0 || 4 || 1 || 5
    if (max_y>=0 && min_y>=0)
    {
        return get_max( get_max(max_x, -min_x), get_max(max_z, -min_z), max_y)*2.1;
    }
    //all on the BACK side - 3 || 7 || 2 || 6
    if (max_y<=0 && min_y<=0)
    {
        return get_max( get_max(max_x, -min_x), get_max(max_z, -min_z), -min_y)*2.1;
    }
    else    //previous conditions exhaustive
    {
        std::cerr << "Error - in 'get_space_width()': body is impossible to locate!" << std::endl;
        std::exit (EXIT_FAILURE);
    }
}
