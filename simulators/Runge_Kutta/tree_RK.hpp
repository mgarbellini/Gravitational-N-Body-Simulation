// F. Borando - francesco.borando@studenti.unimi.it
// M. Garbellini - matteo.garbellini@studenti.unimi.it
// Università degli Studi di Milano
// Dept. of Physics
// January 2019
//
// GRAVITATIONAL N-BODY SIMULATION
// tree_RK.hpp: tree class for Barnes-Hut force calculation algorithm
//

#ifndef TREE_FOR_N_BODY_SIMULATION
#define TREE_FOR_N_BODY_SIMULATION

#include <iostream>
#include <memory>
#include <vector>
#include "algebra.hpp"
#include "body_RK.hpp"

//normalized G
#define G_GRAVITY 1 /* N m^2 Kg^-2 */


template <typename DT>
class oct_tree
{
public:
    //constructors
    oct_tree<DT>();
    oct_tree(std::vector< body<DT> >& m_all_bodies, const DT& tree_space_width);
    
    //destructor
    ~oct_tree();
    
    /////////////Tree construction ////////////////
    void update_center_frame_coord (const int& node_case, std::shared_ptr< node<DT> > parent_node);
    void update_node_masses (std::shared_ptr< node<DT> > son_node, const DT mass_to_add);
    void compute_CoM (std::shared_ptr< node<DT> > node_to_compute_CoM);
    void update_CoM (std::shared_ptr< node<DT> > son_node);
    const int get_node_case (const body<DT>& body_to_insert, const V3<DT>& center_frame) const;
    void insert_in_tree (body<DT>& body_to_insert, std::shared_ptr< node<DT> > inspected_node, const int& level, const bool update_CoM_flag);
    
    //////////////Tree evolution /////////////////
    V3<DT> other_body_contribution (const body<DT>& body_subject_to_force, std::shared_ptr< node<DT> > leaf) const;
    V3<DT> other_CoM_contribution (const body<DT>& body_subject_to_force, std::shared_ptr< node<DT> > parent_node) const;
    std::shared_ptr< node<DT> > find_lucky_leaf (std::shared_ptr< node<DT> > s_node);
    const bool barnes_hut_condition (const body<DT>& selected_body_subject, std::shared_ptr< node<DT> > object) const;
    void levels_contribution (const body<DT>& body_subject, std::shared_ptr< node<DT> > l_leaf, const int& m_up_or_under, V3<DT>& accel_support) const;
    V3<DT> compute_acceleration_on_leaf (std::shared_ptr< node<DT> > subject_leaf);
    void clear_all_ticks (std::shared_ptr< node<DT> > m_node);
    void evolve (const DT& timestep, std::vector< body<DT> >& bodies_vector);
    
private:
    std::shared_ptr< node<DT> > m_root;
};

template <typename DT>
oct_tree<DT>::oct_tree()
{
}

template <typename DT>
oct_tree<DT>::oct_tree(std::vector< body<DT> >& m_all_bodies, const DT& tree_space_width)
{
    m_root = std::make_shared< node<DT> >();
    m_root->its_parent=nullptr;
    m_root->level=0;
    m_root->center_frame_coord.set_coord(0,0,0);
    m_root->square_width=tree_space_width;
    m_root->set_children_null();
    for(typename std::vector< body<DT> >::iterator it = m_all_bodies.begin(); it != m_all_bodies.end(); it++)
    {
        insert_in_tree(*it, m_root, 0, true);
    }
}

template <typename DT>
oct_tree<DT>::~oct_tree()
{
}


template <typename DT>
const int oct_tree<DT>::get_node_case(const body<DT>& body_to_insert, const V3<DT>& center_frame) const
{
    //handling upper_z part
    if (body_to_insert.get_position_z_coord() > center_frame.get_z_coord())
    {
        //in upper_z part: handling eastern_y part
        if (body_to_insert.get_position_y_coord() > center_frame.get_y_coord())
        {
            //in upper_z-eastern_y part: handling eastern_x part
            if(body_to_insert.get_position_x_coord() > center_frame.get_x_coord())
            {
                return 4;
            }
            //in upper_z-eastern_y part: handling western_x part
            else if(body_to_insert.get_position_x_coord() < center_frame.get_x_coord())
            {
                return 5;
            }
            else
            {
                return get_random_number_between(4,5);
            }
        }
        //in upper_z part: handling western_y part
        else if(body_to_insert.get_position_y_coord() < center_frame.get_y_coord())
        {
            //in upper_z-western_y part: handling eastern_x part
            if(body_to_insert.get_position_x_coord() > center_frame.get_x_coord())
            {
                return 0;
            }
            //in upper_z-western_y part: handling western_x part
            else if(body_to_insert.get_position_x_coord() < center_frame.get_x_coord())
            {
                return 1;
            }
            else
            {
                return get_random_number_between(0,1);
            }
        }
        //handling upper_z x-axis (note that upper_z y-axis was already handled in 'else' cases up here!)
        else if(body_to_insert.get_position_x_coord() > center_frame.get_x_coord() &&
                body_to_insert.get_position_y_coord() == center_frame.get_y_coord())
        {
            return get_random_number_between(0,4);
        }
        else if(body_to_insert.get_position_x_coord() < center_frame.get_x_coord() &&
                body_to_insert.get_position_y_coord() == center_frame.get_y_coord())
        {
            return get_random_number_between(1,5);
        }
        else if (body_to_insert.get_position_x_coord() == center_frame.get_x_coord() &&
                 body_to_insert.get_position_y_coord() == center_frame.get_y_coord())
        {
            return get_random_number_between(0,1,4,5);
        }
    }
    
    //handling lower_z part
    if (body_to_insert.get_position_z_coord() < center_frame.get_z_coord())
    {
        //in lower_z part: handling eastern_y part
        if (body_to_insert.get_position_y_coord() > center_frame.get_y_coord())
        {
            //in lower_z-eastern_y part: handling eastern_x part
            if(body_to_insert.get_position_x_coord() > center_frame.get_x_coord())
            {
                return 7;
            }
            //in lower_z-eastern_y part: handling western_x part
            else if(body_to_insert.get_position_x_coord() < center_frame.get_x_coord())
            {
                return 6;
            }
            else
            {
                return get_random_number_between(7,6);
            }
        }
        //in lower_z part: handling western_y part
        else if(body_to_insert.get_position_y_coord() < center_frame.get_y_coord())
        {
            //in lower_z-western_y part: handling eastern_x part
            if(body_to_insert.get_position_x_coord() > center_frame.get_x_coord())
            {
                return 3;
            }
            //in lower_z-western_y part: handling western_x part
            else if(body_to_insert.get_position_x_coord() < center_frame.get_x_coord())
            {
                return 2;
            }
            else
            {
                return get_random_number_between(3,2);
            }
        }
        //handling lower_z x-axis
        else if(body_to_insert.get_position_x_coord() > center_frame.get_x_coord() &&
                body_to_insert.get_position_y_coord() == center_frame.get_y_coord())
        {
            return get_random_number_between(3,7);
        }
        else if(body_to_insert.get_position_x_coord() < center_frame.get_x_coord() &&
                body_to_insert.get_position_y_coord() == center_frame.get_y_coord())
        {
            return get_random_number_between(6,2);
        }
        else if (body_to_insert.get_position_x_coord() == center_frame.get_x_coord() &&
                 body_to_insert.get_position_y_coord() == center_frame.get_y_coord())
        {
            return get_random_number_between(2,3,6,7);
        }
    }
    
    //hanling plane z=0
    else if(body_to_insert.get_position_z_coord() == center_frame.get_z_coord())
    {
        if(body_to_insert.get_position_y_coord() > center_frame.get_y_coord())
        {
            if (body_to_insert.get_position_x_coord() > center_frame.get_x_coord())
            {
                return get_random_number_between(4,7);
            }
            else if(body_to_insert.get_position_x_coord() < center_frame.get_x_coord())
            {
                return get_random_number_between(5,6);
            }
            else
            {
                return get_random_number_between(4,5,6,7);
            }
        }
        else if(body_to_insert.get_position_y_coord() < center_frame.get_y_coord())
        {
            if (body_to_insert.get_position_x_coord() > center_frame.get_x_coord())
            {
                return get_random_number_between(0,3);
            }
            else if(body_to_insert.get_position_x_coord() < center_frame.get_x_coord())
            {
                return get_random_number_between(1,2);
            }
            else
            {
                return get_random_number_between(0,1,2,3);
            }
        }
        else if (body_to_insert.get_position_y_coord() == center_frame.get_y_coord())
        {
            if (body_to_insert.get_position_x_coord() > center_frame.get_x_coord())
            {
                return get_random_number_between(4,7,0,3);
            }
            else if(body_to_insert.get_position_x_coord() < center_frame.get_x_coord())
            {
                return get_random_number_between(5,6,1,2);
            }
            //handling the origin (center of axes)
            else
            {
                return get_random_number_between(0,1,2,3,4,5,6,7);
            }
        }
    }
    
    //previous conditions are exhaustive
    std::cerr << "Error - in 'oct_tree<DT>::get_node_case()': node case identification failed!" << std::endl;
    std::exit (EXIT_FAILURE);
    return EXIT_FAILURE;
}

template <typename DT>
void oct_tree<DT>::update_center_frame_coord(const int& node_case, std::shared_ptr< node<DT> > parent_node)
{
    parent_node->children[node_case]->center_frame_coord.set_coord(parent_node->center_frame_coord.get_x_coord(),
                                                                   parent_node->center_frame_coord.get_y_coord(),
                                                                   parent_node->center_frame_coord.get_z_coord());
    switch (node_case){
        case 0:
            parent_node->children[node_case]->center_frame_coord.add_x_coord(parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_y_coord(- parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_z_coord(parent_node->children[node_case]->square_width/2);
            break;
        case 1:
            parent_node->children[node_case]->center_frame_coord.add_x_coord(- parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_y_coord(- parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_z_coord(parent_node->children[node_case]->square_width/2);
            break;
        case 2:
            parent_node->children[node_case]->center_frame_coord.add_x_coord(- parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_y_coord(- parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_z_coord(- parent_node->children[node_case]->square_width/2);
            break;
        case 3:
            parent_node->children[node_case]->center_frame_coord.add_x_coord(parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_y_coord(- parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_z_coord(- parent_node->children[node_case]->square_width/2);
            break;
        case 4:
            parent_node->children[node_case]->center_frame_coord.add_x_coord(parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_y_coord(parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_z_coord(parent_node->children[node_case]->square_width/2);
            break;
        case 5:
            parent_node->children[node_case]->center_frame_coord.add_x_coord(- parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_y_coord(parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_z_coord(parent_node->children[node_case]->square_width/2);
            break;
        case 6:
            parent_node->children[node_case]->center_frame_coord.add_x_coord(- parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_y_coord(parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_z_coord(- parent_node->children[node_case]->square_width/2);
            break;
        case 7:
            parent_node->children[node_case]->center_frame_coord.add_x_coord(parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_y_coord(parent_node->children[node_case]->square_width/2);
            parent_node->children[node_case]->center_frame_coord.add_z_coord(- parent_node->children[node_case]->square_width/2);
            break;
        default:
            std::cerr<<"Error - in 'oct_tree<DT>::update_center_frame_coord()': node case identification probably failed!"<<std::endl;
            std::exit (EXIT_FAILURE);
    }
}

template <typename DT>
void oct_tree<DT>::update_node_masses(std::shared_ptr< node<DT> > son_node, const DT mass_to_add)
{
    if (son_node==m_root)
    {
        return;
    }
    
    else
    {
        son_node->its_parent->node_mass += mass_to_add;
        update_node_masses(son_node->its_parent, mass_to_add);
    }
}

template <typename DT>
void oct_tree<DT>::compute_CoM(std::shared_ptr< node<DT> > node_to_compute_CoM)
{
    V3<DT> numerator(0,0,0);
    DT denominator=0;;
    
    for(int i=0;i<8;i++)
    {
        //empty node
        if(node_to_compute_CoM->children[i]==nullptr)
        {
            continue;
        }
        //leaf node or parent
        else if(node_to_compute_CoM->children[i]->is_leaf==true || node_to_compute_CoM->children[i]->is_parent==true)
        {
            V3<DT> num_add = mult_by_scalar(node_to_compute_CoM->children[i]->CoM_position,node_to_compute_CoM->children[i]->node_mass);
            numerator += num_add;
            denominator += node_to_compute_CoM->children[i]->node_mass;
        }
        else    //previous conditions are exhaustive
        {
            std::cerr << "Error - in 'oct_tree<DT>::computeCoM()': a node is nor a leaf nor a parent nor nullptr!!" <<std::endl;
            std::exit (EXIT_FAILURE);
        }
    }
    
    V3<DT> fraction = mult_by_scalar(numerator, 1/denominator);
    node_to_compute_CoM->CoM_position=fraction;
}

template <typename DT>
void oct_tree<DT>::update_CoM(std::shared_ptr< node<DT> > son_node)
{
    if (son_node->its_parent==m_root)
    {
        return;
    }
    else
    {
        compute_CoM(son_node->its_parent);
    }
}

//leaf passed is always the father
template <typename DT>
void oct_tree<DT>::insert_in_tree(body<DT>& body_to_insert, std::shared_ptr< node<DT> > inspected_node, const int& level, const bool update_CoM_flag)
{
    int node_case = get_node_case(body_to_insert, inspected_node->center_frame_coord);
    
    //inspected_node->children==nullptr
    if (!inspected_node->children[node_case])
    {
        //heap allocation, for very big trees
        inspected_node->children[node_case]= std::make_shared< node<DT> >();
        if(!inspected_node->children[node_case])
        {
            std::cerr << "Error - in 'oct_tree<DT>::insert_in_tree()': Memory allocation failed" << std::endl;
            std::exit (EXIT_FAILURE);
        }
        
        inspected_node->children[node_case]->m_body=&body_to_insert;
        
        //the old leaf-node turns to a parent permanently (until the destruction of the tree)
        inspected_node->is_parent=true;
        inspected_node->is_leaf=false;
        inspected_node->m_body=NULL;
        inspected_node->children[node_case]->set_children_null();
        
        //the new node becomes a leaf (and remains a leaf until it eventually turns to a parent)
        inspected_node->children[node_case]->its_parent = inspected_node;
        inspected_node->children[node_case]->is_leaf = true;
        inspected_node->children[node_case]->is_parent = false;
        inspected_node->children[node_case]->level=level+1;
        inspected_node->children[node_case]->tick=false;
        inspected_node->children[node_case]->node_mass=body_to_insert.get_mass();
        //CoM of leafs coincides with the body position that they point
        inspected_node->children[node_case]->CoM_position=body_to_insert.get_position();
        inspected_node->children[node_case]->square_width = m_root->square_width/std::pow(2, inspected_node->children[node_case]->level);
    }
    
    else if(inspected_node->children[node_case]->is_parent==true)
    {
        insert_in_tree(body_to_insert, inspected_node->children[node_case], inspected_node->children[node_case]->level, true);
    }
    
    else if (inspected_node->children[node_case]->is_leaf==true)
    {
        //updates center frame of inspected_node->children[node_case]
        update_center_frame_coord(node_case, inspected_node);
        
        insert_in_tree(*(inspected_node->children[node_case]->m_body), inspected_node->children[node_case], inspected_node->children[node_case]->level, false);
        
        insert_in_tree(body_to_insert, inspected_node->children[node_case], inspected_node->children[node_case]->level, true);
    }
    if(update_CoM_flag==true)
    {
        //adds to parents node-mass, so it starts with the newborn leaf
        update_node_masses(inspected_node->children[node_case], body_to_insert.get_mass());
        update_CoM(inspected_node->children[node_case]);
    }
}

////////////////////////////////////////////////
//TREE CREATION

//**********************************************

//TREE EVOLUTION
////////////////////////////////////////////////


template <typename DT>
std::shared_ptr< node<DT> > oct_tree<DT>::find_lucky_leaf(std::shared_ptr< node<DT> > s_node)
{
    for(int i=0;i<8;i++)
    {
        if (s_node->children[i]==nullptr)
        {
            continue;
        }
        else if (s_node->children[i]->is_leaf==true && s_node->children[i]->tick==false)
        {
            s_node->children[i]->tick=true;
            return s_node->children[i];
        }
        else if (s_node->children[i]->is_leaf==true && s_node->children[i]->tick==true)
        {
            continue;
        }
        else if (s_node->children[i]->is_parent==true && s_node->children[i]->tick==false)
        {
            //need to make sure that the all the children of the parent are ticked or not
            //in fact 'ticking' of parents must happen somewhere: HERE!
            if ( s_node->children[i]->children_all_ticked() )
            {
                s_node->children[i]->tick=true;
                continue;
            }
            else
            {
                return find_lucky_leaf(s_node->children[i]);
            }
        }
        else if (s_node->children[i]->is_parent==true && s_node->children[i]->tick==true)
        {
            continue;
        }
        else //previous conditions are exhaustive, last condition to handle errors
        {
            std::cerr << "Error - in 'oct_tree<DT>::find_lucky_leaf()': a node is nor a leaf nor a parent nor nullptr!" << std::endl;
            std::exit (EXIT_FAILURE);
        }
    }
    if(s_node==m_root)
    {
        //gets here iff all the nodes of the tree are ticked
        return m_root;
    }
    else
    {
        return nullptr;
    }
}

template <typename DT>
V3<DT> oct_tree<DT>::other_body_contribution (const body<DT>& body_subject_to_force, std::shared_ptr< node<DT> > leaf) const
{
    const DT d = distance_between(body_subject_to_force.get_position(), leaf->m_body->get_position());
    const DT r_soft = std::pow(d*d + 0.025,0.5);
    const DT accel_tmp =  ( leaf->m_body->get_mass() / std::pow(r_soft,3) ) * G_GRAVITY;
    const V3<DT> acc_contribution(accel_tmp * ( leaf->m_body->get_position_x_coord() - body_subject_to_force.get_position_x_coord() ),
                                  accel_tmp * ( leaf->m_body->get_position_y_coord() - body_subject_to_force.get_position_y_coord() ),
                                  accel_tmp * ( leaf->m_body->get_position_z_coord() - body_subject_to_force.get_position_z_coord() ));
    return acc_contribution;
}

template <typename DT>
V3<DT> oct_tree<DT>::other_CoM_contribution (const body<DT>& body_subject_to_force, std::shared_ptr< node<DT> > parent_node) const
{
    const DT d = distance_between(body_subject_to_force.get_position(),parent_node->CoM_position);
    const DT r_soft = std::pow(d*d + 0.25, 0.5);
    const DT accel_tmp =  ( parent_node->node_mass / std::pow(r_soft,3) ) * G_GRAVITY;
    const V3<DT> acc_contribution(accel_tmp * ( parent_node->CoM_position.get_x_coord() - body_subject_to_force.get_position_x_coord() ),
                                  accel_tmp * ( parent_node->CoM_position.get_y_coord() - body_subject_to_force.get_position_y_coord() ),
                                  accel_tmp * ( parent_node->CoM_position.get_z_coord() - body_subject_to_force.get_position_z_coord() ));
    return acc_contribution;
}

template <typename DT>
const bool oct_tree<DT>::barnes_hut_condition (const body<DT>& selected_body_subject, std::shared_ptr< node<DT> > object) const
{
    if ( ( object->square_width / distance_between(selected_body_subject.get_position(), object->CoM_position) ) < 0.5 )
    {
        return true;
    }
    else
    {
        return false;
    }
}

//levels_contribution is a void function that modifies V3<DT> accel_support, which is in fact always passed by reference to mantain value
//std::shared_ptr< node<DT> > l_leaf is always the "children node" in the level inspected - initially l_leaf=lucky_leaf

enum direction {down,up};

template <typename DT>
void oct_tree<DT>::levels_contribution(const body<DT>& body_subject, std::shared_ptr< node<DT> > l_leaf, const int& m_up_or_under, V3<DT>& accel_support) const
{
    if(l_leaf==m_root)
    {
        return; //returns iff all levels' contributions were taken into account
    }
    
    std::shared_ptr< node<DT> > parent_of_l_leaf = l_leaf->its_parent;
    
    switch (m_up_or_under)
    {
        case 0: //going down
        {
            for(int i=0;i<8;i++)
            {
                if(!parent_of_l_leaf->children[i])
                {
                    continue;
                }
                else if(parent_of_l_leaf->children[i]->is_leaf==true)
                {
                    V3<DT> accel_add = other_body_contribution(body_subject, parent_of_l_leaf->children[i]);
                    accel_support += accel_add;
                }
                else if(parent_of_l_leaf->children[i]->is_parent==true)
                {
                    if ( barnes_hut_condition(body_subject, parent_of_l_leaf->children[i]) )
                    {
                        V3<DT> accel_add = other_CoM_contribution(body_subject, parent_of_l_leaf->children[i]);
                        accel_support += accel_add;
                    }
                    
                    else
                    {
                        //random child of parent_of_l_leaf->children[i] is passed. Must not be nullptr to avoid segmentation fault
                        levels_contribution(body_subject, parent_of_l_leaf->children[i]->children[parent_of_l_leaf->children[i]->find_first_non_NULL_child()], down, accel_support);
                    }
                }
                else
                {
                    std::cerr << "Error - in 'oct_tree<DT>::levels_contribution()': error during node research!" <<std::endl;
                    std::exit (EXIT_FAILURE);
                }
            }
        }
            break; //switch: end case down
            
        case 1: //going up
        {
            int not_to_compute=l_leaf->get_children_id(l_leaf);
            
            for(int i=0;i<8;i++)
            {
                if(!parent_of_l_leaf->children[i])
                {
                    continue;
                }
                
                else if (i==not_to_compute)
                {
                    continue;
                }
                
                else if(parent_of_l_leaf->children[i]->is_leaf==true)
                {
                    V3<DT> accel_add = other_body_contribution(body_subject, parent_of_l_leaf->children[i]);
                    accel_support += accel_add;
                    
                }
                else if(parent_of_l_leaf->children[i]->is_parent==true)
                {
                    if (barnes_hut_condition(body_subject, parent_of_l_leaf->children[i]))
                    {
                        V3<DT> accel_add = other_CoM_contribution(body_subject, parent_of_l_leaf->children[i]);
                        accel_support += accel_add;
                    }
                    
                    else
                    {
                        levels_contribution(body_subject, parent_of_l_leaf->children[i]->children[parent_of_l_leaf->children[i]->find_first_non_NULL_child()], down, accel_support);
                    }
                }
                else
                {
                    std::cerr << "Error - in 'oct_tree<DT>::levels_contribution()': error during node research!" <<std::endl;
                    std::exit (EXIT_FAILURE);
                }
            }
            //'horizontal' family is complete: every contribution of the family under parent_of_l_leaf was considered
            //must go up until the parent of l_leaf is the parent of m_root
            //When stepping up in tree levels "l_leaf" analogically represents the subject_leaf. Operatively it's a node that must not be taken into account when evaluating forces from all the members of a family under a specific parent_of_l_leaf node.
            //When called recursively, parent_of_l_leaf takes place of initial l_leaf node. So, in the recursive call, parent_of_l_leaf would actually be a 'parent_of_parent_of_l_leaf', and so on - until the function calls m_root.
            
            levels_contribution(body_subject, parent_of_l_leaf, up, accel_support);
        }
            break; //switch: end case up
            
        default:
            break;
    }// end of switch
}

template <typename DT>
V3<DT> oct_tree<DT>::compute_acceleration_on_leaf(std::shared_ptr< node<DT> > subject_leaf)
{
    if(subject_leaf==m_root)
    {
        std::cerr << "Error in 'oct_tree<DT>::compute_acceleration_on_leaf()' - subject_leaf==m_root!" << std::endl;
        std::exit (EXIT_FAILURE);
    }
    
    V3<DT> total_acceleration(0,0,0);
    levels_contribution(*(subject_leaf->m_body), subject_leaf, up, total_acceleration);
    return total_acceleration;
}

template <typename DT>
void oct_tree<DT>::clear_all_ticks(std::shared_ptr< node<DT> > m_node)
{
    for (int i=0;i<8;i++)
    {
        if(m_node->children[i]==nullptr)
        {
            continue;
        }
        else if(m_node->children[i]->is_leaf==true)
        {
            m_node->children[i]->tick=false;
        }
        else if (m_node->children[i]->is_parent==true)
        {
            clear_all_ticks(m_node->children[i]);
        }
        else    //previous conditions are exhaustive
        {
            std::cerr<<"Error - in 'oct_tree<DT>::clear_all_ticks()': a node is nor a leaf nor a parent nor nullptr!" <<std::endl;
            std::exit (EXIT_FAILURE);
        }
    }
    m_node->tick=false;
    
}

//Runge-Kutta 4 scheme
template <typename DT>
void oct_tree<DT>::evolve(const DT& t_epsilon, std::vector< body<DT> >& bodies_vector)
{
    std::shared_ptr< node<DT> > lucky_leaf = find_lucky_leaf(m_root);
    //body_shape is a support body-object: initially replicates lucky_leaf->m_body
    body<DT> body_shape;
    
    //Useful for code readability
    double x_0 ,y_0, z_0, v_0x, v_0y, v_0z,
    k1_x, k1_y, k1_z, k1_vx, k1_vy, k1_vz,
    k2_x, k2_y, k2_z, k2_vx, k2_vy, k2_vz,
    k3_x, k3_y, k3_z, k3_vx, k3_vy, k3_vz,
    k4_x, k4_y, k4_z, k4_vx, k4_vy, k4_vz;
    
    while(lucky_leaf!=m_root)
    {
        body_shape.set_mass(lucky_leaf->m_body->get_mass());
        body_shape.set_position(lucky_leaf->m_body->get_position_x_coord(),
                                lucky_leaf->m_body->get_position_y_coord(),
                                lucky_leaf->m_body->get_position_z_coord());
        body_shape.set_velocity(lucky_leaf->m_body->get_velocity_x_coord(),
                                lucky_leaf->m_body->get_velocity_y_coord(),
                                lucky_leaf->m_body->get_velocity_z_coord());
        
        x_0 = body_shape.get_position_x_coord();
        y_0 = body_shape.get_position_y_coord();
        z_0 = body_shape.get_position_z_coord();
        v_0x = body_shape.get_velocity_x_coord();
        v_0y = body_shape.get_velocity_y_coord();
        v_0z = body_shape.get_velocity_z_coord();
        
        //K1
        
        k1_x = v_0x;
        k1_y = v_0y;
        k1_z = v_0z;
        
        body_shape.set_acceleration_by_vector(compute_acceleration_on_leaf(lucky_leaf));
        
        k1_vx = body_shape.get_acceleration_x_coord();
        k1_vy = body_shape.get_acceleration_y_coord();
        k1_vz = body_shape.get_acceleration_z_coord();
        
        // K2
        
        k2_x = v_0x + 0.5*t_epsilon*k1_vx;
        k2_y = v_0y + 0.5*t_epsilon*k1_vy;
        k2_z = v_0z + 0.5*t_epsilon*k1_vz;
        
        body_shape.set_position(x_0 + (t_epsilon*k1_x)/2, y_0 + (t_epsilon*k1_y)/2, z_0 + (t_epsilon*k1_z)/2);
        
        body_shape.set_acceleration_by_vector(compute_acceleration_on_leaf(lucky_leaf));
        
        k2_vx = body_shape.get_acceleration_x_coord();
        k2_vy = body_shape.get_acceleration_y_coord();
        k2_vz = body_shape.get_acceleration_z_coord();
        
        
        // K3
        
        k3_x = v_0x + 0.5*t_epsilon*k2_vx;
        k3_y = v_0y + 0.5*t_epsilon*k2_vy;
        k3_z = v_0z + 0.5*t_epsilon*k2_vz;
        
        body_shape.set_position(x_0 + (t_epsilon*k2_x)/2, y_0 + (t_epsilon*k2_y)/2, z_0 + (t_epsilon*k2_z)/2);
        
        body_shape.set_acceleration_by_vector(compute_acceleration_on_leaf(lucky_leaf));
        
        k3_vx = body_shape.get_acceleration_x_coord();
        k3_vy = body_shape.get_acceleration_y_coord();
        k3_vz = body_shape.get_acceleration_z_coord();
        
        // K4
        
        k4_x = v_0x + t_epsilon*k3_vx;
        k4_y = v_0y + t_epsilon*k3_vy;
        k4_z = v_0z + t_epsilon*k3_vz;
        
        body_shape.set_position(x_0 + t_epsilon*k3_x, y_0 + t_epsilon * k3_y, z_0 + t_epsilon * k3_z);
        
        body_shape.set_acceleration_by_vector(compute_acceleration_on_leaf(lucky_leaf));
        
        k4_vx = body_shape.get_acceleration_x_coord();
        k4_vy = body_shape.get_acceleration_y_coord();
        k4_vz = body_shape.get_acceleration_z_coord();
        
        lucky_leaf->m_body->set_next_pos(
                                         lucky_leaf->m_body->get_position_x_coord() + t_epsilon*(k1_x + 2*k2_x + 2* k3_x + k4_x)/6,
                                         lucky_leaf->m_body->get_position_y_coord() + t_epsilon*(k1_y + 2*k2_y + 2* k3_y + k4_y)/6,
                                         lucky_leaf->m_body->get_position_z_coord() + t_epsilon*(k1_z + 2*k2_z + 2* k3_z + k4_z)/6
                                         );
        lucky_leaf->m_body->set_next_vel(
                                         lucky_leaf->m_body->get_velocity_x_coord() + t_epsilon*(k1_vx + 2*k2_vx + 2* k3_vx + k4_vx)/6,
                                         lucky_leaf->m_body->get_velocity_y_coord() + t_epsilon*(k1_vy + 2*k2_vy + 2* k3_vy + k4_vy)/6,
                                         lucky_leaf->m_body->get_velocity_z_coord() + t_epsilon*(k1_vz + 2*k2_vz + 2* k3_vz + k4_vz)/6);
        
        lucky_leaf=nullptr;
        while(lucky_leaf==nullptr)
        {
            lucky_leaf=find_lucky_leaf(m_root);
        }
    }
    clear_all_ticks(m_root);
    
    //updates all bodies' pos and vel
    for(typename std::vector< body<DT> >::iterator it = bodies_vector.begin(); it != bodies_vector.end(); it++)
    {
        it->set_velocity_by_vector( it->get_next_vel() );
        it->set_position_by_vector( it->get_next_pos() );
    }
}

////////////////////////////////////////////////
//TREE EVOLUTION

#endif //TREE_FOR_N_BODY_SIMULATION
