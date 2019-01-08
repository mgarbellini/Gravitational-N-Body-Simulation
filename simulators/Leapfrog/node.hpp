// F. Borando - francesco.borando@studenti.unimi.it
// M. Garbellini - matteo.garbellini@studenti.unimi.it
// Università degli Studi di Milano
// Dept. of Physics
// January 2019
//
// GRAVITATIONAL N-BODY SIMULATION
// node.hpp: include file for using node class
// this class is needed for tree.hpp

#ifndef NODE_FOR_N_BODY_SIMULATION
#define NODE_FOR_N_BODY_SIMULATION

#include <iostream>
#include <memory>
#include <vector>
#include "algebra.hpp"
#include "body.hpp"

template <typename DT>
class node
{
public:
    node()
    {
    }
    
    ~node()
    {
    }
    
    inline void set_children_null();
    const bool children_all_ticked() const;
    inline const int get_children_id(const std::shared_ptr< node<DT> > this_node) const;
    inline const int find_first_non_NULL_child() const;
    
    std::shared_ptr< node<DT> > its_parent;
    std::vector < std::shared_ptr< node<DT> > > children;
    bool is_leaf;
    bool is_parent;
    bool tick;
    int level;
    body<DT>* m_body;
    DT node_mass;
    DT square_width;
    V3<DT> CoM_position;
    V3<DT> center_frame_coord;
};

template <typename DT>
inline void node<DT>::set_children_null()
{
    std::shared_ptr< node<DT> > null_ptr = std::make_shared< node<DT> >();
    null_ptr=nullptr;
    for(int i=0;i<8;i++)
    {
        this->children.push_back(null_ptr);
    }
}

template <typename DT>
const bool node<DT>::children_all_ticked() const
{
    for(int i=0;i<8;i++)
    {
        if (this->children[i] == nullptr)
        {
            continue;
        }
        
        else if (this->children[i]->tick==false)
        {
            return false;
        }
        
        else if (this->children[i]->tick==true)
        {
            continue;
        }
    }
    return true;
}

template <typename DT>
inline const int node<DT>::get_children_id(const std::shared_ptr< node<DT> > this_node) const
{
    for (int i=0;i<8;i++)
    {
        if( this->its_parent->children[i] == this_node )
        {
            return i;
        }
    }
    //if it gets here, didn't find any correspondance
    std::cerr << "Error - in 'oct_tree<DT>::get_children_id': couldn't find children id!" << std::endl;
    std::exit (EXIT_FAILURE);
}

template <typename DT>
inline const int node<DT>::find_first_non_NULL_child() const
{
    for(int i=0;i<8;i++)
    {
        if(this->children[i] != nullptr)
        {
            return i;
        }
    }
    //gets here if did not find any children.
    std::cerr<<"Error - in 'oct_tree<DT>::find_first_non_NULL_child': did not find any children!" << std::endl;
    std::exit (EXIT_FAILURE);
}

#endif //NODE_FOR_N_BODY_SIMULATION
