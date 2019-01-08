// F. Borando - francesco.borando@studenti.unimi.it
// M. Garbellini - matteo.garbellini@studenti.unimi.it
// Università degli Studi di Milano
// Dept. of Physics
// January 2019
//
// GRAVITATIONAL N-BODY SIMULATION
// algebra.hpp: include file for using V3 vectors, related functions and simple analysis functions.
// this class is needed for body.hpp and tree.hpp

#ifndef ALGEBRAIC_TOOLS_FOR_N_BODY_SIMULATION
#define ALGEBRAIC_TOOLS_FOR_N_BODY_SIMULATION

#include <iostream>
#include <cmath>
#include <cstdlib> // std::rand()
#include <ctime> // std::srand()

template <typename DT>
class V3
 {
 public:
 	V3()
    {
    }

 	V3(const DT& x, const DT& y, const DT& z) :
 		m_x(x), m_y(y), m_z(z)
    {
 	 	m_norm = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
 	 	m_theta = std::acos(z/sqrt(pow(y,2) + pow(z,2)));
 	 	m_phi = std::acos(x/sqrt(pow(x,2) + pow(y,2)));
    }

 	~V3()
    {
    }

 	void set_coord(const DT& x, const DT& y, const DT& z)
    {
 	 	m_x=x;
 	 	m_y=y;
 	 	m_z=z;
 	 	m_norm = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
 	 	m_theta = std::acos(z/m_norm);
 	 	m_phi = std::acos (x/sqrt(pow(x,2) + pow(y,2)));
    }
     
 	inline const DT get_x_coord() const {return m_x;}
 	inline const DT get_y_coord() const {return m_y;}
 	inline const DT get_z_coord() const {return m_z;}
 	inline const DT get_norm() const {return m_norm;}
 	inline const DT get_theta() const {return m_theta;}
 	inline const DT get_phi() const {return m_phi;}
    void set_x_coord(const DT& x)
    {
        m_x=x;
        m_norm = sqrt(pow(m_x,2) + pow(m_y,2) + + pow(m_z,2));
    }
    void set_y_coord(const DT& y)
    {
        m_y=y ;
        m_norm = sqrt(pow(m_x,2) + pow(m_y,2) + pow(m_z,2));
    }
    void set_z_coord(const DT& z)
    {
        m_z=z;
        m_norm = sqrt(pow(m_x,2) + pow(m_y,2) + pow(m_z,2));
    }
    void add_x_coord(const DT& x)
    {
        m_x = m_x + x;
    }
    void add_y_coord(const DT& y)
    {
        m_y = m_y + y;
    }
    void add_z_coord(const DT& z)
    {
        m_z = m_z + z;
    }
     
    V3<DT> operator+ (const V3<DT>& v_1);
    V3<DT> operator- (const V3<DT>& v_1);
    V3<DT> operator+= (const V3<DT>& v_1);
    V3<DT> operator-= (const V3<DT>& v_1);
    DT operator* (const V3<DT>& v_1); //dot product
     
 private:
 	DT m_x;
 	DT m_y;
 	DT m_z;
 	DT m_norm;
 	DT m_theta; //Height angle
 	DT m_phi; //Azimut angle
 };

template <typename DT>
V3<DT> V3<DT>::operator+ (const V3<DT>& v_1)
{
   V3<DT> support_vector(get_x_coord(), get_y_coord(), get_z_coord());
        support_vector.set_coord(support_vector.get_x_coord() + v_1.get_x_coord(),
                                 support_vector.get_y_coord() + v_1.get_y_coord(),
                                 support_vector.get_z_coord() + v_1.get_z_coord());
    return support_vector;
}

template <typename DT>
V3<DT> V3<DT>::operator+= (const V3<DT>& v_1)
{
    this->m_x=m_x+v_1.get_x_coord();
    this->m_y=m_y+v_1.get_y_coord();
    this->m_z=m_z+v_1.get_z_coord();
    return *this;
}

template <typename DT>
V3<DT> V3<DT>::operator-= (const V3<DT>& v_1)
{
    this->m_x=m_x-v_1.get_x_coord();
    this->m_y=m_y-v_1.get_y_coord();
    this->m_z=m_z-v_1.get_z_coord();
    return *this;
}

template <typename DT>
V3<DT> V3<DT>::operator- (const V3<DT>& v_1) //note that order counts
{
   V3<DT> support_vector(get_x_coord(), get_y_coord(), get_z_coord());
        support_vector.set_coord(support_vector.get_x_coord() - v_1.get_x_coord(),
                                 support_vector.get_y_coord() - v_1.get_y_coord(),
                                 support_vector.get_z_coord() - v_1.get_z_coord());
    return support_vector;
}

//dot product
template <typename DT>
DT V3<DT>::operator* (const V3<DT>& v_1)
{
   
    V3<DT> v_2 (get_x_coord(), get_y_coord(), get_z_coord());
        DT dot_product = v_2.get_x_coord() * v_1.get_x_coord() +
        v_2.get_y_coord() * v_1.get_y_coord() +
        v_2.get_z_coord() * v_1.get_z_coord();
    return dot_product;
}

//Useful function
template <typename DT>
const V3<DT> mult_by_scalar (const V3<DT>& point, const DT& scalar)
{
    V3<DT> vector_after_operation;
    vector_after_operation.set_coord(point.get_x_coord() * scalar, point.get_y_coord() * scalar, point.get_z_coord() * scalar);
    return vector_after_operation;
}

//copy to avoid problems with order of the arguments - multiplication by scalar is symmetric
template <typename DT>
const V3<DT> mult_by_scalar (const DT& scalar, const V3<DT>& point)
{
    V3<DT> vector_after_operation;
    vector_after_operation.set_coord(point.get_x_coord() * scalar, point.get_y_coord() * scalar, point.get_z_coord() * scalar);
    return vector_after_operation;
}

template <typename DT>
const DT distance_between(const V3<DT>& point_1, const V3<DT>& point_2)
{
    return std::sqrt( std::pow(point_1.get_x_coord() - point_2.get_x_coord(),2)  + std::pow(point_1.get_y_coord() - point_2.get_y_coord(),2) + std::pow(point_1.get_z_coord() - point_2.get_z_coord(),2) );

}

int get_random_number_between(int n_1, int n_2)
{
    std::srand(time(NULL));
    int random_number=std::rand() % 2;
    if (random_number==0)
    {
        return n_1;
    }
    else if (random_number==1)
    {
        return n_2;
    }
    else
    {
        std::cerr << "Error in get_random_number 2-cases!" << std::endl;
        std::exit (EXIT_FAILURE);
    }
}

int get_random_number_between(int n_1, int n_2, int n_3, int n_4)
{
    std::srand(time(NULL));
    int random_number=std::rand() % 4;
    switch (random_number)
    {
        case 0:
            return n_1;
            break;
        case 1:
            return n_2;
            break;
        case 2:
            return n_3;
            break;
        case 3:
            return n_4;
            break;
        default:
            std::cerr << "Error in get_random_number 4-cases!" << std::endl;
            std::exit (EXIT_FAILURE);
            break;
    }
}

int get_random_number_between(int n_1, int n_2, int n_3, int n_4,
                              int n_5, int n_6, int n_7, int n_8)
{
    std::srand(time(NULL));
    int random_number=std::rand() % 8;
    switch (random_number)
    {
        case 0:
            return n_1;
            break;
        case 1:
            return n_2;
            break;
        case 2:
            return n_3;
            break;
        case 3:
            return n_4;
            break;
        case 4:
            return n_5;
            break;
        case 5:
            return n_6;
            break;
        case 6:
            return n_7;
            break;
        case 7:
            return n_8;
            break;
        default:
            std::cerr << "Error in get_random_number 8-cases!" << std::endl;
            std::exit (EXIT_FAILURE);
            break;
    }
}

double get_max(double a, double b)
{
    return ( (a<b) ? b : a);
}

double get_max(double a, double b, double c)
{
    int max = (a<b) ? b : a;
    return ( (max<c) ? c : max );
}

 #endif //ALGEBRAIC_TOOLS_FOR_N_BODY_SIMULATION
