// F. Borando - francesco.borando@studenti.unimi.it
// M. Garbellini - matteo.garbellini@studenti.unimi.it
// Università degli Studi di Milano
// Dept. of Physics
// January 2019
//
// GRAVITATIONAL N-BODY SIMULATION
// body.hpp: include file for using body class
// the class represents the abstract point-particle of classical mechanics
// this class is needed for tree_RK.hpp and in all main functions (in every *.cpp file) using RK algorithm


#ifndef BODY_CLASS_FOR_N_BODY_SIMULATION
#define BODY_CLASS_FOR_N_BODY_SIMULATION

#include "algebra.hpp"

template <typename DT>
class body
{
public:
    
    body() : m_mass(0)
    {
        m_position.set_coord(0,0,0);
        m_velocity.set_coord(0,0,0);
        m_acceleration.set_coord(0,0,0);
    }
    
    body(const DT& mass, const V3<DT>& initial_position, const V3<DT>& initial_velocity)
    : m_mass(mass)
    {
        m_position = initial_position;
        m_velocity = initial_velocity;
        m_acceleration.set_coord(0,0,0);
    }
    
    body(const DT& mass, const DT& initial_pos_x, const DT& initial_pos_y, const DT& initial_pos_z, const DT& initial_vel_x, const DT& initial_vel_y, const DT& initial_vel_z)
    : m_mass(mass)
    {
        m_position.set_coord(initial_pos_x, initial_pos_y, initial_pos_z);
        m_velocity.set_coord(initial_vel_x, initial_vel_y, initial_vel_z);
        m_acceleration.set_coord(0,0,0);
    }
    
    ~body()
    {
    }
    
    inline const DT get_mass() const {return m_mass;}
    
    inline const V3<DT> get_position() const {return m_position;}
    inline const DT get_position_x_coord() const {return m_position.get_x_coord();}
    inline const DT get_position_y_coord() const {return m_position.get_y_coord();}
    inline const DT get_position_z_coord() const {return m_position.get_z_coord();}
    
    inline const V3<DT> get_velocity() const {return m_velocity;}
    inline const DT get_velocity_x_coord() const {return m_velocity.get_x_coord();}
    inline const DT get_velocity_y_coord() const {return m_velocity.get_y_coord();}
    inline const DT get_velocity_z_coord() const {return m_velocity.get_z_coord();}
    
    inline const V3<DT> get_acceleration() const {return m_acceleration;}
    inline const DT get_acceleration_x_coord() const {return m_acceleration.get_x_coord();}
    inline const DT get_acceleration_y_coord() const {return m_acceleration.get_y_coord();}
    inline const DT get_acceleration_z_coord() const {return m_acceleration.get_z_coord();}
    
    inline const V3<DT> get_next_pos() const {return m_next_pos;}
    inline const DT get_next_pos_x_coord() const {return m_next_pos.get_x_coord();}
    inline const DT get_next_pos_y_coord() const {return m_next_pos.get_y_coord();}
    inline const DT get_next_pos_z_coord() const {return m_next_pos.get_z_coord();}
    
    inline const V3<DT> get_next_vel() const {return m_next_vel;}
    inline const DT get_next_vel_x_coord() const {return m_next_vel.get_x_coord();}
    inline const DT get_next_vel_y_coord() const {return m_next_vel.get_y_coord();}
    inline const DT get_next_vel_z_coord() const {return m_next_vel.get_z_coord();}
    
    void set_mass (const DT& mass) {m_mass=mass;}
    void set_position(const DT& x, const DT& y, const DT& z) {m_position.set_coord(x,y,z);}
    void set_velocity(const DT& x, const DT& y, const DT& z) {m_velocity.set_coord(x,y,z);}
    void set_acceleration(const DT& x, const DT& y, const DT& z) {m_acceleration.set_coord(x,y,z);}
    
    void set_position_by_vector(const V3<DT>& to_assign) {m_position = to_assign;}
    void set_velocity_by_vector(const V3<DT>& to_assign) {m_velocity = to_assign;}
    void set_acceleration_by_vector(const V3<DT>& to_assign) {m_acceleration = to_assign;}
    
    void add_acceleration_contribution(const V3<DT>& accel_to_add) {m_acceleration+=accel_to_add;}
    
    void set_next_pos(const DT& x, const DT& y, const DT& z) {m_next_pos.set_coord(x,y,z);}
    void set_next_vel(const DT& x, const DT& y, const DT& z) {m_next_vel.set_coord(x,y,z);}
    void set_next_pos_by_vector(const V3<DT>& to_assign) {m_next_pos = to_assign;}
    void set_next_vel_by_vector(const V3<DT>& to_assign) {m_next_vel = to_assign;}
private:
    DT m_mass;
    V3<DT> m_position;
    V3<DT> m_velocity;
    V3<DT> m_acceleration;
    //For Runge-Kutta implementation
    V3<DT> m_next_pos;
    V3<DT> m_next_vel;
};

#endif //BODY_CLASS_FOR_N_BODY_SIMULATION
