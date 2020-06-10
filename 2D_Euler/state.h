#include "constants.h"

/* Class State

A State contains 4 doubles corresponding to the 4 conserved physical quantities:
*rho    : density
*mom_x  : momentum in the x direction
*mom_y  : momentum in the y direction
*energy : Total energy (rho.E)


*/

class State {

public:
    /* Variables*/
    double rho, mom_x, mom_y, energy;

    /*constructor*/
    State(){};

    /*operations*/
    State operator+(State);
    State operator-(State);
    State operator*(double);

    /* pres : computes and returns the internal pressure */
    double pres()
    {
        return (gr_gamma - 1.0f) * (energy - 0.5 * (mom_x * mom_x + mom_y * mom_y) / rho);
    };

    /* xvel : computes and returns the x velocity */
    double xvel()
    {
        return mom_x / rho;
    };

    /* yvel : computes and returns the y velocity */
    double yvel()
    {
        return mom_y / rho;
    };

    /* Computes and returns the sound speed
     When we need SS, we usually have computed the pressure before */
    double sound_speed(double pressure)
    {
        return sqrt(gr_gamma * pressure / rho);
    };

    /* compute the time step in a cell according to the CFL condition */
    double self_dt(double CFL, double dx, double dy)
    {
        double vx, vy, c, pressure, carac_time;

        vx = xvel();
        vy = yvel();
        pressure = pres();
        c = sound_speed(pressure);

        carac_time = min(dx / max(abs(vx - c), abs(vx + c)), dy / max(abs(vy - c), abs(vy + c)));

        return CFL * carac_time;
    };

    /* fills the State with 4 givent inputs */
    void fill(double rho_in, double mom_x_in, double mom_y_in, double energy_in)
    {
        rho = rho_in;
        mom_x = mom_x_in;
        mom_y = mom_y_in;
        energy = energy_in;
    }

    /* set the States to u_in */
    void affect(State u_in)
    {
        rho = u_in.rho;
        mom_x = u_in.mom_x;
        mom_y = u_in.mom_y;
        energy = u_in.energy;
    };

    /* Converts the State from p to c*/
    void primitive_to_conservative()
    {

        double vx, vy, pressure;

        vx = mom_x;
        vy = mom_y;
        pressure = energy;

        mom_x = rho * vx;
        mom_y = rho * vy;
        energy = (pressure / (gr_gamma - 1.0f) + 0.5 * rho * (vx * vx + vy * vy));

        /* u(ener)= ( (v(4)/(y-1.)) + 0.5*v(1)*(v(2)**2+v(3)**2) ) */
    }

    /* Converts the State from c to p*/

    void conservative_to_primitive()
    {

        double vx, vy, pressure;
        vx = xvel();
        vy = yvel();
        pressure = pres();

        mom_x = vx;
        mom_y = vy;
        energy = pressure;
    };
};

State State::operator+(State param)
{
    State temp;
    temp.rho = rho + param.rho;
    temp.mom_x = mom_x + param.mom_x;
    temp.mom_y = mom_x + param.mom_y;
    temp.energy = energy + param.energy;

    return (temp);
};

State State::operator-(State param)
{
    State temp;
    temp.rho = rho - param.rho;
    temp.mom_x = mom_x - param.mom_x;
    temp.mom_y = mom_x - param.mom_y;
    temp.energy = energy - param.energy;

    return (temp);
};

State State::operator*(double param)
{
    State temp;
    temp.rho = rho * param;
    temp.mom_x = mom_x * param;
    temp.mom_y = mom_x * param;
    temp.energy = energy * param;

    return (temp);
};
