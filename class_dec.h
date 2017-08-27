// Thomas Goode - 100759942

#include <iostream>

///////// Vector class ///////////
class Vector {
public:
	Vector();
	Vector(double x, double y, double z);
	Vector(const Vector &rhs); 

	~Vector();

	double getX();
	double getY();
	double getZ();
	double getMod();
	Vector multiply(double a);
	
	//overload addition operator
	Vector Vector::operator+(const Vector &rhs);
	//overload subtraction operator
	Vector Vector::operator-(const Vector &rhs);
	//overload plus-equals operator
	Vector Vector::operator+=(const Vector &rhs);
	//overload assignment operator
	Vector Vector::operator=(const Vector &rhs);

private:
	double _x;
	double _y;
	double _z;
	double *parameters;
};

///////////////Body class///////////////
class Body 
{
public:
	Body(double m, Vector pos, Vector vel) 
	{
		_m = m;
		_pos = pos;
		_vel = vel;
	}
	double getM();
	Vector getPos();
	Vector getVel();

private:
	double _m;
	Vector _pos;
	Vector _vel;

};

///////// Euler class /////////
class Euler
{
public:
	Euler();
	~Euler();
	Vector calcPosition(Vector POS, Vector VEL, double time);
	Vector calcVelocity(Vector VEL, Vector FORCE, double time, double mass);
};

////////// Leapfrog class //////
class Leapfrog
{
public:
	Leapfrog();
	~Leapfrog();
	Vector calcPosition(Vector POS, Vector VEL, Vector FORCE, double mass, double time);
	Vector calcVelocity(Vector VEL, Vector FORCE, Vector newFORCE, double time, double mass);
};

/////// Runge-Kutta class ///////////
class RungeKutta
{
public:
	RungeKutta();
	~RungeKutta();
	Vector calcK(Vector FORCE, double mass);
	Vector calcK1(Vector VEL);
	Vector calcPosition(Vector POS, Vector k1, Vector k2, Vector k3, Vector k4);
	Vector calcVelocity(Vector VEL, Vector k1, Vector k2, Vector k3, Vector k4);
	Vector getVelocity(Vector k1, Vector k2, Vector k3, Vector k4, double time);
};
