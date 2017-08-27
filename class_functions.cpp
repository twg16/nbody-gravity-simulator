// Thomas Goode - 100759942

#include "class_dec.h"

//////// Vector Functions //////////
Vector::Vector() : _x(0.0), _y(0.0), _z(0.0)
{
	parameters = new double();
}
Vector::Vector(double x, double y, double z) : _x(x), _y(y), _z(z)
{
	parameters = new double();
} 
Vector::Vector(const Vector &rhs)
{
	_x = rhs._x;
	_y = rhs._y;
	_z = rhs._z;
	parameters = new double();
	*parameters = *(rhs.parameters);
}
Vector::~Vector()
{
	delete parameters;
}

double Vector::getX()
{
	return _x;
}
double Vector::getY()
{
	return _y;
}
double Vector::getZ()
{
	return _z;
}
double Vector::getMod()
{
	return sqrt((_x * _x) + (_y * _y) + (_z * _z));
}
Vector Vector::multiply(double a)
{
	double X = _x * a;
	double Y = _y * a;
	double Z = _z * a;
	return Vector(X, Y, Z);
}

//overload addition operator
Vector Vector::operator+(const Vector &rhs)
{
	Vector temp;
	temp._x = _x + rhs._x;
	temp._y = _y + rhs._y;
	temp._z = _z + rhs._z;
	*(temp.parameters) = *parameters + *(rhs.parameters);
	return temp;
}

//overload subtraction operator
Vector Vector::operator-(const Vector &rhs)
{
	Vector temp;
	temp._x = _x - rhs._x;
	temp._y = _y - rhs._y;
	temp._z = _z - rhs._z;
	*(temp.parameters) = *parameters - *(rhs.parameters);
	return temp;
}
//overload assignment operator
Vector Vector::operator=(const Vector &rhs)
{
	if (this != &rhs)
	{
		_x = rhs._x;
		_y = rhs._y;
		_z = rhs._z;
		*parameters = *(rhs.parameters);
	}
	return *this;
}
//overload plus-equals operator
Vector Vector::operator+=(const Vector &rhs)
{
	if (this != &rhs)
	{
		_x += rhs._x;
		_y += rhs._y;
		_z += rhs._z;
		*parameters += *(rhs.parameters);
	}
	return *this;
}
////////////body functions///////////////
double Body::getM()
{
	return _m;
}
Vector Body::getPos()
{
	return _pos;
}
Vector Body::getVel()
{
	return _vel;
}

//////Euler functions//////
Euler::Euler()
{
}
Euler::~Euler()
{

}
Vector Euler::calcPosition(Vector POS, Vector VEL, double time)
{
	double xn = POS.getX() + (VEL.getX()*time);
	double yn = POS.getY() + (VEL.getY()*time);
	double zn = POS.getZ() + (VEL.getZ()*time);

	return Vector(xn, yn, zn);
}

Vector Euler::calcVelocity(Vector VEL, Vector FORCE, double time, double mass)
{
	double xvn = VEL.getX() + ((FORCE.getX() / mass)*time);
	double yvn = VEL.getY() + ((FORCE.getY() / mass)*time);
	double zvn = VEL.getZ() + ((FORCE.getZ() / mass)*time);

	return Vector(xvn, yvn, zvn);
}
/////Leapfrog Functions/////
Leapfrog::Leapfrog()
{
}
Leapfrog::~Leapfrog()
{
}
Vector Leapfrog::calcPosition(Vector POS, Vector VEL, Vector FORCE, double mass, double time)
{
	double xn = POS.getX() + (VEL.getX()*time) + (0.5*(FORCE.getX() / mass)*time*time);
	double yn = POS.getY() + (VEL.getY()*time) + (0.5*(FORCE.getY() / mass)*time*time);
	double zn = POS.getZ() + (VEL.getZ()*time) + (0.5*(FORCE.getZ() / mass)*time*time);

	return Vector(xn, yn, zn);
}
Vector Leapfrog::calcVelocity(Vector VEL, Vector FORCE, Vector newFORCE, double time, double mass)
{
	double xvn = VEL.getX() + (0.5* ((FORCE.getX() / mass) + (newFORCE.getX() / mass)) * time);
	double yvn = VEL.getY() + (0.5* ((FORCE.getY() / mass) + (newFORCE.getY() / mass)) * time);
	double zvn = VEL.getZ() + (0.5* ((FORCE.getZ() / mass) + (newFORCE.getZ() / mass)) * time);

	return Vector(xvn, yvn, zvn);
}

/////////////Runge Kutta Functions/////////////
RungeKutta::RungeKutta()
{
}
RungeKutta::~RungeKutta()
{
}
Vector RungeKutta::calcK(Vector FORCE, double mass)
{	
	double Kx = ((FORCE.getX() / mass));
	double Ky = ((FORCE.getY() / mass));
	double Kz = ((FORCE.getZ() / mass));
	return Vector(Kx, Ky, Kz);
}
Vector RungeKutta::calcK1(Vector VEL)
{
	double Kx = VEL.getX();
	double Ky = VEL.getY();
	double Kz = VEL.getZ();
	return Vector(Kx, Ky, Kz);
}
Vector RungeKutta::calcPosition(Vector POS, Vector k1, Vector k2, Vector k3, Vector k4)
{
	Vector k2_2 = k2.multiply(2);
	Vector k3_2 = k3.multiply(2);
	Vector k_term = ((k1.operator+(k2_2)).operator+(k3_2)).operator+(k4);
	Vector t_k_term = k_term.multiply(4 * 3600);
	double xn = POS.getX() + t_k_term.getX();
	double yn = POS.getY() + t_k_term.getY();
	double zn = POS.getZ() + t_k_term.getZ();
	return Vector(xn, yn, zn);
}

Vector RungeKutta::calcVelocity(Vector VEL, Vector k1, Vector k2, Vector k3, Vector k4)
{
	Vector k2_2 = k2.multiply(2);
	Vector k3_2 = k3.multiply(2);
	Vector k_term = ((k1.operator+(k2_2)).operator+(k3_2)).operator+(k4);
	Vector t_k_term = k_term.multiply(4 * 3600);
	double xn = VEL.getX() + t_k_term.getX();
	double yn = VEL.getY() + t_k_term.getY();
	double zn = VEL.getZ() + t_k_term.getZ();
	return Vector(xn, yn, zn);
}

Vector RungeKutta::getVelocity(Vector k1, Vector k2, Vector k3, Vector k4, double time)
{
	Vector k2_2 = k2.multiply(2);
	Vector k3_2 = k3.multiply(2);
	Vector k_term = ((k1.operator+(k2_2)).operator+(k3_2)).operator+(k4);
	Vector vel_term = k_term.multiply(0.1666666667);
	Vector vel_term1 = vel_term.multiply(time);
	return vel_term1;
}







