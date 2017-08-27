// Thomas Goode - 100759942

#include <iostream>
#include "class_dec.h"
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
//defining the objects required to implement each method
Euler euler_method;
Leapfrog lf_method;
RungeKutta rk_method;

double AU = 149597870700;
double G = 6.67 * pow(10, -11);

//defining bodies and variables for the 2-body problem
Vector sunPos(0, 0, 0);
Vector earthPos((1*AU), 0, 0);
Vector sunVel(0, 0, 0);
Vector earthVel(0, 15000, 0); //29783,
Body Sun1((1.98892 * pow(10, 30)), sunPos, sunVel);
Body Earth1((5.9742 * pow(10, 24)), earthPos, earthVel);

//defining the bodies for the N-body  solar system problem
Body Sun((1.98892 * pow(10, 30)), Vector(0, 0, 0), Vector(0, 0, 0));
Body Mercury((0.33 * pow(10, 24)), Vector((0.387*AU), 0, 0), Vector(0, 47400, 0));
Body Venus((4.87 * pow(10, 24)), Vector((0.723*AU), 0, 0), Vector(0, 35000, 0));
Body Earth((5.98 * pow(10, 24)), Vector((1 * AU), 0, 0), Vector(0, 29800, 0));
Body Mars((6.42*pow(10, 23)), Vector((1.524*AU), 0, 0), Vector(0, 24100, 0));
Body Jupiter((1.90*pow(10, 27)), Vector((5.203*AU), 0, 0), Vector(0, 13100, 0));
Body Saturn((5.69*pow(10, 26)), Vector((9.582*AU), 0, 0), Vector(0, 9600, 0));
Body Uranus((8.68*pow(10, 25)), Vector((19.20*AU), 0, 0), Vector(0, 6800, 0));
Body Neptune((1.03*pow(10, 26)), Vector((30.05*AU), 0, 0), Vector(0, 5400, 0));

//defining the bodies for the circumbinary orbit problem
Body StarA((4.988 * pow(10, 30)), Vector(0, 0, 0), Vector(0, -5000, 0));
Body StarB((4.988*pow(10, 30)), Vector((15.603 * AU), 0, 0), Vector(0, 5000, 0));
Body Planet((6.42*pow(10, 27)), Vector(0, (-25 * AU), 0), Vector(16000, 0, 0));

//function to calculate force for the 2-body problem
Vector calcForce(Body A, Body B)
{
	Vector R = A.getPos().operator-(B.getPos());	//Vector of the distance between the two bodies
	double r = R.getMod();
	double rrr = pow(r, 3);
	double f = -(G * A.getM() * B.getM()) / (rrr);	//force as a scalr
	return R.multiply(f);							//multiply the two to get a force Vector
}
//function to calculate force for the N-body problem
Vector nBodyCalcForce(int bodyNumber, Body body, std::vector<Body> bodies)
{
	double X = 0;					//defines the variables X, Y and Z
	double Y = 0;
	double Z = 0;					//the for statment iterates through a vector of all the bodies in a system
	for (int j = 0; j <= (bodies.size()-1); j++)
	{
		if (j == bodyNumber)			
		{							//if the iteration number is equal to the number of the body in the 
			continue;				//vector then it skips the iteration using 'continue'
		}
		Vector R = body.getPos().operator-(bodies[j].getPos());	//calculate the distance between the two bodies in question
		double r = R.getMod();
		double rrr = pow(r, 3);
		double f = -(G * body.getM() * bodies[j].getM()) / (rrr);		//calculates the force term as a scalar...
		Vector FORCE = R.multiply(f);					//... then the distance vector is multiplied by the scalar
		X += FORCE.getX();				//each iteration adds to the parameters X, Y and Z
		Y += FORCE.getY();
		Z += FORCE.getZ();
	}
	return Vector(X, Y, Z);	
}

//function to calculate k when manipulating the v=v+at equation for the runge kutta method
Vector calcKForce(Vector kn, Vector kPosE, Vector kVelE, double time)
{
	Vector velStepk = kn.multiply(time);	//calculate a step forward in velocity
	kVelE.operator+=(velStepk);
	Body k_newEarth = Body(Earth1.getM(), kPosE, kVelE);	//update the bodies
	Body k_sun = Body(Sun1.getM(), Sun1.getPos(), Sun1.getVel());
	Vector knplus1Force = calcForce(k_newEarth, k_sun);			//calculate a new force to use to calculate the next k
	Vector knplus1 = rk_method.calcK(knplus1Force, Earth1.getM());		//use the calcK function in the runge kutta class to calculate the next k
	Vector negVelStepk = velStepk.multiply(-1); 	//reset the velocity parameter
	kVelE.operator+=(negVelStepk);
	return knplus1;					//return the next k as a Vector (which in this case is an acceleration)
}
//function to calculate k when manipulating the u=u+vt equation for the runge kutta method
Vector calcKVelocity(Vector kn, Vector kPosE, Vector kVelE, double time)
{
	Vector posStepkv = kn.multiply(time);		//calculate a step forward in position
	kPosE.operator+=(posStepkv);		
	Body k_newEarth = Body(Earth1.getM(), kPosE, kVelE);
	Body k_sun = Body(Sun1.getM(), Sun1.getPos(), Sun1.getVel());
	Vector knplus1 = rk_method.calcK1(kVelE);			//use the K1 function in the runge kutta class to calculate the next K
	Vector negPosStepkv = posStepkv.multiply(-1);		//reset the position parameter
	kPosE.operator+=(negPosStepkv);
	return knplus1;
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
////////////////////////////// main /////////////////////////////////////////////

int main()
{
	//defining the parameters that let the user chose their own intial conditions
	int simulation;
	int initialVelocity;
	int initialPosition;
	int duration;
	int simLength1 = 0;

	//the next part of the code is used to constrcut a user interface
	std::cout << "Please enter number to run the" <<
		" corresponding simulation: " << std::endl;
	std::cout << "\t" << "1 - Orbit of a body with a mass of the Earth around the Sun [euler vs leapfrog vs Runge-Kutta]" << std::endl;
	std::cout << "\t" << "2 - Orbits of the planets in the solar system [Runge-Kutta]" << std::endl;
	std::cout << "\t" << "3 - Circumbinary Orbit [leapfrog]" << std::endl;
	std::cin >> simulation;
	
	if (simulation == 1)
	{	
		std::cout << "Choose an initial tangential velocity of the body orbiting the Sun: " << std::endl;
		std::cout << "\t" << "1 (15,000 m/s)" << std::endl;
		std::cout << "\t" << "2 (29,783 m/s) [Earth's initial velocity]" << std::endl;
		std::cout << "\t" << "3 (40,000 m/s)" << std::endl;

		std::cin >> initialVelocity;
		if (initialVelocity == 1)
		{
			earthVel = Vector(0, 15000, 0);
		}
		if (initialVelocity == 2)
		{
			earthVel = Vector(0, 29873, 0);
		}
		if (initialVelocity == 3)
		{
			earthVel = Vector(0, 40000, 0);
		}
		std::cout << "Choose an initial position for the body orbiting the Sun:" << std::endl;
		std::cout << "\t" << "1 (0.45 AU) [Between orbits of Mercury and Venus]" << std::endl;
		std::cout << "\t" << "2 (1 AU) [Orbit of the Earth]" << std::endl;
		std::cout << "\t" << "3 (5.203 AU) [Orbit of Jupiter]" << std::endl;
		std::cin >> initialPosition;

		if (initialPosition == 1)
		{
			earthPos = Vector((0.45*AU), 0, 0);
		}
		if (initialPosition == 2)
		{
			earthPos = Vector((1 * AU), 0, 0);
		}
		if (initialPosition == 3)
		{
			earthPos = Vector((5.203*AU), 0, 0);
		}
		std::cout << "Choose a duration of time to run the simulation for:" << std::endl;
		std::cout << "\t" << "1 (1 year)" << std::endl;
		std::cout << "\t" << "2 (5 years)" << std::endl;
		std::cout << "\t" << "3 (20 years)" << std::endl;
		std::cin >> duration;
		if (duration == 1)
		{
			simLength1 = (1 * 365);
		}
		if (duration == 2)
		{
			simLength1 = (5 * 365);
		}
		if (duration == 3)
		{
			simLength1 = (20 * 365);
		}

		///////////2-Body Problem - Euler vs leapfrog vs Runge-Kutta////////
		Vector forceE = calcForce(Earth1, Sun1);
		Vector forceS = calcForce(Sun1, Earth1);

		Vector euler_tempPosE = earthPos;
		Vector euler_tempVelE = earthVel;
		Vector euler_tempPosS = sunPos;
		Vector euler_tempVelS = sunVel;
		Vector euler_tempForceE = forceE;
		Vector euler_tempForceS = forceS;

		Vector lf_tempPosE = earthPos;
		Vector lf_tempVelE = earthVel;
		Vector lf_tempPosS = sunPos;
		Vector lf_tempVelS = sunVel;
		Vector lf_tempForceE = forceE;
		Vector lf_tempForceS = forceS;

		Vector rk_tempPosE = earthPos;
		Vector rk_tempVelE = earthVel;
		Vector rk_tempPosS = sunPos;
		Vector rk_tempVelS = sunVel;
		Vector rk_tempForceE = forceE;
		Vector rk_tempForceS = forceS;
		Vector kPosE = earthPos;
		Vector kVelE = earthVel;
		Body k_earth(Earth1.getM(), kPosE, kVelE);
		Body k_newEarth(Earth1.getM(), kPosE, kVelE);
		Vector kPosS = sunPos;
		Vector kVelS = sunVel;
		Body k_sun(Sun1.getM(), kPosS, kVelS);


		std::ofstream EulerMyFileEarth;
		std::ofstream EulerMyFileSun;
		EulerMyFileEarth.open("euler_earth_data1.txt");
		EulerMyFileSun.open("euler_sun_data1.txt");

		std::ofstream LfMyFileEarth;
		std::ofstream LfMyFileSun;
		LfMyFileEarth.open("lf_earth_data1.txt");
		LfMyFileSun.open("lf_sun_data1.txt");

		std::ofstream RkMyFileEarth;
		std::ofstream RkMyFileSun;
		RkMyFileEarth.open("rk_earth_data1.txt");
		RkMyFileSun.open("rk_sun_data1.txt");
	
		std::cout << "Writing files..." << std::endl;
		for (int i = 0; i <= simLength1; i++)
		{
			//////// Euler Method ////////////
			
			Body euler_newEarth(Earth1.getM(), euler_tempPosE, euler_tempVelE);
			Body euler_newSun(Sun1.getM(), euler_tempPosS, euler_tempVelS);
			
			euler_tempPosE = euler_method.calcPosition(euler_tempPosE, euler_tempVelE, (24 * 3600));
			euler_tempVelE = euler_method.calcVelocity(euler_tempVelE, euler_tempForceE, (24 * 3600), Earth1.getM());
			euler_tempPosS = euler_method.calcPosition(euler_tempPosS, euler_tempVelS, (24 * 3600));
			euler_tempVelS = euler_method.calcVelocity(euler_tempVelS, euler_tempForceS, (24 * 3600), Sun1.getM());
			
			euler_tempForceE = calcForce(euler_newEarth, euler_newSun);
			euler_tempForceS = calcForce(euler_newSun, euler_newEarth);
		
			EulerMyFileEarth << euler_tempPosE.getX() << "\t" << euler_tempPosE.getY() << "\t"
				<< euler_tempPosE.getZ() << std::endl;

			EulerMyFileSun << euler_tempPosS.getX() << "\t" << euler_tempPosS.getY() << "\t"
				<< euler_tempPosS.getZ() << std::endl;

			/////////////////////////////////////////////////////////////

			///////// Leapfrog Method ///////////////////////
			
			lf_tempPosE = lf_method.calcPosition(lf_tempPosE, lf_tempVelE, lf_tempForceE, Earth1.getM(), (24 * 3600));
			lf_tempPosS = lf_method.calcPosition(lf_tempPosS, lf_tempVelS, lf_tempForceS, Sun1.getM(), (24 * 3600));
			
			Body lf_newEarth(Earth1.getM(), lf_tempPosE, lf_tempVelE);
			Body lf_newSun(Sun1.getM(), lf_tempPosS, lf_tempVelS);
			
			Vector nextFORCEE = calcForce(lf_newEarth, lf_newSun);
			Vector nextFORCES = calcForce(lf_newSun, lf_newEarth);
			
			lf_tempVelE = lf_method.calcVelocity(lf_tempVelE, lf_tempForceE, nextFORCEE, (24 * 3600), Earth1.getM());
			lf_tempVelS = lf_method.calcVelocity(lf_tempVelS, lf_tempForceS, nextFORCES, (24 * 3600), Sun1.getM());
			
			lf_newEarth = Body(Earth1.getM(), lf_tempPosE, lf_tempVelE);
			lf_newSun = Body(Sun1.getM(), lf_tempPosS, lf_tempVelS);
			
			lf_tempForceE = calcForce(lf_newEarth, lf_newSun);
			lf_tempForceS = calcForce(euler_newSun, euler_newEarth);

			LfMyFileEarth << lf_tempPosE.getX() << "\t" << lf_tempPosE.getY() << "\t"
				<< lf_tempPosE.getZ() << std::endl;

			LfMyFileSun << lf_tempPosS.getX() << "\t" << lf_tempPosS.getY() << "\t"
				<< lf_tempPosS.getZ() << std::endl;

			/////////////Runge-Kutta Method////////////////////
			
			//setting new bodies
			Body rk_newEarth = Body(Earth1.getM(), rk_tempPosE, rk_tempVelE);
			Body rk_newSun = Body(Sun1.getM(), rk_tempPosS, rk_tempVelS);

			//calculate k1, k2, k3, k4
			Vector k1Force = calcForce(k_newEarth, k_sun);
			Vector k1 = rk_method.calcK(k1Force, Earth1.getM());
			Vector k2 = calcKForce(k1, kPosE, kVelE, (12 * 3600));
			Vector k3 = calcKForce(k2, kPosE, kVelE, (12 * 3600));
			Vector k4 = calcKForce(k3, kPosE, kVelE, (24 * 3600));
			//using the get velocity function from the runge kutta class to use in the next equation
			Vector rk_vel = rk_method.getVelocity(k1, k2, k3, k4, (24*3600));
			//add this velocity  change to the velocity used to calculate K 
			kVelE.operator+=(rk_vel);
			//calculate K1, K2, K3, K4
			Vector K1 = rk_method.calcK1(kVelE);
			Vector K2 = calcKVelocity(K1, kPosE, kVelE, (12 * 3600));
			Vector K3 = calcKVelocity(K2, kPosE, kVelE, (12 * 3600));
			Vector K4 = calcKVelocity(K3, kPosE, kVelE, (24 * 3600));

			//calculating positions and velocities
			rk_tempVelE = rk_method.calcVelocity(rk_tempVelE, k1, k2, k3, k4);
			rk_tempPosE = rk_method.calcPosition(rk_tempPosE, K1, K2, K3, K4);
			//resetting the bodies
			rk_newEarth = Body(Earth1.getM(), rk_tempPosE, rk_tempVelE);
			rk_newSun = Body(Sun1.getM(), rk_tempPosS, rk_tempVelS);
			k_newEarth = Body(Earth1.getM(), rk_tempPosE, rk_tempVelE);
			kPosE = rk_tempPosE;
			kVelE = rk_tempVelE;

			
			RkMyFileEarth << rk_tempPosE.getX() << "\t" << rk_tempPosE.getY() << "\t"
				<< rk_tempPosE.getZ() << std::endl;

			RkMyFileSun << rk_tempPosS.getX() << "\t" << rk_tempPosS.getY() << "\t"
				<< rk_tempPosS.getZ() << std::endl;

		}
	
		EulerMyFileEarth.close();
		EulerMyFileSun.close();

		LfMyFileEarth.close();
		LfMyFileSun.close();

		RkMyFileEarth.close();
		RkMyFileSun.close();
	}
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	
	if (simulation == 2)
	{
		////////// N-Body Problem //////////////////////
		std::ofstream RkNBody;
		RkNBody.open("rk_nbody_data.txt");
		
		std::vector<Body> bodies;
		std::vector<Body>::iterator it;

		Vector kPosMercury = Vector((0.387*AU), 0, 0);
		Vector kVelMercury = Vector(0, 47400, 0);
		Body k_mercury(Mercury.getM(), kPosMercury, kVelMercury);
		Body k_newMercury(Mercury.getM(), kPosMercury, kVelMercury);

		Vector kPosVenus = Vector((0.723*AU), 0, 0);
		Vector kVelVenus = Vector(0, 35000, 0);
		Body k_venus(Venus.getM(), kPosVenus, kVelVenus);
		Body k_newVenus(Venus.getM(), kPosVenus, kVelVenus);

		Vector kPosEarth = Vector((1*AU), 0, 0);
		Vector kVelEarth = Vector(0, 29800, 0);
		Body k_earth(Earth.getM(), kPosEarth, kVelEarth);
		Body k_newEarth(Earth.getM(), kPosEarth, kVelEarth);

		Vector kPosMars = Vector((1.524*AU), 0, 0);
		Vector kVelMars = Vector(0, 24100, 0);
		Body k_mars(Mars.getM(), kPosMars, kVelMars);
		Body k_newMars(Mars.getM(), kPosMars, kVelMars);

		Vector kPosJupiter = Vector((5.203*AU), 0, 0);
		Vector kVelJupiter = Vector(0, 13100, 0);
		Body k_jupiter(Jupiter.getM(), kPosJupiter, kVelJupiter);
		Body k_newJupiter(Jupiter.getM(), kPosJupiter, kVelJupiter);

		Vector kPosSaturn = Vector((9.582*AU), 0, 0);
		Vector kVelSaturn = Vector(0, 9600, 0);
		Body k_saturn(Saturn.getM(), kPosSaturn, kVelSaturn);
		Body k_newSaturn(Saturn.getM(), kPosSaturn, kVelSaturn);

		Vector kPosUranus = Vector((19.20*AU), 0, 0);
		Vector kVelUranus = Vector(0, 6800, 0);
		Body k_uranus(Uranus.getM(), kPosUranus, kVelUranus);
		Body k_newUranus(Uranus.getM(), kPosUranus, kVelUranus);

		Vector kPosNeptune = Vector((30.05*AU), 0, 0);
		Vector kVelNeptune = Vector(0, 5400, 0);
		Body k_neptune(Neptune.getM(), kPosNeptune, kVelNeptune);
		Body k_newNeptune(Neptune.getM(), kPosNeptune, kVelNeptune);

		it = bodies.begin();
		bodies.push_back(Sun);
		bodies.push_back(Mercury);
		bodies.push_back(Venus);
		bodies.push_back(Earth);
		bodies.push_back(Mars);
		bodies.push_back(Jupiter);
		bodies.push_back(Saturn);
		bodies.push_back(Uranus);
		bodies.push_back(Neptune);

		Vector nBody_tempPosMercury = Mercury.getPos();
		Vector nBody_tempVelMercury = Mercury.getVel();
		Vector nBody_tempForceMercury = nBodyCalcForce(1, Mercury, bodies);

		Vector nBody_tempPosVenus = Venus.getPos();
		Vector nBody_tempVelVenus = Venus.getVel();
		Vector nBody_tempForceVenus = nBodyCalcForce(2, Venus, bodies);

		Vector nBody_tempPosEarth = Earth.getPos();
		Vector nBody_tempVelEarth = Earth.getVel();
		Vector nBody_tempForceEarth = nBodyCalcForce(3, Earth, bodies);

		Vector nBody_tempPosMars = Mars.getPos();
		Vector nBody_tempVelMars = Mars.getVel();
		Vector nBody_tempForceMars = nBodyCalcForce(4, Mars, bodies);

		Vector nBody_tempPosJupiter = Jupiter.getPos();
		Vector nBody_tempVelJupiter = Jupiter.getVel();
		Vector nBody_tempForceJupiter = nBodyCalcForce(5, Jupiter, bodies);

		Vector nBody_tempPosSaturn = Saturn.getPos();
		Vector nBody_tempVelSaturn = Saturn.getVel();
		Vector nBody_tempForceSaturn = nBodyCalcForce(6, Saturn, bodies);

		Vector nBody_tempPosUranus = Uranus.getPos();
		Vector nBody_tempVelUranus = Uranus.getVel();
		Vector nBody_tempForceUranus = nBodyCalcForce(7, Uranus, bodies);

		Vector nBody_tempPosNeptune = Neptune.getPos();
		Vector nBody_tempVelNeptune = Neptune.getVel();
		Vector nBody_tempForceNeptune = nBodyCalcForce(8, Neptune, bodies);


		std::cout << "Writing files..." << std::endl;
		for (int k = 0; k <= 200 * 365; k++)
		{
			//setting new bodies
			Body rk_newSun = Body(Sun.getM(), Vector(0, 0, 0), Vector(0, 0, 0));
			Body rk_newMercury = Body(Mercury.getM(), nBody_tempPosMercury, nBody_tempVelMercury);
			Body rk_newVenus = Body(Venus.getM(), nBody_tempPosVenus, nBody_tempVelVenus);
			Body rk_newEarth = Body(Earth1.getM(), nBody_tempPosEarth, nBody_tempVelEarth);
			Body rk_newMars = Body(Mars.getM(), nBody_tempPosMars, nBody_tempVelMars);
			Body rk_newJupiter = Body(Jupiter.getM(), nBody_tempPosJupiter, nBody_tempVelJupiter);
			Body rk_newSaturn = Body(Saturn.getM(), nBody_tempPosSaturn, nBody_tempVelSaturn);
			Body rk_newUranus = Body(Uranus.getM(), nBody_tempPosUranus, nBody_tempVelUranus);
			Body rk_newNeptune = Body(Neptune.getM(), nBody_tempPosNeptune, nBody_tempVelNeptune);


			//calculate k1 force
			Vector mercuryk1Force = nBodyCalcForce(1, rk_newMercury, bodies);
			Vector venusk1Force = nBodyCalcForce(2, rk_newVenus, bodies);
			Vector earthk1Force = nBodyCalcForce(3, rk_newEarth, bodies);
			Vector marsk1Force = nBodyCalcForce(4, rk_newMars, bodies);
			Vector jupiterk1Force = nBodyCalcForce(5, rk_newJupiter, bodies);
			Vector saturnk1Force = nBodyCalcForce(6, rk_newSaturn, bodies);
			Vector uranusk1Force = nBodyCalcForce(7, rk_newUranus, bodies);
			Vector neptunek1Force = nBodyCalcForce(8, rk_newNeptune, bodies);

			//MERCURY
			Vector k1mercury = rk_method.calcK(mercuryk1Force, Mercury.getM());
			Vector k2mercury = calcKForce(k1mercury, kPosMercury, kVelMercury, (12 * 3600));
			Vector k3mercury = calcKForce(k2mercury, kPosMercury, kVelMercury, (12 * 3600));
			Vector k4mercury = calcKForce(k3mercury, kPosMercury, kVelMercury, (24 * 3600));
			Vector rk_velMercury = rk_method.getVelocity(k1mercury, k2mercury, k3mercury, k4mercury, (24 * 3600));
			kVelMercury.operator+=(rk_velMercury);
			Vector K1mercury = rk_method.calcK1(kVelMercury);
			Vector K2mercury = calcKVelocity(K1mercury, kPosMercury, kVelMercury, (12 * 3600));
			Vector K3mercury = calcKVelocity(K2mercury, kPosMercury, kVelMercury, (12 * 3600));
			Vector K4mercury = calcKVelocity(K3mercury, kPosMercury, kVelMercury, (24 * 3600));
			//VENUS
			Vector k1venus = rk_method.calcK(venusk1Force, Venus.getM());
			Vector k2venus = calcKForce(k1venus, kPosVenus, kVelVenus, (12 * 3600));
			Vector k3venus = calcKForce(k2venus, kPosVenus, kVelVenus, (12 * 3600));
			Vector k4venus = calcKForce(k3venus, kPosVenus, kVelVenus, (24 * 3600));
			Vector rk_velVenus = rk_method.getVelocity(k1venus, k2venus, k3venus, k4venus, (24 * 3600));
			kVelVenus.operator+=(rk_velVenus);
			Vector K1venus = rk_method.calcK1(kVelVenus);
			Vector K2venus = calcKVelocity(K1venus, kPosVenus, kVelVenus, (12 * 3600));
			Vector K3venus = calcKVelocity(K2venus, kPosVenus, kVelVenus, (12 * 3600));
			Vector K4venus = calcKVelocity(K3venus, kPosVenus, kVelVenus, (24 * 3600));
			//EARTH
			Vector k1earth = rk_method.calcK(earthk1Force, Earth.getM());
			Vector k2earth = calcKForce(k1earth, kPosEarth, kVelEarth, (12 * 3600));
			Vector k3earth = calcKForce(k2earth, kPosEarth, kVelEarth, (12 * 3600));
			Vector k4earth = calcKForce(k3earth, kPosEarth, kVelEarth, (24 * 3600));
			Vector rk_velEarth = rk_method.getVelocity(k1earth, k2earth, k3earth, k4earth, (24 * 3600));
			kVelEarth.operator+=(rk_velEarth);
			Vector K1earth = rk_method.calcK1(kVelEarth);
			Vector K2earth = calcKVelocity(K1earth, kPosEarth, kVelEarth, (12 * 3600));
			Vector K3earth = calcKVelocity(K2earth, kPosEarth, kVelEarth, (12 * 3600));
			Vector K4earth = calcKVelocity(K3earth, kPosEarth, kVelEarth, (24 * 3600));
			//MARS
			Vector k1mars = rk_method.calcK(marsk1Force, Mars.getM());
			Vector k2mars = calcKForce(k1mars, kPosMars, kVelMars, (12 * 3600));
			Vector k3mars = calcKForce(k2mars, kPosMars, kVelMars, (12 * 3600));
			Vector k4mars = calcKForce(k3mars, kPosMars, kVelMars, (24 * 3600));
			Vector rk_velMars = rk_method.getVelocity(k1mars, k2mars, k3mars, k4mars, (24 * 3600));
			kVelMars.operator+=(rk_velMars);
			Vector K1mars = rk_method.calcK1(kVelMars);
			Vector K2mars = calcKVelocity(K1mars, kPosMars, kVelMars, (12 * 3600));
			Vector K3mars = calcKVelocity(K2mars, kPosMars, kVelMars, (12 * 3600));
			Vector K4mars = calcKVelocity(K3mars, kPosMars, kVelMars, (24 * 3600));
			//JUPITER
			Vector k1jupiter = rk_method.calcK(jupiterk1Force, Jupiter.getM());
			Vector k2jupiter = calcKForce(k1jupiter, kPosJupiter, kVelJupiter, (12 * 3600));
			Vector k3jupiter = calcKForce(k2jupiter, kPosJupiter, kVelJupiter, (12 * 3600));
			Vector k4jupiter = calcKForce(k3jupiter, kPosJupiter, kVelJupiter, (24 * 3600));
			Vector rk_velJupiter = rk_method.getVelocity(k1jupiter, k2jupiter, k3jupiter, k4jupiter, (24 * 3600));
			kVelJupiter.operator+=(rk_velJupiter);
			Vector K1jupiter = rk_method.calcK1(kVelJupiter);
			Vector K2jupiter = calcKVelocity(K1jupiter, kPosJupiter, kVelJupiter, (12 * 3600));
			Vector K3jupiter = calcKVelocity(K2jupiter, kPosJupiter, kVelJupiter, (12 * 3600));
			Vector K4jupiter = calcKVelocity(K3jupiter, kPosJupiter, kVelJupiter, (24 * 3600));
			//SATURN
			Vector k1saturn = rk_method.calcK(saturnk1Force, Saturn.getM());
			Vector k2saturn = calcKForce(k1saturn, kPosSaturn, kVelSaturn, (12 * 3600));
			Vector k3saturn = calcKForce(k2saturn, kPosSaturn, kVelSaturn, (12 * 3600));
			Vector k4saturn = calcKForce(k3saturn, kPosSaturn, kVelSaturn, (24 * 3600));
			Vector rk_velSaturn = rk_method.getVelocity(k1saturn, k2saturn, k3saturn, k4saturn, (24 * 3600));
			kVelSaturn.operator+=(rk_velSaturn);
			Vector K1saturn = rk_method.calcK1(kVelSaturn);
			Vector K2saturn = calcKVelocity(K1saturn, kPosSaturn, kVelSaturn, (12 * 3600));
			Vector K3saturn = calcKVelocity(K2saturn, kPosSaturn, kVelSaturn, (12 * 3600));
			Vector K4saturn = calcKVelocity(K3saturn, kPosSaturn, kVelSaturn, (24 * 3600));
			//URANUS
			Vector k1uranus = rk_method.calcK(uranusk1Force, Uranus.getM());
			Vector k2uranus = calcKForce(k1uranus, kPosUranus, kVelUranus, (12 * 3600));
			Vector k3uranus = calcKForce(k2uranus, kPosUranus, kVelUranus, (12 * 3600));
			Vector k4uranus = calcKForce(k3uranus, kPosUranus, kVelUranus, (24 * 3600));
			Vector rk_velUranus = rk_method.getVelocity(k1uranus, k2uranus, k3uranus, k4uranus, (24 * 3600));
			kVelUranus.operator+=(rk_velUranus);
			Vector K1uranus = rk_method.calcK1(kVelUranus);
			Vector K2uranus = calcKVelocity(K1uranus, kPosUranus, kVelUranus, (12 * 3600));
			Vector K3uranus = calcKVelocity(K2uranus, kPosUranus, kVelUranus, (12 * 3600));
			Vector K4uranus = calcKVelocity(K3uranus, kPosUranus, kVelUranus, (24 * 3600));
			//NEPTUNE
			Vector k1neptune = rk_method.calcK(neptunek1Force, Neptune.getM());
			Vector k2neptune = calcKForce(k1neptune, kPosNeptune, kVelNeptune, (12 * 3600));
			Vector k3neptune = calcKForce(k2neptune, kPosNeptune, kVelNeptune, (12 * 3600));
			Vector k4neptune = calcKForce(k3neptune, kPosNeptune, kVelNeptune, (24 * 3600));
			Vector rk_velNeptune = rk_method.getVelocity(k1neptune, k2neptune, k3neptune, k4neptune, (24 * 3600));
			kVelNeptune.operator+=(rk_velNeptune);
			Vector K1neptune = rk_method.calcK1(kVelNeptune);
			Vector K2neptune = calcKVelocity(K1neptune, kPosNeptune, kVelNeptune, (12 * 3600));
			Vector K3neptune = calcKVelocity(K2neptune, kPosNeptune, kVelNeptune, (12 * 3600));
			Vector K4neptune = calcKVelocity(K3neptune, kPosNeptune, kVelNeptune, (24 * 3600));
			//calculating positions and velocities
			nBody_tempVelMercury = rk_method.calcVelocity(nBody_tempVelMercury, k1mercury, k2mercury, k3mercury, k4mercury);
			nBody_tempPosMercury = rk_method.calcPosition(nBody_tempPosMercury, K1mercury, K2mercury, K3mercury, K4mercury);
			
			nBody_tempVelVenus = rk_method.calcVelocity(nBody_tempVelVenus, k1venus, k2venus, k3venus, k4venus);
			nBody_tempPosVenus = rk_method.calcPosition(nBody_tempPosVenus, K1venus, K2venus, K3venus, K4venus);
			
			nBody_tempVelEarth = rk_method.calcVelocity(nBody_tempVelEarth, k1earth, k2earth, k3earth, k4earth);
			nBody_tempPosEarth = rk_method.calcPosition(nBody_tempPosEarth, K1earth, K2earth, K3earth, K4earth);
			
			nBody_tempVelMars = rk_method.calcVelocity(nBody_tempVelMars, k1mars, k2mars, k3mars, k4mars);
			nBody_tempPosMars = rk_method.calcPosition(nBody_tempPosMars, K1mars, K2mars, K3mars, K4mars);
			
			nBody_tempVelJupiter = rk_method.calcVelocity(nBody_tempVelJupiter, k1jupiter, k2jupiter, k3jupiter, k4jupiter);
			nBody_tempPosJupiter = rk_method.calcPosition(nBody_tempPosJupiter, K1jupiter, K2jupiter, K3jupiter, K4jupiter);
		
			nBody_tempVelSaturn = rk_method.calcVelocity(nBody_tempVelSaturn, k1saturn, k2saturn, k3saturn, k4saturn);
			nBody_tempPosSaturn = rk_method.calcPosition(nBody_tempPosSaturn, K1saturn, K2saturn, K3saturn, K4saturn);
			
			nBody_tempVelUranus = rk_method.calcVelocity(nBody_tempVelUranus, k1uranus, k2uranus, k3uranus, k4uranus);
			nBody_tempPosUranus = rk_method.calcPosition(nBody_tempPosUranus, K1uranus, K2uranus, K3uranus, K4uranus);
			
			nBody_tempVelNeptune = rk_method.calcVelocity(nBody_tempVelNeptune, k1neptune, k2neptune, k3neptune, k4neptune);
			nBody_tempPosNeptune = rk_method.calcPosition(nBody_tempPosNeptune, K1neptune, K2neptune, K3neptune, K4neptune);

			//resetting the bodies
			rk_newSun = Body(Sun1.getM(), Sun.getPos(), Sun.getVel());
			rk_newMercury = Body(Mercury.getM(), nBody_tempPosMercury, nBody_tempVelMercury);
			rk_newVenus = Body(Venus.getM(), nBody_tempPosVenus, nBody_tempVelVenus);
			rk_newEarth = Body(Earth.getM(), nBody_tempPosEarth, nBody_tempVelEarth);
			rk_newMars = Body(Mars.getM(), nBody_tempPosMars, nBody_tempVelMars);
			rk_newJupiter = Body(Jupiter.getM(), nBody_tempPosJupiter, nBody_tempVelJupiter);
			rk_newSaturn = Body(Saturn.getM(), nBody_tempPosSaturn, nBody_tempVelSaturn);
			rk_newUranus = Body(Uranus.getM(), nBody_tempPosUranus, nBody_tempVelUranus);
			rk_newNeptune = Body(Neptune.getM(), nBody_tempPosNeptune, nBody_tempVelNeptune);

			k_newMercury = Body(Mercury.getM(), nBody_tempPosMercury, nBody_tempVelMercury);
			k_newVenus = Body(Venus.getM(), nBody_tempPosVenus, nBody_tempVelVenus);
			k_newEarth = Body(Earth.getM(), nBody_tempPosEarth, nBody_tempVelEarth);
			k_newMars = Body(Mars.getM(), nBody_tempPosMars, nBody_tempVelMars);
			k_newJupiter = Body(Jupiter.getM(), nBody_tempPosJupiter, nBody_tempVelJupiter);
			k_newSaturn = Body(Saturn.getM(), nBody_tempPosSaturn, nBody_tempVelSaturn);
			k_newUranus = Body(Uranus.getM(), nBody_tempPosUranus, nBody_tempVelUranus);
			k_newNeptune = Body(Neptune.getM(), nBody_tempPosNeptune, nBody_tempVelNeptune);
			
			kPosMercury = nBody_tempPosMercury;
			kVelMercury = nBody_tempVelMercury;

			kPosVenus = nBody_tempPosVenus;
			kVelVenus = nBody_tempVelVenus;

			kPosEarth = nBody_tempPosEarth;
			kVelEarth = nBody_tempVelEarth;

			kPosMars = nBody_tempPosMars;
			kVelMars = nBody_tempVelMars;

			kPosJupiter = nBody_tempPosJupiter;
			kVelJupiter = nBody_tempVelJupiter;

			kPosSaturn = nBody_tempPosSaturn;
			kVelSaturn = nBody_tempVelSaturn;

			kPosUranus = nBody_tempPosUranus;
			kVelUranus = nBody_tempVelUranus;

			kPosNeptune = nBody_tempPosNeptune;
			kVelNeptune = nBody_tempVelNeptune;

			RkNBody << nBody_tempPosMercury.getX() << "\t" << nBody_tempPosMercury.getY() << "\t" << nBody_tempPosMercury.getZ() << "\t"
				<< nBody_tempPosVenus.getX() << "\t" << nBody_tempPosVenus.getY() << "\t" << nBody_tempPosVenus.getZ() << "\t"
				<< nBody_tempPosEarth.getX() << "\t" << nBody_tempPosEarth.getY() << "\t" << nBody_tempPosEarth.getZ() << "\t"
				<< nBody_tempPosMars.getX() << "\t" << nBody_tempPosMars.getY() << "\t" << nBody_tempPosMars.getZ() << "\t"
				<< nBody_tempPosJupiter.getX() << "\t" << nBody_tempPosJupiter.getY() << "\t" << nBody_tempPosJupiter.getZ() << "\t"
				<< nBody_tempPosSaturn.getX() << "\t" << nBody_tempPosSaturn.getY() << "\t" << nBody_tempPosSaturn.getZ() << "\t"
				<< nBody_tempPosUranus.getX() << "\t" << nBody_tempPosUranus.getY() << "\t" << nBody_tempPosUranus.getZ() << "\t"
				<< nBody_tempPosNeptune.getX() << "\t" << nBody_tempPosNeptune.getY() << "\t" << nBody_tempPosNeptune.getZ() << std::endl;
		}
		RkNBody.close();
	}
	/////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////

	if (simulation == 3)
	{
		//////////////////// Circumbinary orbit ////////////////////
		int bin_initialVelocity;
		int bin_duration;
		int bin_simLength = 0;
		std::cout << "Choose an initial tangential velocity for the planet:" << std::endl;
		std::cout << "\t" << "1 (8,000 m/s)" << std::endl;
		std::cout << "\t" << "2 (12,000 m/s)" << std::endl;
		std::cout << "\t" << "3 (16,000 m/s)" << std::endl;
		std::cin >> bin_initialVelocity;
		std::ofstream LfBinary;
		LfBinary.open("lf_binary_data.txt");
		
		if (bin_initialVelocity == 1)
		{
			Planet = Body((6.42*pow(10, 27)), Vector(0, (-25 * AU), 0), Vector(8000, 0, 0));
		}
		if (bin_initialVelocity == 2)
		{
			Planet = Body((6.42*pow(10, 27)), Vector(0, (-25 * AU), 0), Vector(12000, 0, 0));
		}
		if (bin_initialVelocity == 3)
		{
			Planet = Body((6.42*pow(10, 27)), Vector(0, (-25 * AU), 0), Vector(16000, 0, 0));
		}
		std::cout << "Choose a duration of time to run the simulation for:" << std::endl;
		std::cout << "\t" << "1 (60 years)" << std::endl;
		std::cout << "\t" << "2 (150 years)" << std::endl;
		std::cout << "\t" << "3 (500 years)" << std::endl;
		std::cin >> bin_duration;
		
		if (bin_duration == 1)
		{
			bin_simLength = (60 * 365);
		}
		if (bin_duration == 2)
		{
			bin_simLength = (150 * 365);
		}
		if (bin_duration == 3)
		{
			bin_simLength = (500 * 365);
		}

		std::vector<Body> bin_bodies;
		std::vector<Body>::iterator bin_it;

		bin_it = bin_bodies.begin();
		bin_bodies.push_back(StarA);
		bin_bodies.push_back(StarB);
		bin_bodies.push_back(Planet);

		Vector bin_tempPosStarA = StarA.getPos();
		Vector bin_tempVelStarA = StarA.getVel();
		Vector bin_tempForceStarA = nBodyCalcForce(0, StarA, bin_bodies);
		
		Vector bin_tempPosStarB = StarB.getPos();
		Vector bin_tempVelStarB = StarB.getVel();
		Vector bin_tempForceStarB = nBodyCalcForce(1, StarB, bin_bodies);

		Vector bin_tempPosPlanet = Planet.getPos();
		Vector bin_tempVelPlanet = Planet.getVel();
		Vector bin_tempForcePlanet = nBodyCalcForce(2, Planet, bin_bodies);

		std::cout << "Writing files..." << std::endl;
		for (int o = 0; o <= bin_simLength; o++)
		{
			//setting new positions
			bin_tempPosStarA = lf_method.calcPosition(bin_tempPosStarA,
				bin_tempVelStarA, bin_tempForceStarA, StarA.getM(), (24 * 3600));

			bin_tempPosStarB = lf_method.calcPosition(bin_tempPosStarB,
				bin_tempVelStarB, bin_tempForceStarB, StarB.getM(), (24 * 3600));

			bin_tempPosPlanet = lf_method.calcPosition(bin_tempPosPlanet,
				bin_tempVelPlanet, bin_tempForcePlanet, Planet.getM(), (24 * 3600));

			//setting new bodies
			Body newStarA(StarA.getM(), bin_tempPosStarA, bin_tempVelStarA);
			Body newStarB(StarB.getM(), bin_tempPosStarB, bin_tempVelStarB);
			Body newPlanet(Planet.getM(), bin_tempPosPlanet, bin_tempVelPlanet);

			//replacing bodies in the vector
			bin_bodies.clear();
			bin_bodies.push_back(newStarA);
			bin_bodies.push_back(newStarB);
			bin_bodies.push_back(newPlanet);

			//setting new forces 
			Vector nextSTARAFORCE = nBodyCalcForce(0, newStarA, bin_bodies);
			Vector nextSTARBFORCE = nBodyCalcForce(1, newStarB, bin_bodies);
			Vector nextPLANETFORCE = nBodyCalcForce(2, newPlanet, bin_bodies);

			//setting updated velocities
			bin_tempVelStarA = lf_method.calcVelocity(bin_tempVelStarA,
				bin_tempForceStarA, nextSTARAFORCE, (24 * 3600), StarA.getM());
			bin_tempVelStarB = lf_method.calcVelocity(bin_tempVelStarB,
				bin_tempForceStarB, nextSTARBFORCE, (24 * 3600), StarB.getM());
			bin_tempVelPlanet = lf_method.calcVelocity(bin_tempVelPlanet,
				bin_tempForcePlanet, nextPLANETFORCE, (24 * 3600), Planet.getM());

			//updating the bodies
			newStarA = Body(StarA.getM(), bin_tempPosStarA, bin_tempVelStarA);
			newStarB = Body(StarB.getM(), bin_tempPosStarB, bin_tempVelStarB);
			newPlanet = Body(Planet.getM(), bin_tempPosPlanet, bin_tempVelPlanet);

			//replacing bodies in the vector
			bin_bodies.clear();
			bin_bodies.push_back(newStarA);
			bin_bodies.push_back(newStarB);
			bin_bodies.push_back(newPlanet);

			//updating the force
			bin_tempForceStarA = Vector(nBodyCalcForce(0, newStarA, bin_bodies));
			bin_tempForceStarB = Vector(nBodyCalcForce(1, newStarB, bin_bodies));
			bin_tempForcePlanet = Vector(nBodyCalcForce(2, newPlanet, bin_bodies));


			LfBinary << bin_tempPosStarA.getX() << "\t" << bin_tempPosStarA.getY() << "\t" << bin_tempPosStarA.getZ() << "\t"
				<< bin_tempPosStarB.getX() << "\t" << bin_tempPosStarB.getY() << "\t" << bin_tempPosStarB.getZ() << "\t"
				<< bin_tempPosPlanet.getX() << "\t" << bin_tempPosPlanet.getY() << "\t" << bin_tempPosPlanet.getZ() << std::endl;
		}
		LfBinary.close();
	}
}
