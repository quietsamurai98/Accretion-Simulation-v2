/********************
*
* NBODY GRAVITATION
*
*********************/

#include <math.h>
#include <iostream>
#include <string>

#define uint unsigned int
#define BODY_COUNT	3
#define DELTA_TIME	0.01
#define DT_SQ_HALF	(DELTA_TIME * DELTA_TIME * 0.5)
#define DT_HALF		(DELTA_TIME * 0.5)
#define GRAV_CONST	1

struct Vector {
	double x;
	double y;
};

struct Mass {
 private:
	double value;
	double inverse;
 public:
	void set(double new_mass){
		value = new_mass;
		inverse = 1/value;
	}
	double get(){
		return value;
	}
	double inv(){
		return inverse;
	}
};
struct Body {
	bool   alive;
	Mass   mass;
	Vector pos;
	Vector vel;
	Vector acc;
	Vector new_pos;
	Vector new_vel;
	Vector new_acc;
	Vector force;
	
	void calc_pos(){
		new_pos = pos;
		new_pos.x += vel.x*DELTA_TIME + acc.x*DT_SQ_HALF;
		new_pos.y += vel.y*DELTA_TIME + acc.y*DT_SQ_HALF;
	}
	void calc_acc(){
		new_acc = force;
		new_acc.x*=mass.inv();
		new_acc.y*=mass.inv();
	}
	void calc_vel(){
		new_vel = vel;
		new_vel.x += (acc.x + new_acc.x)*DT_HALF;
		new_vel.y += (acc.y + new_acc.y)*DT_HALF;
	}
	void calc_force(Body (&universe)[BODY_COUNT], uint idx){
		for(uint i = idx+1; i < BODY_COUNT; i++){
			Body &other = universe[i];
			if(!other.alive)
				continue; //Skip particle if it isn't "alive"
			Vector disp = {new_pos.x - other.new_pos.x, new_pos.y - other.new_pos.y}; //Displacement vector
			double dist = sqrt(disp.x * disp.x + disp.y * disp.y); //Displacement scalar (distance between bodies)
			double scalar_force = -GRAV_CONST * mass.get() * other.mass.get() / (dist*dist*dist); //Muliply dist by scalar_force to get force vector
			
			//disp is now used as the force vector, despite the name.
			disp.x*=scalar_force;
			disp.y*=scalar_force;
			
			force.x += disp.x;
			force.y += disp.y;
			other.force.x -= disp.x;
			other.force.y -= disp.y;
		}
	}
	void update(){
		pos = new_pos;
		acc = new_acc;
		vel = new_vel;
		force = {0,0};
	}
	
	std::string to_string(){
		std::string out = "mass=";
		out+=std::to_string(mass.get());
		out+=" pos=(";
		out+=std::to_string(pos.x);
		out+=",";
		out+=std::to_string(pos.y);
		out+=") vel=(";
		out+=std::to_string(vel.x);
		out+=",";
		out+=std::to_string(vel.y);
		out+=") acc=(";
		out+=std::to_string(acc.x);
		out+=",";
		out+=std::to_string(acc.y);
		out+=")";
		return out;
	}
};

int main(int argc, char *argv[]) {
	Body universe[BODY_COUNT] = { 0 };
	universe[0].mass.set(1);
	universe[1].mass.set(1);
	universe[2].mass.set(0.01);
	universe[0].pos = { 1,0};
	universe[1].pos = {-1,0};
	universe[2].pos = {10,0};
	universe[0].vel = {0,0.4};
	universe[1].vel = {0,-.4};
	universe[2].vel = {0,std::stod(argv[1])};
	
	Body barycenter = { 0 };
	
	for(uint i = 0; i < BODY_COUNT; i++){
		universe[i].alive = true;
		std::cout << "Body " << i << " X Position";
		std::cout << ',';
		std::cout << "Body " << i << " Y Position";
		/*
		std::cout << ',';
		std::cout << "Body " << i << " X Velocity";
		std::cout << ',';
		std::cout << "Body " << i << " Y Velocity";
		std::cout << ',';
		std::cout << "Body " << i << " X Acceleration";
		std::cout << ',';
		std::cout << "Body " << i << " Y Acceleration";
		*/
		std::cout << ',';
		std::cout << ',';
	}
	std::cout << "Barycenter X Position (absolute)";
	std::cout << ',';
	std::cout << "Barycenter Y Position (absolute)";
	std::cout << ',';
	std::cout << "Barycenter X Velocity (absolute)";
	std::cout << ',';
	std::cout << "Barycenter Y Velocity (absolute)";
	std::cout << ',';
	std::cout << "Barycenter X Acceleration (absolute)";
	std::cout << ',';
	std::cout << "Barycenter Y Acceleration (absolute)";
	std::cout << std::endl;
	
	for(uint tick = 0; tick < (unsigned)std::stoi(argv[2]); tick++){
		
		barycenter = { 0 };
		for(uint i = 0; i < BODY_COUNT; i++){
			Body &b = universe[i];
			double m = b.mass.get();
			barycenter.mass.set(m+barycenter.mass.get());
			barycenter.pos.x+=m*b.pos.x;
			barycenter.pos.y+=m*b.pos.y;
			barycenter.vel.x+=m*b.vel.x;
			barycenter.vel.y+=m*b.vel.y;
			barycenter.acc.x+=m*b.acc.x;
			barycenter.acc.y+=m*b.acc.y;
		}
		barycenter.pos.x*=barycenter.mass.inv();
		barycenter.pos.y*=barycenter.mass.inv();
		barycenter.vel.x*=barycenter.mass.inv();
		barycenter.vel.y*=barycenter.mass.inv();
		barycenter.acc.x*=barycenter.mass.inv();
		barycenter.acc.y*=barycenter.mass.inv();
		
		for(uint i = 0; i < BODY_COUNT; i++){
			Body &b = universe[i];
			std::cout << b.pos.x - barycenter.pos.x;
			std::cout << ',';
			std::cout << b.pos.y - barycenter.pos.y;
			/*
			std::cout << ',';
			std::cout << b.vel.x - barycenter.vel.x;
			std::cout << ',';
			std::cout << b.vel.y - barycenter.vel.y;
			std::cout << ',';
			std::cout << b.acc.x - barycenter.acc.x;
			std::cout << ',';
			std::cout << b.acc.y - barycenter.acc.y;
			*/
			std::cout << ',';
			std::cout << ',';
		}
		std::cout << barycenter.pos.x << ',' << barycenter.pos.y << ',' << barycenter.vel.x << ',' << barycenter.vel.y << ',' << barycenter.acc.x << ',' << barycenter.acc.y;
		std::cout << std::endl;
		for(uint i = 0; i < BODY_COUNT; i++)
			universe[i].calc_pos();
		for(uint i = 0; i < BODY_COUNT; i++)
			universe[i].calc_force(universe, i);
		for(uint i = 0; i < BODY_COUNT; i++){
			universe[i].calc_acc();
			universe[i].calc_vel();
			universe[i].update();
		}
	}
}