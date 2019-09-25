/********************
*
* NBODY GRAVITATION
*
*********************/

/*
g++ mainV3.cpp -o nbodyV3 -O2 -Wall -std=c++17 -pthread -funroll-loops
./nbodyV3 <FILE TO SAVE HISTORY IN> 1.0 0.5 36120
*/

#include <math.h>
#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <cstring>
#include <stdio.h>
#include <stdint.h>
#include "ThreadPool.h"

#define uint uint64_t
#define BODY_COUNT	 1000
#define DELTA_TIME	 0.01
#define DT_SQ_HALF	 (DELTA_TIME * DELTA_TIME * 0.5)
#define DT_HALF		 (DELTA_TIME * 0.5)
#define GRAV_CONST	 1
#define DIMENSIONS	 2
#define SERIAL_BODY_SIZE (((__SIZEOF_DOUBLE__ * DIMENSIONS) * 3) + (2 * __SIZEOF_DOUBLE__))

#ifndef PI
#define PI (3.14159265358979323846)
#endif

std::string dtos(double x){
	char *buf;
	
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wunused-result"
	asprintf(&buf, "%.*e", 16, x); //Allocate suitably large cstring, populate with text, ignore the returned length.
	#pragma GCC diagnostic pop
	
	std::string out = buf;
	free(buf);
	return out;
}

struct UnivIdx {
	size_t arr[BODY_COUNT];
	size_t len;
};

struct Vector {
	double x;
	double y;
	Vector operator+= (const Vector &obj){
		x+=obj.x;
		y+=obj.y;
		return *this;
	}
	Vector operator-= (const Vector &obj){
		x-=obj.x;
		y-=obj.y;
		return *this;
	}
	Vector operator*= (double scalar){
		x*=scalar;
		y*=scalar;
		return *this;
	}
	Vector operator/= (double scalar){
		x/=scalar;
		y/=scalar;
		return *this;
	}
	std::string to_string(){
		return dtos(x) + ',' + dtos(y);
	}
};

Vector operator+ (const Vector &lhs, const Vector &rhs){
	Vector out = {
		lhs.x+rhs.x,
		lhs.y+rhs.y
	};
	return out;
}
Vector operator- (const Vector &lhs, const Vector &rhs){
	Vector out = {
		lhs.x-rhs.x,
		lhs.y-rhs.y
	};
	return out;
}
Vector operator* (const Vector &lhs, double rhs){
	Vector out = {
		lhs.x*rhs,
		lhs.y*rhs
	};
	return out;
}
Vector operator/ (const Vector &lhs, double rhs){
	Vector out = {
		lhs.x/rhs,
		lhs.y/rhs
	};
	return out;
}
Vector operator* (double lhs, const Vector &rhs){
	return rhs*lhs;
}
Vector operator/ (double lhs, const Vector &rhs){
	return rhs/lhs;
}

struct Mass {
 private:
	double value;
	double inverse;
	double radius;
 public:
	void set(double new_mass){
		value = new_mass;
		radius = sqrt(new_mass)*.25;
		if(value == 0){
			inverse = 0;
		} else {
			inverse = 1/value;
		}
	}
	double get(){
		return value;
	}
	double inv(){
		return inverse;
	}
	double rad(){
		return radius;
	}
};
struct UnivIdx;
struct Body {
	bool   alive;
	bool   collide;
	Mass   mass;
	Vector pos;
	Vector vel;
	Vector acc;
	Vector new_pos;
	Vector new_vel;
	Vector new_acc;
	Vector force;
	
	void calc_pos(){
		new_pos = pos + (vel*DELTA_TIME) + (acc*DT_SQ_HALF);
	}
	void calc_force(Body (&universe)[BODY_COUNT], uint idx, UnivIdx &univ_idx){
		for(uint i = 0; i < univ_idx.len; ++i){
			if(i == idx)
				continue;
			Body &other = universe[univ_idx.arr[i]];
			if(!other.alive)
				continue; //Skip particle if it isn't "alive"
			
			Vector disp = new_pos - other.new_pos; //Displacement vector
			double dist_sq = disp.x * disp.x + disp.y * disp.y;
			double dist = sqrt(dist_sq); //Displacement scalar (distance between bodies)
			double padded_divisor = (dist_sq*dist) + 0.001; //Add an epsilon to prevent bodies that get too close from flinging eachother away at ludicrous speed
			double scalar_force = -GRAV_CONST * mass.get() * other.mass.get() / padded_divisor; //Note: if you muliply dist by scalar_force, you get the force vector
			
			//disp is now used as the force vector, despite the name.
			disp*=scalar_force;
			
			force += disp;
			
			if((other.mass.rad()+mass.rad()) > dist){
				collide = true;
			}
		}
	}
	void calc_acc(){
		new_acc = force * mass.inv();
	}
	void calc_vel(){
		new_vel = vel + (acc + new_acc)*DT_HALF;
	}
	void update(){
		pos = new_pos;
		acc = new_acc;
		vel = new_vel;
		force = {0,0};
	}
	
	void serialize(char (&dest_buf)[SERIAL_BODY_SIZE]){
		unsigned short s = sizeof(double);
		double tmp_m = mass.get();
		double tmp_r = mass.rad();
		std::memcpy(&dest_buf[s*0], &tmp_m, s);	//MASS
		std::memcpy(&dest_buf[s*1], &tmp_r, s);	//RADIUS
		std::memcpy(&dest_buf[s*2], &pos.x, s);	//POS_X
		std::memcpy(&dest_buf[s*3], &pos.y, s);	//POS_Y
		std::memcpy(&dest_buf[s*4], &vel.x, s);	//VEL_X
		std::memcpy(&dest_buf[s*5], &vel.y, s);	//VEL_Y
		std::memcpy(&dest_buf[s*6], &acc.x, s);	//ACC_X
		std::memcpy(&dest_buf[s*7], &acc.y, s);	//ACC_Y
	}
};

void update_barycenter(Body &barycenter, Body (&universe)[BODY_COUNT], UnivIdx &univ_idx){
	
	if(barycenter.mass.get() == 0){
		double m = 0;
		for(uint i = 0; i < univ_idx.len; ++i){
			Body &b = universe[univ_idx.arr[i]];
			m += b.mass.get();
		}
		barycenter.mass.set(m);
	}
	
	barycenter.pos = { 0 };
	barycenter.vel = { 0 };
	barycenter.acc = { 0 };
	for(uint i = 0; i < univ_idx.len; ++i){
		Body &b = universe[univ_idx.arr[i]];
		double m = b.mass.get();
		barycenter.pos+=m*b.pos;
		barycenter.vel+=m*b.vel;
		barycenter.acc+=m*b.acc;
	}
	barycenter.pos*=barycenter.mass.inv();
	barycenter.vel*=barycenter.mass.inv();
	barycenter.acc*=barycenter.mass.inv();
}

void write_csv_header(){
	for(uint i = 0; i < BODY_COUNT; ++i){
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
}

void write_csv_frame(Body &barycenter, Body (&universe)[BODY_COUNT]){
	for(uint i = 0; i < BODY_COUNT; ++i){
		Body &b = universe[i];
		if(b.alive)
			std::cout << (b.pos - barycenter.pos).to_string();
		else
			std::cout << ",";
		/*
		std::cout << ',';
		std::cout << (b.vel - barycenter.vel).to_string();
		std::cout << ',';
		std::cout << (b.acc - barycenter.acc).to_string();
		*/
		std::cout << ',';
		std::cout << ',';
	}
	std::cout << barycenter.pos.to_string() << ',' << barycenter.vel.to_string() << ',' << barycenter.acc.to_string() << std::endl;
}

void write_bin_header(uint tick_limit, FILE *bout){
	unsigned short s = sizeof(uint);
	uint count = BODY_COUNT;
	uint ticks = tick_limit;
	char blurb1[]="NBODY SIMULATION";
	char blurb2[]="UNIVERSE HISTORY";
	char dest_buf[32+sizeof(uint)*2];
	
	std::memcpy(&dest_buf[00+s*0], &blurb1,16);
	std::memcpy(&dest_buf[16+s*0], &count, s);
	std::memcpy(&dest_buf[16+s*1], &ticks, s);
	std::memcpy(&dest_buf[16+s*2], &blurb2,16);
	fwrite(dest_buf, sizeof(char), 32 + sizeof(uint)*2, bout);
	fflush(bout);
}

void write_bin_frame(Body &barycenter, Body (&universe)[BODY_COUNT], FILE *bout){
	char buf[SERIAL_BODY_SIZE];
	for(uint i = 0; i < BODY_COUNT; ++i){
		universe[i].serialize(buf);
		fwrite(buf, sizeof(char), SERIAL_BODY_SIZE, bout);
	}
	barycenter.serialize(buf);
	fwrite(buf, sizeof(char), SERIAL_BODY_SIZE, bout);
	fflush(bout);
}

void create_universe(Body (&universe)[BODY_COUNT], Body &barycenter, int argc, char *argv[]){
	double DISK_RADIUS = 10.0;
	double INIT_MASS   = 0.001;
	double VEL_MEAN     = std::stod(argv[2]);
	double VEL_STDDEV   = std::stod(argv[3]);
	
	std::uniform_real_distribution<double> rand_u(0.0,1.0);
	std::normal_distribution<double> rand_n(VEL_MEAN,VEL_STDDEV);
	std::default_random_engine rand_engn;
	rand_engn.seed(std::chrono::system_clock::now().time_since_epoch().count());
	auto rand_unif = [&rand_u, &rand_engn](){return rand_u(rand_engn);};
	auto rand_nrml = [&rand_n, &rand_engn](){return rand_n(rand_engn);};
	
	/*
	universe[0].mass.set(1);
	universe[1].mass.set(1);
	universe[2].mass.set(0.01);
	universe[0].pos = { 1,0};
	universe[1].pos = {-1,0};
	universe[2].pos = {10,0};
	universe[0].vel = {0,0.4};
	universe[1].vel = {0,-.4};
	universe[2].vel = {0,std::stod(argv[2])}; //.45 for stable orbit
	*/
	
	universe[0].mass.set(4);
	universe[0].alive = true;
	for(uint i = 1; i < BODY_COUNT; i++){
		Body &b = universe[i];
		b.alive = true;
		b.mass.set(INIT_MASS);
		double radial_dist = pow(rand_unif(), 0.75) * DISK_RADIUS; 
		double theta = rand_unif()*PI*2.0;
		b.pos = { cos(theta), sin(theta) };
		b.pos*= radial_dist;
		
		//theta = rand_unif()*PI*2.0;
		theta = atan2(b.pos.y,b.pos.x) + 0.5*PI;
		b.vel = { cos(theta), sin(theta) };
		b.vel*= rand_nrml() * tanh(radial_dist*PI/DISK_RADIUS); //Slow in middle, faster near edge
	}
}

void rebuild_idx(Body (&universe)[BODY_COUNT], UnivIdx &univ_idx){
	univ_idx.len = 0;
	for(size_t i = univ_idx.arr[0]; i < BODY_COUNT; ++i){
		if(universe[i].alive){
			univ_idx.arr[univ_idx.len++]=i;
		}
	}
}

void collide_universe(Body (&universe)[BODY_COUNT], UnivIdx &univ_idx){
	
	//Get indicies of live bodies with collision flag set
	size_t idx_arr[BODY_COUNT];
	size_t idx_len = 0;
	for(uint i = 0; i < BODY_COUNT; ++i){
		if(universe[i].alive && universe[i].collide){
			idx_arr[idx_len]=i;
			idx_len++;
		}
	}
	
	bool idx_dirty = false;
		
	for(uint i = 0; i < idx_len; ++i){
		Body &a = universe[idx_arr[i]];
		if(!(a.alive && a.collide))
			continue; //Skip dead and non-colliding particles
		for(uint j = i+1; j < idx_len; ++j){
			Body &b = universe[idx_arr[j]];
			if(!(b.alive && b.collide))
				continue; //Skip dead and non-colliding particles
			
			Vector disp = a.pos - b.pos;
			double dist_sq = disp.x * disp.x + disp.y * disp.y;
			double m_ab	   = a.mass.get()+b.mass.get();
			double r_ab	   = a.mass.rad()+b.mass.rad();
			if(r_ab*r_ab > dist_sq){
				//a and b are colliding!
				double a_m = a.mass.get();
				double b_m = b.mass.get();
				
				a.mass.set(m_ab);
				a.pos = ((a.pos*a_m)+(b.pos*b_m))/(m_ab);
				a.vel = ((a.vel*a_m)+(b.vel*b_m))/(m_ab);
				a.acc = ((a.acc*a_m)+(b.acc*b_m))/(m_ab); //AFAIK, averaging the accelerations between two colliding bodies makes little sense, but ¯\_(ツ)_/¯
				
				b = { 0 };
				b.alive=false;
				b.collide = false;
				idx_dirty = true;
			}
		}
		a.collide = false;
	}
	if(idx_dirty)
		rebuild_idx(universe, univ_idx);
}

int main(int argc, char *argv[]) {
	FILE *bout = fopen(argv[1], "wb"); //Binary output file
	char *bbuf = (char*) malloc((BODY_COUNT+1)*SERIAL_BODY_SIZE);
	setbuf(bout, bbuf);
	
	bool PRINT_CSV = argc > 5;
	
	progschj::ThreadPool pool;
	
	uint tick_limit = (unsigned)std::stoull(argv[4]);
	
	if(PRINT_CSV)
		write_csv_header();
	
	write_bin_header(tick_limit, bout);
	
	Body universe[BODY_COUNT] = { 0 };
	Body barycenter = { 0 };
	
	if(!PRINT_CSV)
		printf("Creating universe...\r\n");
	create_universe(universe, barycenter, argc, argv);
	if(!PRINT_CSV)
		printf("Universe created!\r\n");
	
	UnivIdx univ_idx = { 0 };
	rebuild_idx(universe, univ_idx);
	
	update_barycenter(barycenter, universe, univ_idx);
	//Ensure the universe is using barycentric coordinates and reference frame
	for(uint i = 0; i < BODY_COUNT; i++){
		universe[i].pos-=barycenter.pos;
		universe[i].vel-=barycenter.vel;
	}
	update_barycenter(barycenter, universe, univ_idx);
	
	int pad_len = (int)(0.5+log10(tick_limit))+1;
	
	int csv_skip_factor = 1;
	if(tick_limit > 25000)
		csv_skip_factor = (tick_limit/25000)+1;

	for(uint tick = 0; tick < tick_limit; ++tick){
		collide_universe(universe, univ_idx);
		
		update_barycenter(barycenter, universe, univ_idx);
		write_bin_frame(barycenter, universe, bout);
		
		if(PRINT_CSV && !(tick%csv_skip_factor)){
			write_csv_frame(barycenter, universe);			
		}
		
		for(uint i = 0; i < univ_idx.len; ++i){
			pool.enqueue([i, &universe, &univ_idx]{
				universe[univ_idx.arr[i]].calc_pos();
			});
		}
		pool.wait_until_empty();
		pool.wait_until_nothing_in_flight();
		
		for(uint i = 0; i < univ_idx.len; ++i){
			pool.enqueue([i, &universe, &univ_idx]{
				universe[univ_idx.arr[i]].calc_force(universe, i, univ_idx);
			});
		}
		pool.wait_until_empty();
		pool.wait_until_nothing_in_flight();
		
		for(uint i = 0; i < univ_idx.len; ++i){
			pool.enqueue([i, &universe, &univ_idx]{
				universe[univ_idx.arr[i]].calc_acc();
				universe[univ_idx.arr[i]].calc_vel();
				universe[univ_idx.arr[i]].update();
			});
		}
		pool.wait_until_empty();
		pool.wait_until_nothing_in_flight();
		
		if(!PRINT_CSV){// && tick%10==0){
			printf("%0*lu/%lu\r",pad_len,tick,tick_limit);
			std::cout << std::flush;
		}
	}
	
	fflush(bout);
	fclose(bout);
	free(bbuf);
}