/*
g++ render.cpp -o render  -O3 -Wall -std=c++17 -L/usr/X11R6/lib -lm -lpthread -lX11
./render <HISTORY FILE> 1 50000 | ffmpeg -framerate 60 -r 60 -y -f rawvideo -pixel_format gbrp -video_size 1920x1080 -i - <OUTPUT VIDEO FILE>
*/

#include <math.h>
#include <iostream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "CImg.h"

#define uint uint64_t

namespace CImg = cimg_library;

static const unsigned char WHITE[] = {255, 255, 255};
//static const unsigned char WHITE[] = {255};
static const unsigned char BLACK[] = {  0,   0,   0};
//static const unsigned char BLACK[] = {  0};

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

#pragma pack(push, 1)
struct Header{
	char blurb1[16];
	uint body_count;
	uint tick_count;
	char blurb2[16];
};
#pragma pack(pop)

#pragma pack(push, 1)
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
#pragma pack(pop)

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

#pragma pack(push, 1)
struct Body{
	double mass;
	double radius;
	Vector pos;
	Vector vel;
	Vector acc;
};
#pragma pack(pop)


struct View {
	int h;
	int w;
	double x_min;
	double x_max;
	double y_min;
	double y_max;
 private:
	size_t px_count;
	double dX;
	double dXi;
	double dY;
	double dYi;
 public:
	size_t pixels(){
		if(!px_count)
			return px_count;
		return px_count = w*h;
	}
	double deltaX(){
		if(!dX)
			return dX;
		return dX = (x_max-x_min)/((double)w);
	}
	double invDeltaX(){
		if(!dXi)
			return dXi;
		return dXi = ((double)w)/(x_max-x_min);
	}
	double deltaY(){
		if(!dY)
			return dY;
		return dY = (y_max-y_min)/((double)h);
	}
	double invDeltaY(){
		if(!dYi)
			return dYi;
		return dYi = ((double)h)/(y_max-y_min);
	}
	int x2c(double x){
		return (int)((x-x_min)*invDeltaX());
	}
	int y2r(double y){
		return (int)((y_max-y)*invDeltaY());
	}
};

void draw_body (CImg::CImg<unsigned char> &img, View &view, Body &body){
	if(!body.radius)
		return;
	int col = view.x2c(body.pos.x);
	int row = view.y2r(body.pos.y);
	int rad = (int)(body.radius * view.invDeltaX());
	if(!rad){
		img.draw_point(col, row, WHITE);
	} else {
		img.draw_circle(col, row, rad, WHITE);
	}
}

/***
*
* Read the binary header.
*
* Returns 0 on success, 1 on failure
*
***/
int read_header(Header &head, FILE *bin){
	if(1 != fread(&head, sizeof(Header), 1, bin)){
		//Could not read header at all.
		return 1;
	}
	if(memcmp(head.blurb1, "NBODY SIMULATION", 16) || memcmp(head.blurb2, "UNIVERSE HISTORY", 16)){
		//Blurbs do not match expected values, indicating a malformed simulation history binary.
		return 1;
	}
	return 0;
}

/***
*
* Read the next frame of the simulation into the universe array and barycenter body.
*
* Returns 0 on success, 1 on proper EOF, and 2 in the case of an error.
*	(Proper EOF := EOF occurs at the end of a simulation frame)
*
***/
int next_frame(Body *universe, Body &barycenter, Header &head, FILE *bin){
	size_t read_count = fread(universe, sizeof(Body), head.body_count, bin);
	if(head.body_count!=read_count){
		if(read_count == 0){
			return 1; //Proper EOF
		} else {
			return 2; //Couldn't read, or an EOF occurs mid-frame
		}
	}
	if(1!=fread(&barycenter, sizeof(Body), 1, bin)){
		return 2; //Couldn't read, or an EOF occurs mid-frame
	}
	return 0;
}

void output_frame(CImg::CImg<unsigned char> &image, View &view){
	char* s=reinterpret_cast<char*>(image.data()+view.pixels());	// Get start of G plane
	std::cout.write(s,view.pixels());                               // Output it
	s=reinterpret_cast<char*>(image.data()+2*view.pixels());        // Get start of B plane
	std::cout.write(s,view.pixels());                               // Output it
	s=reinterpret_cast<char*>(image.data());                        // Get start of R plane
	std::cout.write(s,view.pixels());                               // Output it
}

void process_frame(Body *universe, Body &barycenter, uint current_tick, CImg::CImg<unsigned char> &image, Header &head, View &view){
	image.fill(0);
	for(uint i = 0; i < head.body_count; ++i){
		//std::cout << "TICK " << current_tick << " BODY " << i << " POS: (" << universe[i].pos.to_string() << ')' << " RADIUS: " << universe[i].radius << std::endl;
		//std::cout << "TICK " << current_tick << " BODY " << i << " VEL: (" << universe[i].vel.to_string() << ')' << std::endl;
		//std::cout << "TICK " << current_tick << " BODY " << i << " ACC: (" << universe[i].acc.to_string() << ')' << std::endl;
		draw_body(image, view, universe[i]);
	}
	//char buffer[40];
	//sprintf(buffer, "./frames/%020lu.bmp", current_tick);
	//printf("%s\n", buffer);
	//image.save(buffer);
	output_frame(image, view);
	
	//std::cout << "TICK " << current_tick << " BARY POS: (" << barycenter.pos.to_string() << ')' << std::endl;
	//std::cout << "TICK " << current_tick << " BARY VEL: (" << barycenter.vel.to_string() << ')' << std::endl;
	//std::cout << "TICK " << current_tick << " BARY ACC: (" << barycenter.acc.to_string() << ')' << std::endl;
}

int main(int argc, char *argv[]) {
	FILE *bin = fopen(argv[1], "rb"); //Binary input file;
	
	View view;
	view.h = 1080;
	view.w = 1920;
	view.x_min = -176.0/9.0;
	view.x_max =  176.0/9.0;
	view.y_min = -11;
	view.y_max =  11;
	
	CImg::CImg<unsigned char> image(view.w,view.h,1,3);
	image.fill(0);
	
	Header head;
	if(read_header(head, bin)){
		std::cerr << "Could not read header!" << std::endl;
		return EXIT_FAILURE;
	}
	
	Body *universe = (Body*) malloc(head.body_count * sizeof(Body));
	Body barycenter;
	int  status; //Tracks the status of the simulation readback. 0=good to go, 1=expected EOF, 2=error
	uint current_tick = 0;
	
	uint decimation_rate = argc > 2 ? std::stoull(argv[2]) : 1;
	uint max_tick		 = argc > 3 ? std::stoull(argv[3]) : 1;
	
	/*Read in and process all the frames sequentially*/
	while(!(status=next_frame(universe, barycenter, head, bin))){
		if(decimation_rate==1 || current_tick%decimation_rate == 0)
			process_frame(universe, barycenter, current_tick, image, head, view);
		current_tick++;
		if(current_tick == max_tick)
			break;
	}
	
	/*Check if read/process loop ended due to an error when reading*/
	if(status == 2){
		std::cerr << "Malformed frame!" << std::endl;
		return EXIT_FAILURE;
	}
	
	image.save("./orbit.bmp");
	
	return EXIT_SUCCESS;
}