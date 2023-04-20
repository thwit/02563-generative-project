#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#include <random>
#include <glm/glm.hpp>
#include <chrono>
//Libnoise
#include <libnoise/noise.h>

class World{
public:
     void generate();
    void erode(int cycles);
    glm::vec3 surfaceNormal(int i, int j);

    double scale = 60.0;
    double heightmap[256][256] = {0.0};
    glm::vec2 dim = glm::vec2(256);

    int SEED = 0;

    float dt = 1.2;
    float density = 1.0;
    float evapRate = 0.001;
    float depositionRate = 0.1;
    float minVol = 0.01;
    float friction = 0.05;
};

struct Particle{
    Particle(glm::vec2 _pos){ pos = _pos; }

    glm::vec2 pos;
    glm::vec2 speed = glm::vec2(0.0);

    float volume = 1.0;
    float sediment = 0.0;
};

void World::generate(){
    std::cout<<"Generating New World"<<std::endl;
    SEED = time(NULL);
    std::cout<<"Seed: "<<SEED<<std::endl;
    //Seed the Random Generator
    srand(SEED);

    std::cout<<"... generating height ..."<<std::endl;

    //Initialize Heightmap
    noise::module::Perlin perlin;
    perlin.SetOctaveCount(4);
    perlin.SetFrequency(1.0);
    perlin.SetPersistence(0.6);

    float min, max = 0.0;
    for(int i = 0; i < dim.x; i++){
        for(int j = 0; j < dim.y; j++){
            heightmap[i][j] = perlin.GetValue(i*(1.0/dim.x), j*(1.0/dim.y), SEED);
            if(heightmap[i][j] > max) max = heightmap[i][j];
            if(heightmap[i][j] < min) min = heightmap[i][j];
        }
    }

    //Normalize
    for(int i = 0; i < dim.x; i++){
        for(int j = 0; j < dim.y; j++){
            heightmap[i][j] = (heightmap[i][j] - min)/(max - min);
        }
    }
}


glm::vec3 World::surfaceNormal(int i, int j){
  /*
    Note: Surface normal is computed in this way, because the square-grid surface is meshed using triangles.
    To avoid spatial artifacts, you need to weight properly with all neighbors.
  */

  glm::vec3 n = glm::vec3(0.15) * glm::normalize(glm::vec3(scale*(heightmap[i][j]-heightmap[i+1][j]), 1.0, 0.0));  //Positive X
  n += glm::vec3(0.15) * glm::normalize(glm::vec3(scale*(heightmap[i-1][j]-heightmap[i][j]), 1.0, 0.0));  //Negative X
  n += glm::vec3(0.15) * glm::normalize(glm::vec3(0.0, 1.0, scale*(heightmap[i][j]-heightmap[i][j+1])));    //Positive Y
  n += glm::vec3(0.15) * glm::normalize(glm::vec3(0.0, 1.0, scale*(heightmap[i][j-1]-heightmap[i][j])));  //Negative Y

  //Diagonals! (This removes the last spatial artifacts)
  n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale*(heightmap[i][j]-heightmap[i+1][j+1])/sqrt(2), sqrt(2), scale*(heightmap[i][j]-heightmap[i+1][j+1])/sqrt(2)));    //Positive Y
  n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale*(heightmap[i][j]-heightmap[i+1][j-1])/sqrt(2), sqrt(2), scale*(heightmap[i][j]-heightmap[i+1][j-1])/sqrt(2)));    //Positive Y
  n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale*(heightmap[i][j]-heightmap[i-1][j+1])/sqrt(2), sqrt(2), scale*(heightmap[i][j]-heightmap[i-1][j+1])/sqrt(2)));    //Positive Y
  n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale*(heightmap[i][j]-heightmap[i-1][j-1])/sqrt(2), sqrt(2), scale*(heightmap[i][j]-heightmap[i-1][j-1])/sqrt(2)));    //Positive Y

  return n;
}

void World::erode(int cycles){

  /*
    Note: Everything is properly scaled by a time step-size "dt"
  */

  //Do a series of iterations! (5 Particles)
  for(int i = 0; i < cycles; i++){

    //Spawn New Particle
    glm::vec2 newpos = glm::vec2(rand()%(int)dim.x, rand()%(int)dim.y);
    Particle drop(newpos);

    //As long as the droplet exists...
    while(drop.volume > minVol){

      glm::ivec2 ipos = drop.pos;                   //Floored Droplet Initial Position
      glm::vec3 n = surfaceNormal(ipos.x, ipos.y);  //Surface Normal at Position

      //Accelerate particle using newtonian mechanics using the surface normal.
      drop.speed += dt*glm::vec2(n.x, n.z)/(drop.volume*density);//F = ma, so a = F/m
      drop.pos   += dt*drop.speed;
      drop.speed *= (1.0-dt*friction);       //Friction Factor

      /*
        Note: For multiplied factors (e.g. friction, evaporation)
        time-scaling is correctly implemented like above.
      */

      //Check if Particle is still in-bounds
      if(!glm::all(glm::greaterThanEqual(drop.pos, glm::vec2(0))) ||
         !glm::all(glm::lessThan(drop.pos, dim))) break;

      //Compute sediment capacity difference
      float maxsediment = drop.volume*glm::length(drop.speed)*(heightmap[ipos.x][ipos.y]-heightmap[(int)drop.pos.x][(int)drop.pos.y]);
      if(maxsediment < 0.0) maxsediment = 0.0;
      float sdiff = maxsediment - drop.sediment;

      //Act on the Heightmap and Droplet!
      drop.sediment += dt*depositionRate*sdiff;
      heightmap[ipos.x][ipos.y] -= dt*drop.volume*depositionRate*sdiff;

      //Evaporate the Droplet (Note: Proportional to Volume! Better: Use shape factor to make proportional to the area instead.)
      drop.volume *= (1.0-dt*evapRate);
    }
  }
}


void saveHeightmap(double hmap[][256], const char *fname) {
    std::fstream of(fname, std::ios::out | std::ios::trunc);
    for (int x = 0; x < 256; ++x) {
        for (int y = 0; y < 256; ++y) {
            of << hmap[x][y] << " ";
        }
        of << "\n";
    }
    of.close();
}

int main( int argc, char* args[] ) {
	World world;
	world.generate();
    saveHeightmap(world.heightmap, "originalmap.txt");
    
    int stepsLeft = 250000;
    int stepSize = 10000;
	while(stepsLeft){
        std::cout << stepsLeft << std::endl;
        world.erode(stepSize);
        stepsLeft -= stepSize;
	}
    
    std::fstream of("Map.txt", std::ios::out | std::ios::trunc);

    saveHeightmap(world.heightmap, "map.txt");

	return 0;
}