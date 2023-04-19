#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#include <random>
#include "eigen3/Eigen/Dense"
#include <chrono>

#define DB_PERLIN_IMPL
#include "db_perlin.hpp"

struct Particle {
    Eigen::Vector2d pos = Eigen::Vector2d::Zero();
    Eigen::Vector2d speed = Eigen::Vector2d::Zero();
    double volume = 1.;
    double sediment = 0.;
};

class World {
public:
    double heightmap[256][256] = {0.};
    double dt = 1.2;
    double min_volume = 0.01;
    double density = 1.;
    double friction = 0.05;
    double deposition_rate = 0.1;
    double evaporation_rate = 0.001;
    double scale = 60.;
    void generate();
    Eigen::Vector3d surface_normal(int i, int j);
    void erode(int cycles);
};

void World::generate() {
    float min, max = 0.;
	for (int x = 0; x < 256; ++x)
	{
		for (int y = 0; y < 256; ++y)
		{
			const double noise = db::perlin((x * 0.0075), (y * 0.0075));
			
			heightmap[x][y] = noise;
            if(heightmap[x][y] > max) max = heightmap[x][y];
            if(heightmap[x][y] < min) min = heightmap[x][y];
		}
	}

      //Normalize
    for(int x = 0; x < 256; x++){
        for(int y = 0; y < 256; y++){
            heightmap[x][y] = (heightmap[x][y] - min)/(max - min);
    }
  }
}

Eigen::Vector3d World::surface_normal(int i, int j) {
    // Compute surface normal for x neighbors
    Eigen::Vector3d nx_pos(scale * (heightmap[i][j] - heightmap[i+1][j]), 1.0, 0.0);
    Eigen::Vector3d nx_neg(scale * (heightmap[i-1][j] - heightmap[i][j]), 1.0, 0.0);
    // Compute surface normal for y neighbors
    Eigen::Vector3d ny_pos(0.0, 1.0, scale * (heightmap[i][j] - heightmap[i][j+1]));
    Eigen::Vector3d ny_neg(0.0, 1.0, scale * (heightmap[i][j-1] - heightmap[i][j]));
    // Compute surface normal for diagonal neighbors
    Eigen::Vector3d n_diag1(scale * (heightmap[i][j] - heightmap[i+1][j+1]) / std::sqrt(2), std::sqrt(2), scale * (heightmap[i][j] - heightmap[i+1][j+1]) / std::sqrt(2));
    Eigen::Vector3d n_diag2(scale * (heightmap[i][j] - heightmap[i+1][j-1]) / std::sqrt(2), std::sqrt(2), scale * (heightmap[i][j] - heightmap[i+1][j-1]) / std::sqrt(2));
    Eigen::Vector3d n_diag3(scale * (heightmap[i][j] - heightmap[i-1][j+1]) / std::sqrt(2), std::sqrt(2), scale * (heightmap[i][j] - heightmap[i-1][j+1]) / std::sqrt(2));
    Eigen::Vector3d n_diag4(scale * (heightmap[i][j] - heightmap[i-1][j-1]) / std::sqrt(2), std::sqrt(2), scale * (heightmap[i][j] - heightmap[i-1][j-1]) / std::sqrt(2));
    // Compute total surface normal vector by summing all surface normal vectors with appropriate weights
    Eigen::Vector3d n = 0.15 * (nx_pos.norm() * nx_pos + nx_neg.norm() * nx_neg + ny_pos.norm() * ny_pos + ny_neg.norm() * ny_neg);
    n += 0.1 * (n_diag1.norm() * n_diag1 + n_diag2.norm() * n_diag2 + n_diag3.norm() * n_diag3 + n_diag4.norm() * n_diag4);
    // Normalize the total surface normal vector and return
    return n.normalized();
}

void World::erode(int cycles) {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::microseconds;

    World world = World();
    
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(2, 254); // define the range

    int prev_finished = -1;
    int prev_i = -1;

            
    auto t1 = high_resolution_clock::now();
    auto now_ms = std::chrono::time_point_cast<std::chrono::microseconds>(t1);
    auto epoch = now_ms.time_since_epoch();
    auto value = std::chrono::duration_cast<std::chrono::microseconds>(epoch);
    long t1_ms = value.count();

    for (int i = 0; i < cycles; i++) {        
        int finished = 100 * (i / (double) cycles);
        if (prev_finished != finished) {
            std::cout << finished << "% ";
            prev_finished = finished;

            auto t2 = high_resolution_clock::now();
            auto ms_int = duration_cast<microseconds>(t2 - t1);
            duration<double, std::milli> ms_double = t1 - t2;

            auto now_ms = std::chrono::time_point_cast<std::chrono::microseconds>(t2);
            auto epoch = now_ms.time_since_epoch();
            auto value = std::chrono::duration_cast<std::chrono::microseconds>(epoch);
            long t2_ms = value.count();


            auto t = high_resolution_clock::now();
            std::cout << "Avg. us/cycle " << 1.0 * (t2_ms-t1_ms) / (i-prev_i) << std::endl;
            prev_i = i;
            
            now_ms = std::chrono::time_point_cast<std::chrono::microseconds>(t);
            epoch = now_ms.time_since_epoch();
            value = std::chrono::duration_cast<std::chrono::microseconds>(epoch);
            t1_ms = value.count();

        }


        Eigen::Vector2d pos = Eigen::Vector2d(distr(gen), distr(gen));
        Particle drop = Particle();
        drop.pos = pos;

        int ipos[2] = {static_cast<int>(drop.pos[0]), static_cast<int>(drop.pos[1])};

        while (drop.volume > world.min_volume) {
            ipos[0] = static_cast<int>(drop.pos[0]);
            ipos[1] = static_cast<int>(drop.pos[1]);
            Eigen::Vector2d n(surface_normal(ipos[0], ipos[1])(0), surface_normal(ipos[0], ipos[1])(2));
            drop.speed += world.dt * n / (drop.volume * world.density);
            drop.pos += world.dt * drop.speed;
            drop.speed *= (1 - world.dt * world.friction);

            if (drop.pos[0] < 1 || drop.pos[1] < 1 || drop.pos[0] >= 255 || drop.pos[1] >= 255) {
                break;
            }

            float max_sediment = drop.volume * std::sqrt(std::pow(drop.speed[0], 2) + std::pow(drop.speed[1], 2)) * (heightmap[ipos[0]][ipos[1]] - heightmap[static_cast<int>(drop.pos[0])][static_cast<int>(drop.pos[1])]);
            if (max_sediment > 0) {
                max_sediment = 0;
            }
            float sdiff = max_sediment - drop.sediment;

            drop.sediment += world.dt * world.deposition_rate * sdiff;
            heightmap[ipos[0]][ipos[1]] -= world.dt * drop.volume * world.deposition_rate * sdiff;
            drop.volume *= (1 - world.dt * world.evaporation_rate);
        }
    }
    std::cout << 100 << "% Done" << std::endl;
}

int main()
{
    

    World world = World();
    world.generate();

    std::fstream of_("originalmap.txt", std::ios::out | std::ios::trunc);

    for (int x = 0; x < 256; ++x) {
        for (int y = 0; y < 256; ++y) {
            of_ << world.heightmap[x][y] << " ";
        }
        of_ << "\n";
    }

    of_.close();

    world.erode(5000000);

    std::fstream of("Map.txt", std::ios::out | std::ios::trunc);

    for (int x = 0; x < 256; ++x) {
        for (int y = 0; y < 256; ++y) {
            of << world.heightmap[x][y] << " ";
        }
        of << "\n";
    }

    of.close();
}