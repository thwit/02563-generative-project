#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#include <random>
#include <glm/glm.hpp>
#include <chrono>
// Libnoise
#include <libnoise/noise.h>

class World
{
public:
    void generate();
    void erode(int cycles);
    glm::vec3 surfaceNormal(int i, int j);

    double scale = 60.0;
    double heightmap[256][256] = {0.0};
    double stream[256][256] = {0.0};
    double pool[256][256] = {0.0};
    glm::ivec2 dim = glm::vec2(256);

    int SEED = 0;

    float dt = 1.2;
    float density = 1.0;
    float evapRate = 0.001;
    float depositionRate = 0.1;
    float minVol = 0.01;
    float friction = 0.05;
};

struct Particle
{
    Particle(glm::vec2 _pos) { pos = _pos; }
    void descend(World &world, std::vector<std::vector<bool>> &track);
    void flood(World &world);

    glm::vec2 pos;
    glm::vec2 speed = glm::vec2(0.0);

    float volume = 1.0;
    float sediment = 0.0;
    const double volumeFactor = 100.0; //"Water Deposition Rate"
};

void Particle::descend(World &world, std::vector<std::vector<bool>> &track)
{
    // As long as the droplet exists...
    while (volume > world.minVol)
    {
        glm::ivec2 ipos = pos; // Floored Droplet Initial Position
        track[ipos.x][ipos.y] = true;

        double effD = world.depositionRate;
        double effF = world.friction * (1.0 - 0.5 * world.stream[ipos.x][ipos.y]);
        double effR = world.evapRate * (1.0 - 0.2 * world.stream[ipos.x][ipos.y]);

        glm::vec3 n = world.surfaceNormal(ipos.x, ipos.y); // Surface Normal at Position

        // Accelerate particle using newtonian mechanics using the surface normal.
        glm::vec2 acc = glm::vec2(n.x, n.z) / (volume * world.density);
        speed += world.dt * acc; // F = ma, so a = F/m
        pos += world.dt * speed;
        speed *= (1.0 - world.dt * effF); // Friction Factor

        /*
          Note: For multiplied factors (e.g. friction, evaporation)
          time-scaling is correctly implemented like above.
        */

        // Check if Particle is still in-bounds
        if (!glm::all(glm::greaterThanEqual(pos, glm::vec2(0))) ||
            !glm::all(glm::lessThan(pos, glm::vec2(world.dim))))
        {
            volume = 0.;
            break;
        }

        // Slow-Down
        if (world.stream[ipos.x][ipos.y] > 0.5 && length(acc) < 0.01)
        {
            break;
        }

        // Enter Pool
        if (world.pool[ipos.x][ipos.y] > 0.0)
        {
            break;
        }

        // Compute sediment capacity difference
        float maxsediment = volume * glm::length(speed) * (world.heightmap[ipos.x][ipos.y] - world.heightmap[(int)pos.x][(int)pos.y]);
        if (maxsediment < 0.0)
            maxsediment = 0.0;
        float sdiff = maxsediment - sediment;

        // Act on the Heightmap and Droplet!
        sediment += world.dt * effD * sdiff;
        world.heightmap[ipos.x][ipos.y] -= world.dt * volume * world.depositionRate * sdiff;

        // Evaporate the Droplet (Note: Proportional to Volume! Better: Use shape factor to make proportional to the area instead.)
        volume *= (1.0 - world.dt * effR);
    }
}

void Particle::flood(World &world)
{
    glm::ivec2 ipos = pos;
    double plane = world.heightmap[ipos.x][ipos.y] + world.pool[ipos.x][ipos.y]; // Testing Plane
    double initialplane = plane;                                                 // Water Level

    // Flood Set
    std::vector<glm::vec2> set;
    int fail = 10; // Just in case...

    // Iterate while particle still has volume
    while (volume > world.minVol && fail)
    {

        set.clear();
        std::vector<std::vector<bool>> tried(world.dim.x,
                                             std::vector<bool>(world.dim.y, false));

        // Lowest Drain
        glm::ivec2 drain;
        bool drainfound = false;

        // Recursive Flood-Fill Function
        std::function<void(glm::ivec2)> fill = [&](glm::ivec2 v)
        {
            int i = v.x;
            int j = v.y;
            // Out of Bounds
            if (i >= world.dim.x || i < 0)
                return;
            if (j >= world.dim.y || j < 0)
                return;

            // Position has been tried
            if (tried[i][j])
                return;
            tried[i][j] = true;

            // Wall / Boundary of the Pool
            if (plane < world.heightmap[i][j] + world.pool[i][j])
                return;

            // Drainage Point
            if (initialplane > world.heightmap[i][j] + world.pool[i][j])
            {

                // No Drain yet
                if (!drainfound)
                    drain = v;

                // Lower Drain
                else if (world.pool[drain.x][drain.y] + world.heightmap[drain.x][drain.y] < world.pool[i][j] + world.heightmap[i][j])
                    drain = v;

                drainfound = true;
                return; // No need to flood from here
            }

            // Part of the Pool
            set.push_back(v);
            // Fill Neighbors
            fill(glm::vec2(i, j + 1));
            fill(glm::vec2(i, j - 1));
            fill(glm::vec2(i + 1, j));
            fill(glm::vec2(i - 1, j));
            // Diagonals (Improves Drainage)
            fill(glm::vec2(i - 1, j - 1));
            fill(glm::vec2(i - 1, j + 1));
            fill(glm::vec2(i + 1, j - 1));
            fill(glm::vec2(i + 1, j + 1));
        };

        // Perform Flood
        fill(ipos);

        // Drainage Point
        if (drainfound)
        {

            // Set the Particle Position
            pos = drain;

            // Set the New Waterlevel (Slowly)
            double drainage = 0.001;
            plane = (1.0 - drainage) * initialplane + drainage * (world.heightmap[drain.x][drain.y] + world.pool[drain.x][drain.y]);

            // Compute the New Height
            for (auto &s_ : set)
            { // Iterate over Set
                glm::ivec2 s = s_;
                world.pool[s.x][s.y] = (plane > world.heightmap[s.x][s.y]) ? (plane - world.heightmap[s.x][s.y]) : 0.0;
            }
            // Remove some sediment
            sediment *= 0.1;
            break;
        } // Get Volume under Plane
        double tVol = 0.0;
        for (auto &s_ : set)
        {
            glm::ivec2 s = s_;
            tVol += volumeFactor * (plane - (world.heightmap[s.x][s.y] + world.heightmap[s.x][s.y]));
        }

        // We can partially fill this volume
        if (tVol <= volume && initialplane < plane)
        {

            // Raise water level to plane height
            for (auto &s_ : set)
            {
                glm::ivec2 s = s_;
                world.pool[s.x][s.y] = plane - world.heightmap[s.x][s.y];
            }

            // Adjust Drop Volume
            volume -= tVol;
            tVol = 0.0;
        }

        // Plane was too high and we couldn't fill it
        else
        {
            fail--;
        }

        // Adjust Planes
        float approach = 0.5;
        initialplane = (plane > initialplane) ? plane : initialplane;
        plane += approach * (volume - tVol) / (double)set.size() / volumeFactor;
    }

    // Couldn't place the volume (for some reason)- so ignore this drop.
    if (fail == 0)
        volume = 0.0;

} // End of Flood Algorithm

void World::generate()
{
    std::cout << "Generating New World" << std::endl;
    SEED = time(NULL);
    std::cout << "Seed: " << SEED << std::endl;
    // Seed the Random Generator
    srand(SEED);

    std::cout << "... generating height ..." << std::endl;

    // Initialize Heightmap
    noise::module::Perlin perlin;
    perlin.SetOctaveCount(4);
    perlin.SetFrequency(1.0);
    perlin.SetPersistence(0.6);

    float min, max = 0.0;
    for (int i = 0; i < dim.x; i++)
    {
        for (int j = 0; j < dim.y; j++)
        {
            heightmap[i][j] = perlin.GetValue(i * (1.0 / dim.x), j * (1.0 / dim.y), SEED);
            if (heightmap[i][j] > max)
                max = heightmap[i][j];
            if (heightmap[i][j] < min)
                min = heightmap[i][j];
        }
    }

    // Normalize
    for (int i = 0; i < dim.x; i++)
    {
        for (int j = 0; j < dim.y; j++)
        {
            heightmap[i][j] = (heightmap[i][j] - min) / (max - min);
        }
    }
}

glm::vec3 World::surfaceNormal(int i, int j)
{
    /*
      Note: Surface normal is computed in this way, because the square-grid surface is meshed using triangles.
      To avoid spatial artifacts, you need to weight properly with all neighbors.
    */

    glm::vec3 n = glm::vec3(0.15) * glm::normalize(glm::vec3(scale * (heightmap[i][j] - heightmap[i + 1][j]), 1.0, 0.0)); // Positive X
    n += glm::vec3(0.15) * glm::normalize(glm::vec3(scale * (heightmap[i - 1][j] - heightmap[i][j]), 1.0, 0.0));          // Negative X
    n += glm::vec3(0.15) * glm::normalize(glm::vec3(0.0, 1.0, scale * (heightmap[i][j] - heightmap[i][j + 1])));          // Positive Y
    n += glm::vec3(0.15) * glm::normalize(glm::vec3(0.0, 1.0, scale * (heightmap[i][j - 1] - heightmap[i][j])));          // Negative Y

    // Diagonals! (This removes the last spatial artifacts)
    n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale * (heightmap[i][j] - heightmap[i + 1][j + 1]) / sqrt(2), sqrt(2), scale * (heightmap[i][j] - heightmap[i + 1][j + 1]) / sqrt(2))); // Positive Y
    n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale * (heightmap[i][j] - heightmap[i + 1][j - 1]) / sqrt(2), sqrt(2), scale * (heightmap[i][j] - heightmap[i + 1][j - 1]) / sqrt(2))); // Positive Y
    n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale * (heightmap[i][j] - heightmap[i - 1][j + 1]) / sqrt(2), sqrt(2), scale * (heightmap[i][j] - heightmap[i - 1][j + 1]) / sqrt(2))); // Positive Y
    n += glm::vec3(0.1) * glm::normalize(glm::vec3(scale * (heightmap[i][j] - heightmap[i - 1][j - 1]) / sqrt(2), sqrt(2), scale * (heightmap[i][j] - heightmap[i - 1][j - 1]) / sqrt(2))); // Positive Y

    return n;
}

void World::erode(int cycles)
{

    /*
      Note: Everything is properly scaled by a time step-size "dt"
    */

    // Do a series of iterations! (5 Particles)
    for (int i = 0; i < cycles; i++)
    { // Spawn Particle
        glm::vec2 newpos = glm::vec2(rand() % (int)dim.x, rand() % (int)dim.y);
        Particle drop(newpos);

        // Allow to enter a new pool 5 times before dying.. Prevents a potentially infinite loop
        int spill = 5;

        std::vector<std::vector<bool>> track(dim.x,
                                             std::vector<bool>(dim.y, false));

        while (drop.volume > minVol && spill != 0)
        {

            drop.descend(*this, track);

            if (drop.volume > minVol)
            {
                drop.flood(*this);
            }
            spill--;
        }
        // Update Path
        double lrate = 0.01; // Adaptation Rate
        for (int i = 0; i < dim.x; i++)
        {
            for (int j = 0; j < dim.y; j++)
            {
                stream[i][j] = (1.0 - lrate) * stream[i][j] + lrate * ((track[i][j]) ? 1.0 : 0.0);
            }
        }
    }
}

void saveMap(double hmap[][256], const char *fname)
{
    std::fstream of(fname, std::ios::out | std::ios::trunc);
    for (int x = 0; x < 256; ++x)
    {
        for (int y = 0; y < 256; ++y)
        {
            of << hmap[x][y] << " ";
        }
        of << "\n";
    }
    of.close();
}

int main(int argc, char *args[])
{
    World world;
    world.generate();
    saveMap(world.heightmap, "originalheightmap.txt");

    int stepsLeft = 500000;
    int stepSize = 10000;
    while (stepsLeft)
    {
        std::cout << stepsLeft << std::endl;
        world.erode(stepSize);
        stepsLeft -= stepSize;
    }

    saveMap(world.heightmap, "heightmap.txt");
    saveMap(world.stream, "streammap.txt");
    saveMap(world.pool, "poolmap.txt");

    return 0;
}