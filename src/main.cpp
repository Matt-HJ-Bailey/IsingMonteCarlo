#include <iostream>
#include <random>
#include <vector>
#include "ising_grid.h"

int constexpr DEFAULT_ATOMS_PER_SIDE = 30;
double FIELD_STRENGTH = 0;
double COUPLING_PARAMETER = -1.0;
double constexpr BOLTZ = 1;
double constexpr TEMP = 0.01;
double constexpr BETA = 1.0 / (BOLTZ * TEMP);
int constexpr STEPS_PER_OUTPUT = 100;
int constexpr STEPS_PER_CALCULATE = 100;
int constexpr STEPS = 50000;

std::random_device seed;
std::default_random_engine generator(seed());
std::uniform_real_distribution<double> distr(0.0,1.0);


int main(int argc, char** argv) {
	std::cout << "Welcome to Ising Monte Carlo\n";
	int atoms_per_side = DEFAULT_ATOMS_PER_SIDE;
	if (argc == 2) {
		atoms_per_side = std::atoi(argv[1]);
	}
    std::cout << "Making a grid of " << atoms_per_side << " per side\n";

    std::uniform_int_distribution<int> int_distr(0,atoms_per_side);
    IsingGrid main_grid = IsingGrid(atoms_per_side, atoms_per_side);

    for (int step = 0; step < STEPS; ++step) {
        int x = int_distr(generator);
        int y = int_distr(generator);
        
        double energy_change = (2 * COUPLING_PARAMETER * main_grid.get(x, y) * (main_grid.get(x+1, y) + main_grid.get(x-1, y) + main_grid.get(x, y+1) + main_grid.get(x, y-1))) + (2 * FIELD_STRENGTH * main_grid.get(x, y));
        if (energy_change < 0) {
            main_grid.flip(x, y);
        } else {
            double probability    = std::exp(- energy_change * BETA);
            double random_chance  = distr(generator);
            if (random_chance < probability) {
                main_grid.flip(x, y);
            }
        }
        
        if (step % STEPS_PER_OUTPUT == 0) {
            std::cout << main_grid << std::endl;
        }
        if (step % STEPS_PER_CALCULATE == 0) {
            std::cout << main_grid.magnetisation() << "\n";
        }
    }
    std::cout << "Done. Cleaning up.\n";
    return 0;
}

