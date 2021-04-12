#pragma once
#include <vector>
#include <iostream>

extern std::default_random_engine generator;
extern std::uniform_real_distribution<double> distr;
extern double COUPLING_PARAMETER;
extern double FIELD_STRENGTH;

class IsingGrid {
    private:
        std::vector<int> _grid;
    public:
        int m_height;
        int m_width;
        int get(int x, int y) const;
        void flip(int x, int y);
        double magnetisation() const;
        double magnetisation_squared() const;
        double energy() const;
        double energy_squared() const;
        IsingGrid(int m_height_, int m_width_);
};

std::ostream& operator<<(std::ostream& os, const IsingGrid& ising_grid);