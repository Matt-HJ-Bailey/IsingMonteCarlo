#include <vector>
#include <random>
#include "ising_grid.h"



IsingGrid::IsingGrid(int m_height_, int m_width_) {
    // Constructs the lattice of size (m_height_, m_width_)
    // and assigns values of +1 or -1 to a vector.
    m_height = m_height_;
    m_width = m_width_;
    _grid.reserve(m_height * m_width);
    for (int i = 0; i < (m_height * m_width); ++i) {
        auto random = distr(generator);
        if (random < 0.5) {
            _grid.push_back(1);
        } else {
            _grid.push_back(-1);
        }
    }
 }


int IsingGrid::get(int x, int y) const {
    // Retrieves the value at the lattice site (x, y)
    // applying periodic boundary conditions to it.
    if (x > m_width) {
        x = x - m_width;
    } else if (x < 0) {
        x = x + m_width;
    }
    
    if (y > m_height) {
        y = y - m_height;
    } else if (y < 0) {
        y = y + m_height;
    }
    return _grid[(y * m_height) + x];
}

void IsingGrid::flip(int x, int y) {
    // Flips the value at the lattice site (x, y)
    // from +1 to -1 or vice-versa,
    // applying periodic boundary conditions to it.
    if (x > m_width) {
        x = x - m_width;
    } else if (x < 0) {
        x = x + m_width;
    }
    
    if (y > m_height) {
        y = y - m_height;
    } else if (y < 0) {
        y = y + m_height;
    }
    _grid[(y * m_height) + x] = _grid[(y * m_height) + x] * -1;
}

double IsingGrid::magnetisation() const {
    // Calculates the magnetisation of the grid,
    // using the equation:
    // <M> = 1/N Sum_N m_i
    int excess = 0;
    for (auto spin: _grid) {
        excess += spin;
    }
    return excess / static_cast<double>(m_height * m_width);
}

double IsingGrid::magnetisation_squared() const {
    // Calculates the magnetisation of the grid,
    // using the equation:
    // <M> = 1/N Sum_N m_i^2
    int excess = 0;
    for (auto spin: _grid) {
        excess += (spin * spin);
    }
    return excess / static_cast<double>(m_height * m_width);
}

double IsingGrid::energy() const {
    // Calculates the magnetisation of the grid,
    // using the equation:
    // M = 1/2 Sum_N H_i
    double energy_temp = 0.0;
    for (int y = 0; y < m_height; ++y) {
        for (int x = 0; x < m_width; ++x) {
            energy_temp += - (COUPLING_PARAMETER * get(x, y) * (get(x-1, y) + get(x+1, y) + get(x, y-1) + get(x, y+1))) - FIELD_STRENGTH * get(x, y);
        }  
    }
    return 0.5 * energy_temp;
}

double IsingGrid::energy_squared() const {
    double energy = 0.0;
    for (int y = 0; y < m_height; ++y) {
        for (int x = 0; x < m_width; ++x) {
            double energy_temp = - (COUPLING_PARAMETER * get(x, y) * (get(x-1, y) + get(x+1, y) + get(x, y-1) + get(x, y+1))) - FIELD_STRENGTH * get(x, y);
            energy += energy_temp * energy_temp;
        }  
    }
    return 0.5 * energy;
}

std::ostream& operator<<(std::ostream& os, const IsingGrid& ising_grid) {
    for (int y = 0; y < ising_grid.m_height; ++y) {
        for (int x = 0; x < ising_grid.m_width; ++x) {
            switch (ising_grid.get(x, y))
            {
                case -1:
                {
                    os << "█";
                    break;
                }
                case 1:
                {
                    os << "░";
                    break;
                }
                default:
                {
                    os << "!";
                    break;
                }
            }
        }
        os << "\n";
    }
    return os;
}