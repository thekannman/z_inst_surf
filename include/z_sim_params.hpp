//Copyright (c) 2015 Zachary Kann
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

// ---
// Author: Zachary Kann

#include <iostream>
#include <fstream>
#include <string>
#include "z_constants.hpp"
#include "xdrfile_xtc.h"
#include <armadillo>

using namespace arma;

#ifndef Z_LIB_Z_SIMPARAMS_HPP_
#define Z_LIB_Z_SIMPARAMS_HPP_

enum Comb_rule {C6C12, LB, GEOMETRIC};

//Provides convenient access to simulation parameters as well as the ability
//to read those parameters from an intuitive parameter file.
class Sim_params
{
    public:
        unsigned numMols, numChromos, wAtoms,
            numCations, numAnions,
            numIons, numOthers;
        double boxLength, zboxLength, zmult;
        double rminOO, rminOO2;
        double rminC, rminA, rminAH, rminC2, rminA2, rminAH2,
            rmin2C, rmin2A, rmin2C2, rmin2A2,
            rmin3C, rmin3A, rmin3C2, rmin3A2;
        double gamm, gammC, gammA;

        inline void set_filename(const std::string& filename)
        {
            filename_ = filename;
        }
        void read_params();
        inline void read_params(const std::string& filename)
        {
            set_filename(filename);
            read_params();
        }
        inline void set_temperature(const double temperature)
        {
            temperature_ = temperature;
            kT_ = KB*temperature;
            beta_ = 1.0/kT_;
        }
        inline double temperature() const { return temperature_;}
        inline double kT() const { return kT_;}
        inline double beta() const { return beta_;}
        inline double dt() const { return dt_; }
        inline unsigned steps() const { return steps_; }
        inline int num_atoms() const { return num_atoms_; }
        inline unsigned max_steps() const { return max_steps_; }
        inline unsigned max_time() const { return (max_steps_*dt_); }
        inline void set_max_time(const double max_time)
        {
            max_steps_ = max_time / dt_;
        }
        inline void set_comb_rule(const Comb_rule& comb_rule)
        {
            comb_rule_ = comb_rule;
        }
        void extract_traj_metadata(char *traj, rvec **x_in, rowvec& box);

    private:
        std::string filename_;
        std::ifstream file_;
        double temperature_;
        double kT_;
        double beta_;
        double dt_;
        unsigned steps_;
        int num_atoms_;
        Comb_rule comb_rule_;
        unsigned max_steps_;

};
#endif
