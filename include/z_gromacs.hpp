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

#include "xdrfile_xtc.h"
#include <map>
#include <armadillo>
#include "z_sim_params.hpp"

using namespace arma;

#ifndef _Z_GROMACS_HPP_
#define _Z_GROMACS_HPP_


extern void read_atom_types(std::ifstream& top_file,
    std::vector<Atom>& atom_types);
extern std::vector<Molecule> gen_molecules(const std::string& top,
    Sim_params& params);
extern std::map<std::string, std::vector<unsigned> > read_ndx(
    const std::string& ndx);
extern std::vector<unsigned> select_group(
    std::map<std::string, std::vector<unsigned> >& groups,
    const std::string& group_name);
extern rowvec rvec_to_row(const rvec& vec1);
#endif
