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

#include "boost/multi_array.hpp"
#include <armadillo>
#include "z_constants.hpp"


using namespace arma;

typedef boost::multi_array<double, 4> hypercube;

#ifndef _Z_VEC_HPP_
#define _Z_VEC_HPP_

class Grid
{
    private:
        static const int dims = 3;
        rowvec length;
        urowvec points;

    public:
        double spacing;
        cube rho;
        hypercube grid;
        cube upper_mesh;
        cube lower_mesh;
        Grid(rowvec box, double grid_spacing = 0.1);
        void zero_rho();
        void read_mesh(std::fstream &mesh_file);
        double avg_mesh();
        unsigned  size_x() const { return points(0); }
        unsigned  size_y() const { return points(1); }
        unsigned  size_z() const { return points(2); }
        urowvec size() const { return points; }
};

extern void setupDX(cube& xoo, cube& xio, cube& xii, mat& roo2, mat& rio2, mat& rii2, const mat& x, const mat& xion, icube& shiftoo,
                    icube& shiftio, icube& shiftii, const unsigned  numMols, const unsigned  numIons, const rowvec& box);
extern void setupDX_io(cube& xio, const mat& x, const mat& xion, icube& shiftio, const unsigned  numMols, const unsigned  numIons, const rowvec& box);
extern void setupShift_oo(const mat& x, icube& shiftoo, const unsigned  numMols, const rowvec& box);
extern void findDx_noShift(rowvec& result, const rowvec& vec1, const rowvec& vec2, const rowvec& box);
extern void findDx(rowvec& result, const rowvec& vec1, const rowvec& vec2, const rowvec& box, irowvec& shift);

template <typename type1>
rowvec veclike_to_row(const type1& vec1, const unsigned  length)
{
    rowvec result(length);
    for (int i=0; i<length; i++)
    {
        result[i] = vec1[i];
    }
    return result;
}

inline rowvec useDx(const rowvec& vec1, const rowvec& vec2, const rowvec& box, const irowvec& shift) { return vec1 - vec2 + shift%box; }

#endif
