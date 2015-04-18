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

#include "z_vec.hpp"

using namespace arma;

Grid::Grid(rowvec box, double grid_spacing)
{
    length = box;
    spacing = grid_spacing;
    points.resize(dims);
    for (unsigned  i=0; i<dims; i++)
        points(i) = length(i)/spacing;
    upper_mesh = zeros<cube>(points(0),points(1),dims);
    for (unsigned  i_x=0; i_x<points(0); i_x++)
    {
        for (unsigned  i_y=0; i_y<points(0); i_y++)
        {
            upper_mesh(i_x,i_y,0) = grid_spacing*(i_x+0.5);
            upper_mesh(i_x,i_y,1) = grid_spacing*(i_y+0.5);
        }
    }
    lower_mesh = upper_mesh;
    grid.resize(boost::extents[points(0)][points(1)][points(2)][dims]);
    for (unsigned  i_x=0; i_x<points(0); i_x++)
        for (unsigned  i_y=0; i_y<points(1); i_y++)
            for (unsigned  i_z=0; i_z<points(2); i_z++)
            {
                grid[i_x][i_y][i_z][0] = grid_spacing*(i_x+0.5);
                grid[i_x][i_y][i_z][1] = grid_spacing*(i_y+0.5);
                grid[i_x][i_y][i_z][2] = grid_spacing*(i_z+0.5);
            }
}
void Grid::zero_rho()
{
    rho = zeros<cube>(points(0),points(1),points(2));
}
void Grid::read_mesh(std::fstream &mesh_file)
{
    for (unsigned  i_x=0; i_x<points(0); i_x++)
    {
        for (unsigned  i_y=0; i_y<points(1); i_y++)
        {
            // Skip x,y coords. May remove and rewrite surf files to contain only z coords and grid spacing.
            mesh_file >> upper_mesh(i_x,i_y,2);// >> lower_mesh(i_x,i_y,2);
            mesh_file >> lower_mesh(i_x,i_y,2);
        }
    }
}
double Grid::avg_mesh()
{
    // Only upper for now
    double avg;
    for (unsigned  i_x=0; i_x<points(0); i_x++)
    {
        for (unsigned  i_y=0; i_y<points(1); i_y++)
        {
            avg += upper_mesh(i_x,i_y,2);
        }
    }
    avg /= (points(0)*points(1));
    return avg;
}

void setupDX(cube& xoo, cube& xio, cube& xii, mat& roo2, mat& rio2, mat& rii2, const mat& x, const mat& xion, icube& shiftoo,
                    icube& shiftio, icube& shiftii, const unsigned numMols, const unsigned numIons, const rowvec& box)
{
    rowvec x_row(DIMS);
    irowvec shift_row(DIMS);

    for (unsigned i1=0; i1<numMols; i1++)
        for (unsigned i2=i1+1; i2<numMols; i2++)
        {
            findDx(x_row, x.row(i1), x.row(i2), box, shift_row);
            xoo.tube(i1,i2) = x_row;
            shiftoo.tube(i1,i2) = shift_row;
            xoo.tube(i2,i1) = -1.0*x_row;
            shiftoo.tube(i2,i1) = -1*shift_row;
            roo2(i2,i1) = roo2(i1,i2) = dot(x_row,x_row);
        }
    for (unsigned i1=0; i1<numIons; i1++)
    {
        for (unsigned i2=0; i2<numMols; i2++)
        {
            findDx(x_row, xion.row(i1), x.row(i2), box, shift_row);
            xio.tube(i1,i2) = x_row;
            shiftio.tube(i1,i2) = shift_row;
            rio2(i2,i1) = rio2(i1,i2) = dot(x_row,x_row);
        }
        for (unsigned i2=i1+1; i2<numIons; i2++)
        {
            findDx(x_row, xion.row(i1), xion.row(i2), box, shift_row);
            xii.tube(i1,i2) = x_row;
            shiftii.tube(i1,i2) = shift_row;
            xii.tube(i2,i1) = -1.0*x_row;
            shiftii.tube(i2,i1) = -1*shift_row;
            rii2(i2,i1) = rii2(i1,i2) = dot(x_row,x_row);
        }
    }
}

void setupDX_io(cube& xio, const mat& x, const mat& xion, icube& shiftio, const unsigned numMols, const unsigned numIons, const rowvec& box)
{
    rowvec x_row(DIMS);
    irowvec shift_row(DIMS);

    for (unsigned i1=0; i1<numIons; i1++)
        for (unsigned i2=0; i2<numMols; i2++)
        {
            findDx(x_row, xion.row(i1), x.row(i2), box, shift_row);
            xio.tube(i1,i2) = x_row;
            shiftio.tube(i1,i2) = shift_row;
        }
}

void setupShift_oo(const mat& x, icube& shiftoo, const unsigned numMols, const rowvec& box)
{
    rowvec x_row(DIMS);
    irowvec shift_row(DIMS);

    for (unsigned i1=0; i1<numMols; i1++)
        for (unsigned i2=i1+1; i2<numMols; i2++)
        {
            findDx(x_row, x.row(i1), x.row(i2), box, shift_row);
            shiftoo.tube(i1,i2) = shift_row;
            shiftoo.tube(i2,i1) = -1*shift_row;
        }
}

void findDx_noShift(rowvec& result, const rowvec& vec1, const rowvec& vec2, const rowvec& box)
{
    result = vec1 - vec2;
    for (unsigned i1=0; i1<DIMS; i1++)
    {
        if (result(i1)>box(i1)/2.0)
            result(i1) -= box(i1);
        else if (result(i1)<-box(i1)/2.0)
            result(i1) += box(i1);
    }
}

void findDx(rowvec& result, const rowvec& vec1, const rowvec& vec2, const rowvec& box, irowvec& shift)
{
    result = vec1 - vec2;
    for (unsigned i1=0; i1<DIMS; i1++)
    {
        if (result(i1)>box(i1)/2.0)
        {
            result(i1) -= box(i1);
            shift(i1) = -1;
        }
        else if (result(i1)<-box(i1)/2.0)
        {
            result(i1) += box(i1);
            shift(i1) = 1;
        }
        else shift(i1) = 0;
    }
}
