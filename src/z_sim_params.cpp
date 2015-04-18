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

#include "z_sim_params.hpp"
#include "z_string.hpp"
#include "z_conversions.hpp"
#include <vector>
#include <cassert>
#include <boost/lexical_cast.hpp>

void Sim_params::read_params()
{
    std::string line;
    std::string subline;
    file_.open(filename_.c_str());
    while (getline(file_, line))
    {
        if(line[0] == '#' || line[0] == '@') continue;
        std::vector<std::string> split_line = split(line, '=');
        std::string variable = trim(split_line[0]);
        std::string value = trim(split_line[1]);
        assert(!variable.empty() && !value.empty());
        if (variable == "boxLength")
            boxLength = boost::lexical_cast<double>(value.c_str());
        else if (variable == "zboxLength")
            zboxLength = boost::lexical_cast<double>(value.c_str());
        else if (variable == "numMols")
            numMols = boost::lexical_cast<int>(value.c_str());
        else if (variable == "numCations")
            numCations = boost::lexical_cast<int>(value.c_str());
        else if (variable == "numAnions")
            numAnions = boost::lexical_cast<int>(value.c_str());
        else if (variable == "numOthers")
            numOthers = boost::lexical_cast<int>(value.c_str());
        else if (variable == "numSteps")
            max_steps_ = boost::lexical_cast<int>(value.c_str());
        else if (variable == "rminOO")
            rminOO = boost::lexical_cast<double>(value.c_str());
        else if (variable == "rminC")
            rminC = boost::lexical_cast<double>(value.c_str());
        else if (variable == "rminA")
            rminA = boost::lexical_cast<double>(value.c_str());
        else if (variable == "rminAH")
            rminAH = boost::lexical_cast<double>(value.c_str());
        else if (variable == "rmin2C")
            rmin2C = boost::lexical_cast<double>(value.c_str());
        else if (variable == "rmin2A")
            rmin2A = boost::lexical_cast<double>(value.c_str());
        else if (variable == "rmin3C")
            rmin3C = boost::lexical_cast<double>(value.c_str());
        else if (variable == "rmin3A")
            rmin3A = boost::lexical_cast<double>(value.c_str());
        else if (variable == "dt")
            dt_ = boost::lexical_cast<double>(value.c_str());
        else if (variable == "temp")
            set_temperature(boost::lexical_cast<double>(value.c_str()));
        else if (variable == "gamm")
            gamm = boost::lexical_cast<double>(value.c_str())*THZ_CONVERT;
        else if (variable == "gammC")
            gammC = boost::lexical_cast<double>(value.c_str())*THZ_CONVERT;
        else if (variable == "gammA")
            gammA = boost::lexical_cast<double>(value.c_str())*THZ_CONVERT;
//        else if (variable == "hist size")
//            hist.size = boost::lexical_cast<int>(value.c_str());
//        else if (variable == "hist interval")
//            hist.interval = boost::lexical_cast<double>(value.c_str());
        else
        {
            std::cerr << "Param option not recognized" << std::endl;
        }
    }
    numChromos = 2.0*(numMols);
    numIons = (numCations)+(numAnions);
    rminOO2 = (rminOO)*(rminOO);
    rminC2 = (rminC)*(rminC);
    rminA2 = (rminA)*(rminA);
    rminAH2 = (rminAH)*(rminAH);
    rmin2C2 = (rmin2C)*(rmin2C);
    rmin2A2 = (rmin2A)*(rmin2A);
    rmin3C2 = (rmin3C)*(rmin3C);
    rmin3A2 = (rmin3A)*(rmin3A);
    if(zboxLength)
        zmult = (zboxLength)/(boxLength);
    file_.close();
}

void Sim_params::extract_traj_metadata(char *traj, rvec **x_in, rowvec& box)
{
    int st;
    matrix box_mat;
    float t1, t2, prec;
    XDRFILE *traj_file;
    read_xtc_natoms(traj, &num_atoms_);
    traj_file = xdrfile_open(traj, "r");
    *x_in = new rvec [num_atoms_];
    read_xtc(traj_file, num_atoms_, &st, &t1, box_mat, *x_in, &prec);
    read_xtc(traj_file, num_atoms_, &st, &t2, box_mat, *x_in, &prec);
    dt_ = t2 - t1;
    for(unsigned i=0; i<DIMS; i++) box(i) = box_mat[i][i];
    xdrfile_close(traj_file);
}
