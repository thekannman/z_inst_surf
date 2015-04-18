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

#include "z_atom.hpp"
#include "z_molecule.hpp"
#include "z_gromacs.hpp"
#include "z_string.hpp"
#include <cassert>
#include <boost/lexical_cast.hpp>

using namespace arma;

void read_atom_types(std::ifstream& top_file, std::vector<Atom>& atom_types)
{
    std::string line;
    while (getline(top_file, line))
    {
        if(line == "[ atomtypes ]") break;
    }
    while (getline(top_file, line))
    {
        if (line[0] == ';') continue;
        if (line.empty() || line[0] == '[') break;
        Atom atom_type;
        std::vector<std::string> split_line = split(line, ' ');
        atom_type.set_type(split_line[0]);
        atom_type.set_mass(boost::lexical_cast<double>(split_line[2]));
        atom_type.set_sigma(boost::lexical_cast<double>(split_line[5]));
        atom_type.set_epsilon(boost::lexical_cast<double>(split_line[6]));
        atom_types.push_back(atom_type);
    }
    do
    {
        if(line == "[ nonbond_params ]") break;
    } while (getline(top_file, line));
    while (getline(top_file, line))
    {
        if (line[0] == ';') continue;
        if (line.empty() || line[0] == '[') break;
        std::vector<std::string> split_line = split(line, ' ');
        for (std::vector<Atom>::iterator i_type = atom_types.begin();
            i_type != atom_types.end(); i_type++)
        {
            if ((*i_type).name() == split_line[0])
            {
                (*i_type).set_cross_lj(split_line[1],
                    boost::lexical_cast<double>(split_line[2]),
                    boost::lexical_cast<double>(split_line[3]));
            }
            else if ((*i_type).name() == split_line[1])
            {
                (*i_type).set_cross_lj(split_line[0],
                    boost::lexical_cast<double>(split_line[2]),
                    boost::lexical_cast<double>(split_line[3]));
            }
        }
    }
}

std::vector<Molecule> gen_molecules(const std::string& top, Sim_params& params)
{
    std::ifstream top_file;
    std::string line;
    top_file.open(top.c_str());
    std::vector<Molecule> molecules;
    std::vector<Atom> atom_types;

    assert(top_file.is_open());
    while (getline(top_file, line))
    {
        if(line == "[ defaults ]") break;
    }
    while (getline(top_file, line))
    {
        if (line[0] == ';') continue;
        if (line.empty() || line[0] == '[') break;
        std::vector<std::string> split_line = split(line, ' ');
        params.set_comb_rule(
            static_cast<Comb_rule>(boost::lexical_cast<int>(split_line[1])));
        break;
    }
    read_atom_types(top_file, atom_types);
    while (!top_file.eof())
    {
        do
        {
            if(line == "[ moleculetype ]") break;
        } while (getline(top_file, line));
        if (top_file.eof()) break;
        Molecule molecule;
        while (getline(top_file, line))
        {
            if (line[0] == ';') continue;
            assert(!line.empty() && line[0] != '[');
            std::vector<std::string> split_line = split(line, ' ');
            molecule.set_name(split_line[0]);
            break;
        }
        while (getline(top_file, line))
        {
            if(line == "[ atoms ]") break;
        }
        while (getline(top_file, line))
        {
            if (line[0] == ';') continue;
            if (line.empty() || line[0] == '[') break;
            std::vector<std::string> split_line = split(line, ' ');
            Atom atom;
            for (std::vector<Atom>::iterator i_type = atom_types.begin();
                i_type != atom_types.end(); i_type++)
            {
                if ((*i_type).type() == split_line[1])
                {
                    atom = *i_type;
                    break;
                }
            }
            atom.set_name(split_line[4]);
            atom.set_charge(boost::lexical_cast<double>(split_line[6]));
            molecule.add_atom(atom);
        }
        molecules.push_back(molecule);
    }
    return molecules;
}

std::map<std::string, std::vector<unsigned> > read_ndx(const std::string& ndx)
{
    std::map<std::string, std::vector<unsigned> > groups;
    std::ifstream ndx_file;
    std::string line;
    std::string subline, group_name;
    std::vector<unsigned> group;

    ndx_file.open(ndx.c_str());
    while (getline(ndx_file, line))
    {
        std::istringstream iline(line);
        iline >> subline;
        if(subline == "[")
        {
            if (!group.empty())
            {
                groups[group_name] = group;
                group.clear();
            }
            iline >> group_name;
        }
        else
        {
            do
            {
                group.push_back(boost::lexical_cast<int>(subline.c_str())-1);
            } while (iline >> subline);
        }
    }
    if (!group.empty())
    {
        groups[group_name] = group;
    }
    ndx_file.close();

    return groups;
}

std::vector<unsigned> select_group(
    std::map<std::string, std::vector<unsigned> >& groups,
    const std::string& group_name)
{
    std::map<std::string, std::vector<unsigned> >::iterator i_group =
        groups.find(group_name);
    if (i_group == groups.end())
    {
        std::cout << "Group name '" << group_name <<
            "' not found. Try one of:" << std::endl;
        for(i_group = groups.begin(); i_group != groups.end(); ++i_group)
        {
            std::cout << i_group->first << std::endl;
        }
        exit(EXIT_FAILURE);
    }
    return groups[group_name];
}

rowvec rvec_to_row(const rvec& vec1)
{
    rowvec result(3);
    for (unsigned i=0; i<3; i++)
        result(i) = vec1[i];
    return result;
}
