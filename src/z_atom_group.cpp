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

#include "z_atom_group.hpp"
#include <cassert>
#include <boost/lexical_cast.hpp>
#include "z_string.hpp"
#include <boost/format.hpp>

using namespace arma;

void Atom_group::set_indices(std::string name, std::vector<unsigned> indices,
    std::vector<double> index_to_mass,
    std::vector<unsigned> index_to_molecule)
{
    name_ = name;
    indices_ = indices;
    group_size_ = indices_.size();
    positions_ = zeros(group_size_, DIMS);
    velocities_ = zeros(group_size_, DIMS);
    for (std::vector<unsigned>::iterator i_index = indices_.begin();
        i_index != indices_.end(); i_index++)
    {
        index_to_mass_.push_back(index_to_mass[*i_index]);
        index_to_molecule_.push_back(index_to_molecule[*i_index]);
    }
    std::vector<double>::iterator i_mass = index_to_mass_.begin();
    for (std::vector<unsigned>::iterator i_mol = index_to_molecule_.begin();
        i_mol != index_to_molecule_.end(); i_mol++, i_mass++)
    {
        std::vector<unsigned>::iterator mol_index;
        mol_index = std::find(molecular_index_to_molecule_.begin(),
            molecular_index_to_molecule_.end(), *i_mol);
        if ( mol_index == molecular_index_to_molecule_.end())
        {
            molecular_index_to_mass_.push_back(0.0);
            index_to_molecular_index_.push_back(molecular_index_to_molecule_.size());
            molecular_index_to_molecule_.push_back(*i_mol);
        }
        else
            index_to_molecular_index_.push_back(std::distance(molecular_index_to_molecule_.begin(),mol_index));
        molecular_index_to_mass_.back() += *i_mass;
    }
    num_molecules_ = molecular_index_to_molecule_.size();
    positions_.set_size(group_size_, DIMS);
    velocities_.set_size(group_size_, DIMS);
    com_positions_.set_size(num_molecules_, DIMS);
    com_velocities_.set_size(num_molecules_, DIMS);
}
void Atom_group::zero_com()
{
    com_positions_ = zeros(num_molecules_, DIMS);
    com_velocities_ = zeros(num_molecules_, DIMS);
}

void Atom_group::remove_molecule(const unsigned molecule_id, rowvec& position, rowvec& velocity)
{
    std::vector<unsigned> mark_for_rem;
    for (std::vector<unsigned>::iterator i_mol = index_to_molecule_.begin();
        i_mol != index_to_molecule_.end(); i_mol++)
    {
        if (*i_mol == molecule_id)
            mark_for_rem.push_back(std::distance(index_to_molecule_.begin(),i_mol));
    }
    assert(!mark_for_rem.empty());
    molecular_index_to_mass_.erase(molecular_index_to_mass_.begin() + index_to_molecular_index_[mark_for_rem.front()]);
    molecular_index_to_molecule_.erase(molecular_index_to_molecule_.begin() + index_to_molecular_index_[mark_for_rem.front()]);
    com_positions_.shed_row(index_to_molecular_index_[mark_for_rem.front()]);
    com_velocities_.shed_row(index_to_molecular_index_[mark_for_rem.front()]);
    for (std::vector<unsigned>::reverse_iterator i_atom = mark_for_rem.rbegin();
        i_atom != mark_for_rem.rend(); i_atom++)
    {
        indices_.erase(indices_.begin() + *i_atom);
        index_to_mass_.erase(index_to_mass_.begin() + *i_atom);
        index_to_molecule_.erase(index_to_molecule_.begin() + *i_atom);
        index_to_molecular_index_.erase(index_to_molecular_index_.begin() + *i_atom);
        positions_.shed_row(*i_atom);
        velocities_.shed_row(*i_atom);
        group_size_--;
    }
    num_molecules_--;
}

void Atom_group::remove_molecule(const unsigned molecule_id)
{
    rowvec position(DIMS);
    rowvec velocity(DIMS);
    remove_molecule(molecule_id, position, velocity);
}
void Atom_group::add_molecule(Molecule mol_to_add,
    const rowvec& position, const rowvec& velocity)
{
    //ADD FUNCTION!!!
}
void Atom_group::replace_molecule(const unsigned mol_to_remove,
    Molecule mol_to_add)
{
    rowvec position (DIMS);
    rowvec velocity(DIMS);
    remove_molecule(mol_to_remove, position, velocity);
    add_molecule(mol_to_add, position, velocity);
}

void Atom_group::write_gro(const std::string& gro,
    rowvec box, std::string description)
{
    std::ofstream gro_file;
    gro_file.open(gro.c_str());
    std::cout << description << std::endl;
    std::cout << group_size_ << std::endl;
    for (unsigned i_atom = 0; i_atom < group_size_; i_atom++)
    {
        std::cout << boost::format("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f") %
            index_to_molecular_index_[i_atom] %
            molecule_names_[i_atom] % atom_names_[i_atom] %
            indices_[i_atom] % positions_(i_atom,0) %
            positions_(i_atom,1) % positions_(i_atom,2) %
            velocities_(i_atom,0) % velocities_(i_atom,1) %
            velocities_(i_atom,2) << std::endl;
    }
    std::cout << box(0) << " " << box(1) << " " << box(2);
}

void Atom_group::update_com()
{
  // Actually com_position*mass and com_velocity*mass
    zero_com();
    for (unsigned i_atom = 0; i_atom < group_size_; i_atom++)
    {
        com_positions_.row(index_to_molecular_index_[i_atom]) +=
            index_to_mass_[i_atom]*positions_.row(i_atom);
        com_velocities_.row(index_to_molecular_index_[i_atom]) +=
            index_to_mass_[i_atom]*velocities_.row(i_atom);
    }
}

void Atom_group::read_gro(const std::string& gro,
    const std::vector<Molecule> molecules)
{
    std::ifstream gro_file;
    std::string line;
    gro_file.open(gro.c_str());
    getline(gro_file, line);
    getline(gro_file, line);
    std::vector<double> index_to_mass;
    std::vector<unsigned> index_to_molecule;
    const unsigned num_atoms = boost::lexical_cast<int>(trim(line));
    std::vector<unsigned> indices;
    Molecule mol_match;
    positions_.zeros(num_atoms,DIMS);
    velocities_.zeros(num_atoms,DIMS);
    unsigned mol_number = std::numeric_limits<unsigned>::max(), mol_number_old;
    for (unsigned i_atom = 0; i_atom < num_atoms; i_atom++)
    {
        indices.push_back(i_atom);
        std::string mol_name;
        assert(!gro_file.eof());
        getline(gro_file, line);
        std::vector<std::string> split_line = split(line, ' ');
        std::istringstream iline(split_line[0]);
        mol_number_old = mol_number;
        iline >> mol_number >> mol_name;
        mol_number--;
        atom_names_.push_back(split_line[1]);
        if (mol_number_old != mol_number)
            molecule_names_.push_back(mol_name);
        index_to_molecule.push_back(mol_number);
        for (std::vector<Molecule>::const_iterator i_mol = molecules.begin();
            i_mol != molecules.end(); i_mol++)
        {
            if (mol_name == (*i_mol).name())
            {
                mol_match = *i_mol;
                break;
            }
        }
        for (std::vector<Atom>::const_iterator i_at = mol_match.atoms().begin();
            i_at != mol_match.atoms().end(); i_at++)
        {
            (*i_at).print();
            if (split_line[1] == (*i_at).name())
            {
                index_to_mass.push_back((*i_at).mass());
                break;
            }
        }
        for (unsigned i_dim = 0; i_dim<DIMS; i_dim++)
        {
            positions_(i_atom,i_dim) =
                boost::lexical_cast<double>(split_line[3+i_dim]);
        }
        if (split_line.size() == 9)
            for (unsigned i_dim = 0; i_dim<DIMS; i_dim++)
            {
                velocities_(i_atom,i_dim) =
                    boost::lexical_cast<double>(split_line[6+i_dim]);
            }
    }
    set_indices("all", indices, index_to_mass, index_to_molecule);
}
