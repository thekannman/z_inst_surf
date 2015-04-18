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

#include <armadillo>
#include "z_molecule.hpp"
#include "z_constants.hpp"
#include "z_gromacs.hpp"

using namespace arma;

#ifndef _Z_ATOM_GROUP_HPP_
#define _Z_ATOM_GROUP_HPP_

class Atom_group
{
    public:
        Atom_group(): num_molecules_(0), group_size_(0) {}
        inline Atom_group(const std::string& gro,
            const std::vector<Molecule> molecules):
            num_molecules_(0), group_size_(0)
        {
            read_gro(gro, molecules);
        }

        inline std::string const name() { return name_; }
        void set_indices(std::string name, std::vector<unsigned> indices,
            std::vector<double> index_to_mass,
            std::vector<unsigned> index_to_molecule);
        inline void set_indices(std::string name, std::vector<unsigned> indices,
            const Atom_group& all_atoms)
        {
            set_indices(name, indices, all_atoms.index_to_mass_,
                all_atoms.index_to_molecule_);
        }
        void zero_com();
        void remove_molecule(const unsigned molecule_id, rowvec& position,
            rowvec& velocity);
        void remove_molecule(const unsigned molecule_id);
        void add_molecule(Molecule mol_to_add,
            const rowvec& position = zeros<rowvec>(DIMS),
            const rowvec& velocity = zeros<rowvec>(DIMS));
        void replace_molecule(const unsigned mol_to_remove, Molecule mol_to_add);
        std::vector<unsigned>::iterator begin() { return indices_.begin(); }
        std::vector<unsigned>::iterator end() { return indices_.end(); }
        inline unsigned size() const { return group_size_; }
        inline std::vector<unsigned> indices() const { return indices_; }
        inline unsigned indices(int i) const { return indices_[i]; }
        void write_gro(const std::string& gro,
            rowvec box, std::string description);
        inline double num_molecules() const { return num_molecules_; }
        inline void set_position(const unsigned i, const rvec& position)
        {
            positions_.row(i) = rvec_to_row(position);
        }
        inline void set_velocity(const unsigned i, const rvec& velocity)
        {
            velocities_.row(i) = rvec_to_row(velocity);
        }
        void update_com();
        inline mat positions() const { return positions_; }
        inline mat velocities() const { return velocities_; }
        inline mat com_positions() const { return com_positions_; }
        inline mat com_velocities() const { return com_velocities_; }
        inline rowvec position(unsigned atom) const
        {
            return positions_.row(atom);
        }
        inline rowvec velocity(unsigned atom) const
        {
            return velocities_.row(atom);
        }
        inline rowvec com_position(unsigned molecule) const
        {
            return com_positions_.row(molecule);
        }
        inline rowvec com_velocity(unsigned molecule) const
        {
            return com_velocities_.row(molecule);
        }
        inline double position(unsigned atom, unsigned dim) const
        {
            return positions_(atom, dim);
        }
        inline double velocity(unsigned atom, unsigned dim) const
        {
            return velocities_(atom, dim);
        }
        inline double com_position(unsigned molecule, unsigned dim) const
        {
            return com_positions_(molecule, dim);
        }
        inline double com_velocity(unsigned molecule, unsigned dim) const
        {
            return com_velocities_(molecule, dim);
        }
        inline double mass(unsigned atom) const { return index_to_mass_[atom]; }
        inline double molecule_mass(unsigned molecule) const
        {
            return molecular_index_to_mass_[molecule];
        }

    private:
        std::string name_;
        std::vector<std::string> atom_names_;
        std::vector<std::string> molecule_names_;
        std::vector<double> index_to_mass_;
        std::vector<unsigned> index_to_molecule_;
        std::vector<unsigned> index_to_molecular_index_;
        std::vector<double> molecular_index_to_mass_;
        std::vector<unsigned>molecular_index_to_molecule_;
        unsigned num_molecules_;
        mat positions_;
        mat velocities_;
        mat com_positions_;
        mat com_velocities_;
        std::vector<unsigned> indices_;
        unsigned group_size_;
        void read_gro(const std::string& gro,
            const std::vector<Molecule> molecules);

};

#endif
