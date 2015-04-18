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

#include <vector>
#include <fstream>
#include "z_atom.hpp"

#ifndef _Z_MOLECULE_HPP_
#define _Z_MOLECULE_HPP_

class Molecule
{
    public:
        Molecule(): num_atoms_(0) {}
        inline std::vector<Atom> atoms() const { return atoms_; }
        inline unsigned  num_atoms() const { return num_atoms_; }
        inline std::string name() const { return name_; }
        inline void set_name(const std::string& name) { name_ = name; }
        void add_atom(Atom atom);
        friend std::ostream &operator<<( std::ostream &output,
            Molecule &molecule );
        void print() const;

    private:
        unsigned num_atoms_;
        std::vector<Atom> atoms_;
        std::string name_;
};

#endif
