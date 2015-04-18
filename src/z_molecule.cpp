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

#include "z_molecule.hpp"
#include <iomanip>
#include <iostream>

void Molecule::add_atom(Atom atom)
{
    atoms_.push_back(atom);
    num_atoms_++;
}
std::ostream &operator<<( std::ostream &output, Molecule &molecule )
{
    output << molecule.name_ << std::endl;
    output << std::setw(5) << "";
    output << std::left << std::setw(9) << "Name ";
    output << std::left << std::setw(9) << "Type ";
    output << std::left << std::setw(9) << "Mass ";
    output << std::left << std::setw(9) << "Charge " << std::endl;
    for (std::vector<Atom>::iterator i_atom = molecule.atoms_.begin();
        i_atom != molecule.atoms_.end(); i_atom++)
    {
        output << std::setw(5) << "";
        output << (*i_atom);
    }
    return output;
}

void Molecule::print() const
{
    std::cout << name_ << std::endl;
    std::cout << std::setw(5) << "";
    std::cout << std::left << std::setw(9) << "Name ";
    std::cout << std::left << std::setw(9) << "Type ";
    std::cout << std::left << std::setw(9) << "Mass ";
    std::cout << std::left << std::setw(9) << "Charge " << std::endl;
    for (std::vector<Atom>::const_iterator i_atom = atoms_.begin();
        i_atom != atoms_.end(); i_atom++)
    {
        std::cout << std::setw(5) << "";
        std::cout << (*i_atom);
    }
}
