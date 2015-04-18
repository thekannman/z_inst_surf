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

#include <iomanip>
#include "z_atom.hpp"
#include <iostream>

void Atom::print() const
{
    std::cout << std::left << std::setw(9) << "Name ";
    std::cout << std::left << std::setw(9) << "Type ";
    std::cout << std::left << std::setw(9) << "Mass ";
    std::cout << std::left << std::setw(9) << "Charge " << std::endl;
    std::cout << (*this);
}

std::ostream &operator<<( std::ostream &output,
    const Atom &atom )
{
    output << std::left << std::setw(9) << atom.name_;
    output << std::left << std::setw(9) << atom.type_;
    output << std::left << std::setw(9) << std::setprecision(6) <<
        atom.mass_;
    output << std::left << std::setw(9) << std::setprecision(4) <<
        atom.charge_ << std::endl;
    return output;
}
