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
#include "z_vec.hpp"
#include "z_conversions.hpp"
#include "z_tcf.hpp"
#include "z_molecule.hpp"
#include "z_atom_group.hpp"
#include "z_gromacs.hpp"
#include "xdrfile_trr.h"
#include "boost/program_options.hpp"

namespace po = boost::program_options;
using namespace arma;
// Units are nm, ps.

int main (int argc, char *argv[])
{
    int st;
    Sim_params params;
    TCF<double> z_dev_tcf(6000);
    double zavg = 0.0;
    unsigned max_steps = std::numeric_limits<unsigned>::max();
    bool surfCheck = false, dev = false;
    enum temp_type_t {MOL_COM, ATOM, NONE};

    po::options_description desc("Options");
    desc.add_options()
        ("help,h",  "Print help messages")
        ("group,g", po::value<std::string>()->default_value("He"),
            "Group for density/temperature profiles")
        ("liquid,l", po::value<std::string>()->default_value("OW"),
            "Group to use for calculation of surface")
        ("index,n", po::value<std::string>()->default_value("index.ndx"),
            ".ndx file containing atomic indices for groups")
        ("gro", po::value<std::string>()->default_value("conf.gro"),
            ".gro file containing list of atoms/molecules")
        ("top", po::value<std::string>()->default_value("topol.top"),
            ".top file containing atomic/molecular properties")
        ("temp_type,T", po::value<std::string>()->default_value("none"),
            "Use atom or mol_com velocities for temperature profile")
        ("max_time,t", po::value<double>()->default_value(
            std::numeric_limits<double>::max()),
            "Maximum simulation time to use in calculations")
        ("avg,a", "Calculate average Z position of surface")
        ("density,d", "Calculate density profile for selected group");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        exit(EXIT_SUCCESS);
    }

    temp_type_t temp_type;
    const std::string& vm_temp_type = vm["temp_type"].as<std::string>();
    if (vm_temp_type == "m" ||
        vm_temp_type == "mol_com")
        temp_type = MOL_COM;
    else if (vm_temp_type == "a" ||
        vm_temp_type == "atom")
        temp_type = ATOM;
    else if (vm_temp_type == "n" ||
        vm_temp_type == "none" )
        temp_type = NONE;
    else
    {
        std::cerr << std::endl << "Unrecognized temp_type option: " << vm_temp_type << std::endl;
        exit(EXIT_FAILURE);
    }

    bool density = vm.count("density") ? true : false;
    bool avg = vm.count("avg") ? true : false;
    std::map<std::string, std::vector<unsigned> > groups;
    groups = read_ndx(vm["index"].as<std::string>());

    std::vector<Molecule> molecules = gen_molecules(vm["top"].as<std::string>(), params);
    Atom_group all_atoms(vm["gro"].as<std::string>(), molecules);

    Atom_group selected_group;
    selected_group.set_indices(vm["group"].as<std::string>(),
        select_group(groups, vm["group"].as<std::string>()),
        all_atoms);
    Atom_group liquid_group;
    liquid_group.set_indices(vm["liquid"].as<std::string>(),
        select_group(groups, vm["liquid"].as<std::string>()),
        all_atoms);

    rvec *x_in = NULL;
    matrix box_mat;
    rowvec box = zeros<rowvec>(DIMS);
    std::string xtc_filename = "prod.xtc";
    std::string trr_filename = "prod.trr";
    XDRFILE *xtc_file, *trr_file;
    params.extract_traj_metadata(strdup(xtc_filename.c_str()), (&x_in), box);
    trr_file = xdrfile_open(strdup(trr_filename.c_str()), "r");
    xtc_file = xdrfile_open(strdup(xtc_filename.c_str()), "r");
    Grid grid(box);
    params.boxLength = box[0];
    params.zboxLength = box[2];
    params.set_max_time(vm["max_time"].as<double>());

    std::string surf_filename = "surf.dat";
    std::fstream surf_file;
    surf_file.open(surf_filename.c_str(),  std::fstream::in);
    if (!surf_file.is_open())
        surf_file.open(surf_filename.c_str(),  std::fstream::out);
    else
    {
        surfCheck = true;
    }
    double cgLength = 0.24, cgLength2 = cgLength*cgLength;
    double rmax = 3.0*cgLength, rmax2 = rmax*rmax;
    double rho_pre_factor = 1.0/(2.0*M_PI*cgLength2);
    rho_pre_factor = sqrt(rho_pre_factor*rho_pre_factor*rho_pre_factor);
    double rho_exp_factor = -1.0/2.0/cgLength2;
    double phi_rmax = rho_pre_factor*exp(-9.0/2.0);
    unsigned zcount = 0;

    cube zdev = zeros(grid.size_x(),grid.size_y(), z_dev_tcf.length());

    int densSize, velSize = 16;
    double densSpacing = 0.01;
    rowvec dx, dx2;
    double vMax = 4.0;

    densSize = 2.0*static_cast<unsigned>(params.boxLength/densSpacing/2);
    irowvec counter(densSize);
    rowvec zdensity(densSize);
    rowvec tempxyprofile(densSize);
    rowvec tempzprofile(densSize);
    mat zvdensity(densSize,velSize);
    rvec *v_in = NULL;
    v_in = new rvec [params.num_atoms()];
    float time, lambda, prec;
    unsigned steps = 0;
    for (unsigned step=0, steps=0; step<max_steps; step++)
    {
        steps++;
        unsigned mod = step%z_dev_tcf.length();
        unsigned start = (step+1)%z_dev_tcf.length();

        if(read_xtc(xtc_file, params.num_atoms(), &st, &time, box_mat, x_in,
                    &prec))
            break;
        if(read_trr(trr_file, params.num_atoms(), &st, &time, &lambda, box_mat,
                    NULL, v_in, NULL))
            break;
        if (selected_group.num_molecules()) selected_group.zero_com();
        unsigned i = 0;
        for (std::vector<unsigned>::iterator i_atom = selected_group.begin();
            i_atom != selected_group.end(); i_atom++, i++)
        {
            selected_group.set_position(i, x_in[*i_atom]);
            selected_group.set_velocity(i, v_in[*i_atom]);
        }
        selected_group.update_com();
        i = 0;
        for (std::vector<unsigned>::iterator i_atom = liquid_group.begin();
            i_atom != liquid_group.end(); i_atom++, i++)
        {
            liquid_group.set_position(i, x_in[*i_atom]);
        }
        if (surfCheck)
        {
            grid.read_mesh(surf_file);

            if (avg) zavg += grid.avg_mesh();
            if (dev)
                for (unsigned i_x=0; i_x<grid.size_x(); i_x++)
                    for (unsigned i_y=0; i_y<grid.size_y(); i_y++)
                    {
                        zdev(i_x,i_y,mod) = (grid.upper_mesh(i_x,i_y,2) - zavg);
                        if (step>=z_dev_tcf.length()-1)
                        {
                            z_dev_tcf.correlate_1D_one_direction(
                                zdev.tube(i_x,i_y), start);
                            //Previously was
                            //correlation1D_mod(zdevCorr, zdev.tube(i_x,i_y),
                            //    zdev.tube(i_x,i_y), start, corr.length);
                            zcount++; // Should include count in TCF.
                        }
                    }
        }
        else
        {
            grid.zero_rho();
            unsigned grid_cutoff = rmax/grid.spacing + 1;
            for (unsigned i_atom=0; i_atom<liquid_group.size(); i_atom++)
            {
                irowvec closest_grid_point(DIMS);
                closest_grid_point = conv_to<irowvec>::from(
                    floor(liquid_group.position(i_atom)/grid.spacing)+
                    grid.size());

                irowvec min_grid_point = closest_grid_point - grid_cutoff;
                irowvec max_grid_point = closest_grid_point + grid_cutoff;
                for (int i_x=min_grid_point(0); i_x<max_grid_point(0); i_x++)
                {
                    int i_x_shift = i_x%grid.size_x();
                    for (int i_y=min_grid_point(1); i_y<max_grid_point(1);
                        i_y++)
                    {
                        int i_y_shift = i_y%grid.size_y();
                        for (int i_z=min_grid_point(2); i_z<max_grid_point(2);
                            i_z++)
                        {
                            int i_z_shift = i_z%grid.size_z();
                            rowvec grid_point = veclike_to_row(
                                grid.grid[i_x_shift][i_y_shift][i_z_shift],
                                DIMS);
                            findDx_noShift(dx, grid_point,
                                liquid_group.position(i_atom), box);
                            double r2 = dot(dx,dx);
                            if (r2>rmax2) continue;
                            grid.rho(i_x_shift,i_y_shift,i_z_shift) +=
                                rho_pre_factor*exp(r2*rho_exp_factor) -
                                phi_rmax;
                        }
                    }
                }
            }

            for (unsigned i_x=0; i_x<grid.size_x(); i_x++)
                for (unsigned i_y=0; i_y<grid.size_y(); i_y++)
                {
                    bool FOUND_UPPER = false, FOUND_LOWER= false;
                    for (unsigned i_z=0; i_z<grid.size_z(); i_z++)
                    {
                        double rhoCutoff = 16.0;// 14.75;
                        int i_z_shift = (i_z == 0) ? (grid.size_z()-1) : (i_z-1);
                        if (!FOUND_UPPER && grid.rho(i_x,i_y,i_z) < rhoCutoff &&
                            grid.rho(i_x,i_y,i_z_shift) > rhoCutoff)
                        {
                            grid.upper_mesh(i_x,i_y,2) =
                                grid.grid[i_x][i_y][i_z][2] -
                                (rhoCutoff-grid.rho(i_x,i_y,i_z))/
                                (grid.rho(i_x,i_y,i_z_shift) -
                                grid.rho(i_x,i_y,i_z))*grid.spacing;
                            if (i_z == 0) grid.lower_mesh(i_x,i_y,2) += box[2];
                            FOUND_UPPER = true;
                        }
                        if (!FOUND_LOWER && grid.rho(i_x,i_y,i_z) > rhoCutoff &&
                            grid.rho(i_x,i_y,i_z_shift) < rhoCutoff)
                        {
                            grid.lower_mesh(i_x,i_y,2) =
                                grid.grid[i_x][i_y][i_z][2] -
                                (grid.rho(i_x,i_y,i_z)-rhoCutoff)/
                                (grid.rho(i_x,i_y,i_z) -
                                grid.rho(i_x,i_y,i_z_shift))*grid.spacing;
                            if (grid.lower_mesh(i_x,i_y,2) < 0.0)
                                grid.lower_mesh(i_x,i_y,2) += box[2];
                            FOUND_LOWER = true;
                        }
                        if (FOUND_UPPER && FOUND_LOWER) break;
                    }
                    surf_file << grid.upper_mesh(i_x,i_y,2) << " " <<
                        grid.lower_mesh(i_x,i_y,2) << std::endl;
                }
        }
        if (temp_type==ATOM)
            for (unsigned i_atom=0; i_atom<selected_group.size(); i_atom++)
            {
                unsigned min_x = 0, min_y = 0;
                double rmin = 1000.0;
                for (unsigned i_x=0; i_x<grid.size_x(); i_x++)
                {
                    for (unsigned i_y=0; i_y<grid.size_y(); i_y++)
                    {
                        findDx_noShift(dx,
                            vec(grid.upper_mesh.tube(i_x,i_y)).t(),
                            selected_group.position(i_atom), box);
                        double r2 = dot(dx,dx);
                        if (r2<rmin)
                        {
                            min_x = i_x;
                            min_y = i_y;
                            rmin = r2;
                        }
                    }
                }
                int min_x_shiftl = (min_x == 0) ? (grid.size_x()-1) : (min_x-1);
                int min_y_shiftl = (min_y == 0) ? (grid.size_y()-1) : (min_y-1);
                int min_x_shiftr = (min_x == grid.size_x()-1) ? 0 : (min_x+1);
                int min_y_shiftr = (min_y == grid.size_y()-1) ? 0 : (min_y+1);
                findDx_noShift(dx, grid.upper_mesh.tube(min_x_shiftr,min_y),
                    grid.upper_mesh.tube(min_x_shiftl,min_y), box);
                findDx_noShift(dx2, grid.upper_mesh.tube(min_x,min_y_shiftr),
                    grid.upper_mesh.tube(min_x,min_y_shiftl), box);
                rowvec norm = normalise(cross(dx,dx2));
                findDx_noShift(dx, selected_group.position(i_atom),
                    vec(grid.upper_mesh.tube(min_x,min_y)).t(), box);
                double r2 = dot(dx,norm);
                if (abs(r2) > params.boxLength) continue;
                int which_bin =
                    static_cast<unsigned>(floor(r2/params.boxLength*densSize/2)) +
                    densSize/2;
                if (temp_type==ATOM)
                {
                    double mass = selected_group.mass(i_atom);
                    if (which_bin < 0 || which_bin > densSize) continue;
                    tempxyprofile[which_bin] += mass*
                        (selected_group.velocity(i_atom,0)*
                        selected_group.velocity(i_atom,0) +
                        selected_group.velocity(i_atom,1)*
                        selected_group.velocity(i_atom,1));
                    tempzprofile[which_bin] +=
                        mass*selected_group.velocity(i_atom,2)*
                        selected_group.velocity(i_atom,2);
                    counter[which_bin] += 1;
                }
                zdensity[which_bin] += 1.0;
                //zvdensity[(int)floor(r2/params.boxLength*densSize/2)+
                    //densSize/2][(int)(sqrt(v2)/vMax*velSize)] += 1.0;
            }
        if (temp_type == MOL_COM)
            for (unsigned i_mol=0; i_mol<selected_group.num_molecules(); i_mol++)
            {
                selected_group.com_position(i_mol) /=
                    selected_group.molecule_mass(i_mol);
                selected_group.com_velocity(i_mol) /=
                    selected_group.molecule_mass(i_mol);
                unsigned min_x = 0, min_y = 0;
                double rmin = 1000.0;
                for (unsigned i_x=0; i_x<grid.size_x(); i_x++)
                {
                    for (unsigned i_y=0; i_y<grid.size_y(); i_y++)
                    {
                        findDx_noShift(dx,
                            vec(grid.upper_mesh.tube(i_x,i_y)).t(),
                            selected_group.com_position(i_mol), box);
                        double r2 = dot(dx,dx);
                        if (r2<rmin)
                        {
                            min_x = i_x;
                            min_y = i_y;
                            rmin = r2;
                        }
                    }
                }
                unsigned min_x_shiftl = (min_x == 0) ?
                    (grid.size_x()-1) : (min_x-1);
                unsigned min_y_shiftl = (min_y == 0) ?
                    (grid.size_y()-1) : (min_y-1);
                unsigned min_x_shiftr = (min_x == grid.size_x()-1) ?
                    0 : (min_x+1);
                unsigned min_y_shiftr = (min_y == grid.size_y()-1) ?
                    0 : (min_y+1);
                findDx_noShift(dx, grid.upper_mesh.tube(min_x_shiftr,min_y),
                    grid.upper_mesh.tube(min_x_shiftl,min_y), box);
                findDx_noShift(dx2, grid.upper_mesh.tube(min_x,min_y_shiftr),
                    grid.upper_mesh.tube(min_x,min_y_shiftl), box);
                rowvec norm = normalise(cross(dx,dx2));
                findDx_noShift(dx, selected_group.com_position(i_mol),
                    vec(grid.upper_mesh.tube(min_x,min_y)).t(), box);
                double r2 = dot(dx,norm);
                if (abs(r2) > params.boxLength) continue;
                zdensity[static_cast<unsigned>(floor(r2/params.boxLength*densSize/2)) +
                    densSize/2] += 1.0;
                //zvdensitystatic_cast<unsigned>(floor(r2/params.boxLength*densSize/2)) +
                //    densSize/2][(int)(sqrt(v2)/vMax*velSize)] += 1.0;

                if (temp_type==MOL_COM)
                {
                    int which_bin =
                        static_cast<unsigned>(floor(r2/params.boxLength*densSize/2)) +
                            densSize/2;
                    if (which_bin < 0 || which_bin > densSize) continue;
                    // Don't need to multiply by mass since com velocity above
                    // didn't divide by mass
                    tempxyprofile[which_bin] +=
                        selected_group.molecule_mass(i_mol)*
                        dot(selected_group.com_velocities()(i_mol, span(0,1)),
                            selected_group.com_velocities()(i_mol, span(0,1)));
                    tempzprofile[which_bin] +=
                        selected_group.molecule_mass(i_mol)*
                        selected_group.com_velocity(i_mol,2)*
                        selected_group.com_velocity(i_mol,2);
                    counter[which_bin] += 1;
                }
            }
    }
    xdrfile_close(xtc_file);
    xdrfile_close(trr_file);
    surf_file.close();

    if (avg) std::cout << "Avg is " << zavg/steps << std::endl;

    if (dev)
    {
        std::ofstream dev_file;
        dev_file.open("zdev.txt");
        for (unsigned i=0; i<z_dev_tcf.length(); i++)
            dev_file << i*params.dt() << " " << z_dev_tcf.tcf(i)/zcount/2 <<
            std::endl;
        dev_file.close();
    }


    if (density)
    {
        std::ofstream density_file;
        density_file.open((vm["group"].as<std::string>() + "_density.txt").c_str());
        for (int i=0; i<densSize; i++)
            density_file << (i-densSize/2)*params.boxLength/(densSize/2) <<
            " " << zdensity[i]/steps/params.boxLength/params.boxLength <<
            std::endl;
        density_file.close();
    }

    //for (int i=0; i<velSize; i++)
    //    zvdensity[densSize/2][i] /=2.0;
    for (int i_d=0; i_d<densSize; i_d++)
    {
        //double sum = 0.0;
        //for (int i_v=0; i_v<velSize; i_v++)
        //{
        //    sum += zvdensity[i_d][i_v];
        //}
        //for (int i_v=0; i_v<velSize; i_v++)
        //{
        //    if (sum>0.0) zvdensity[i_d][i_v]/=sum;
        //}
        if (counter(i_d))
        {
            tempxyprofile[i_d] *= AMU_TO_KG*1000*1000/counter(i_d)/2.0/KB;
            tempzprofile[i_d] *= AMU_TO_KG*1000*1000/counter(i_d)/KB;
        }
    }

    std::string xy_temperature_filename;
    std::ofstream xy_temperature_file;
    if(temp_type==MOL_COM)
        xy_temperature_filename = vm["group"].as<std::string>() + "_tempxyinstprofile_com.txt";
    else if(temp_type==ATOM)
        xy_temperature_filename = vm["group"].as<std::string>() + "_tempxyinstprofile_atom.txt";
    xy_temperature_file.open(xy_temperature_filename.c_str());
    for (int i=0; i<densSize; i++)
        xy_temperature_file <<
            (i-densSize/2.0)*params.boxLength/(densSize/2.0) << " " <<
            tempxyprofile[i] << std::endl;
    xy_temperature_file.close();

    std::string z_temperature_filename;
     std::ofstream z_temperature_file;
    if(temp_type==MOL_COM)
        z_temperature_filename = vm["group"].as<std::string>() + "_tempzinstprofile_com.txt";
    else if(temp_type==ATOM)
        z_temperature_filename = vm["group"].as<std::string>() + "_tempzinstprofile_atom.txt";
    z_temperature_file.open(z_temperature_filename.c_str());
    for (int i=0; i<densSize; i++)
        z_temperature_file <<
            (i-densSize/2.0)*params.boxLength/(densSize/2.0) << " " <<
            tempzprofile[i] << std::endl;
    z_temperature_file.close();
}
