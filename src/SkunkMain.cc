
#include <vector>
#include <string>

#include "SkunkBase.hh"
#include "libFeathersUtils/Permute.hh"
#include "libFeathersSchemes/SkunkFvSolver.hh"

#include <chrono>
#include <fstream>
#include <iomanip>

inline std::string my_to_string(feathers::uint_t i) {
    std::string s = std::to_string(i);
    std::string z(5-s.size(), '0');
    return z + s;
}

#if 0
static void print(
        feathers::int_t nn,
        const cMesh & m,
        const std::array<real_t, 5>* u) {
    std::cout << nn << std::endl;
    std::ofstream file("out/fields-" + my_to_string(nn) + ".csv");
    file << std::setprecision(std::numeric_limits<real_t>::digits10 + 1);
    file << "x,y,z,r,p,vx,vy,vz" << std::endl;
    for (uint_t cell_ind = m.begin_marked_cell(0); cell_ind != m.end_marked_cell(0); ++cell_ind) {
        const auto& c = m.get_cell(cell_ind);
        MhdHydroVars v({}, u[cell_ind].data());
        auto p = m.get_cell_center_coords(cell_ind);
        file << p.x << ',' << p.y << ',' << p.z-0.5 << ','
             << v.rho << ',' << v.p << ','
             << v.V.x << ',' << v.V.y << ',' << v.V.z << ','
             << std::endl;
    }
}
#endif

template<unsigned N>
static void print_vtk(feathers::uint_t nn,
                      const cMesh& m,
                      const std::array<real_t, N>* u) {
    using namespace feathers;
    std::ofstream file("out/fields-" + my_to_string(nn) + ".vtk");
    file << std::setprecision(std::numeric_limits<real_t>::digits10 + 1);
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "kek" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    file << "POINTS " << m.num_nodes() << " double" << std::endl;
    std::for_each(begin_node(m), end_node(m), [&](tNodeIter node) {
        const vec3_t& coords = node.get_coords();
        file << coords.x << " " << coords.y << " " << coords.z << std::endl;
    });

    file << "CELLS " << m.num_marked_cells(0) << " " << m.num_marked_cells(0)*4 << std::endl;
    std::for_each(begin_interior_cell(m), end_interior_cell(m), [&](tCellIter cell) {
        file << "3 ";
        cell.for_each_node([&](uint_t node_index) {
            file << node_index << " ";
        });
        file << std::endl;
    });
    file << "CELL_TYPES " << m.num_marked_cells(0) << std::endl;
    std::for_each(begin_interior_cell(m), end_interior_cell(m), [&](tCellIter cell) {
        file << "5" << std::endl;
    });

    file << "CELL_DATA " << m.num_marked_cells(0) << std::endl;
    for (uint_t i = 0; i < 1; ++i) {
        const char* const names[]{"rho", "p", "vx", "vy", "vz", "bx", "by", "bz"};
        file << "SCALARS " << names[i] << " double 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        std::for_each(begin_interior_cell(m), end_interior_cell(m), [&](tCellIter cell) {
            MhdHydroVars v({}, u[cell].data());
            file << v.prim[i] << std::endl;
        });
    }
}


#if 1

int main(int argc, char** argv) {
    set_max_num_threads(10);

    std::shared_ptr<cMesh> mesh(new cMesh(2));

#if 1
    mesh->read_triangle("mesh/step_.1.");

    tScalarField<5> uc(mesh->num_cells());
    tScalarField<5> up(mesh->num_cells());
    MhdFvSolverT<MhdPhysicsIdealGas> solver(mesh);
    for (uint_t cell_ind = 0; cell_ind < mesh->num_cells(); ++cell_ind) {
        /*double d = 1.0, p = 1.0, u = 1.0, v = 1.0;
        double x = mesh->get_cell(cell_ind).get_center_coords().x;
        double y = mesh->get_cell(cell_ind).get_center_coords().y;
        d += 0.2*std::sin(M_PI*(x+y));
        std::array<real_t, 5> q{ d, p, u, v, 0.0 };
        MhdPhysicsIdealGas::tFluidState qq({}, nullptr, q.data());
        uc[cell_ind] = qq.cons;*/

        /*double d, p;
        if (mesh->get_cell(cell_ind).get_center_coords().x < 1.0) {
            d = 1.0, p = 1.0;
        } else {
            d = 0.125, p = 0.1;
        }
        std::array<real_t, 5> q{ d, p, 0.0, 0.0, 0.0 };
        MhdHydroVars v({}, nullptr, q.data());
        uc[cell_ind] = v.cons;*/

        /*double d, p, w;
        if (mesh->get_cell(cell_ind).get_center_coords().x < 1.0) {
            d = 3.857134, p = 10.33333, w = 2.629369;
        } else {
            auto x = mesh->get_cell(cell_ind).get_center_coords().x;
            d = 1.0 + 0.2*std::sin(5.0*x);
            p = 1, w = 0;
        }
        std::array<real_t, 5> q{d, p, w, 0.0, 0.0 };
        MhdHydroVars v({}, nullptr, q.data());
        uc[cell_ind] = v.cons;*/

        /*double d, p, u, v;
        double X = mesh->get_cell(cell_ind).get_center_coords().x;
        double Y = mesh->get_cell(cell_ind).get_center_coords().y;
        if (X < 0.5) {
            if (Y < 0.5) {
                //LB
                d = 1.0625, p = 0.4, u = 0, v = 0.2145;
            } else {
                //LT
                d = 2, p = 1, u = 0, v = -0.3;
            }
        } else {
            if (Y < 0.5) {
                //RB
                d = 0.51917, p = 0.4, u = 0, v = -1.1259;
            } else {
                //RT
                d = 1, p = 1, u = 0, v = -0.4;
            }
        }
        std::array<real_t, 5> q{ d, p, u, v, 0.0 };
        MhdPhysicsIdealGas::tFluidState qq({}, nullptr, q.data());
        uc[cell_ind] = qq.cons;*/

        /*std::array<real_t, 5> q{ 1.4, 1.0, 3.0, 0.0, 0.0 };
        MhdHydroVars v({}, nullptr, q.data());
        uc[cell_ind] = v.cons;*/

        std::array<real_t, 5> q{ 1.0, 1.0, 1.0, 0.0, 0.0 };
        MhdHydroVars v({}, nullptr, q.data());
        uc[cell_ind] = v.cons;
    }
    real_t dt = 1e-4;
#endif

    system("rm out/*");
    const uint_t freq = 500;
    print_vtk<5>(0, *mesh, &uc[0]);
    real_t tt = 0.0;
    {
        std::chrono::high_resolution_clock::time_point t0, t1;
        for (uint_t l = 1; l <= 2000000; ++l) {
            t0 = std::chrono::high_resolution_clock::now();
            solver.calc_step(dt, uc, up);
            t1 = std::chrono::high_resolution_clock::now();
            tt += real_t(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) * 1e-9;

            if (l%freq == 0) {
                std::cout << l/freq << "\t" << tt << "\t" << std::endl;
                print_vtk<5>(l/freq, *mesh, &up[0]);
                tt = 0.0;
            }

            uc.swap(up);
        }
    }
    return 0;
}
#endif
