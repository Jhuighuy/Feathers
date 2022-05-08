
#include <vector>
#include <string>

#include "SkunkBase.hh"
#include "libFeathersUtils/ChickenThoughts.hh"
#include "libFeathersUtils/Image.hh"
#include "libFeathersSchemes/SkunkFvSolver.hh"

#include <chrono>
#include <fstream>
#include <iomanip>

inline std::string my_to_string(feathers::uint_t i) {
    std::string s = std::to_string(i);
    std::string z(5-s.size(), '0');
    return z + s;
}

#if 1
int main(int argc, char** argv) {
    set_max_num_threads(10);
    print_cockatiel();

    std::shared_ptr<cMesh> mesh(new cMesh(2));

    //cImage2D image;
    //image.init(2048, 2048);
    //image.load("mesh/img/Domain-318.ppm");
    //image.store("mesh/img/Domain-318.jpg");

#if 1
    //mesh->read_image2D("mesh/img/Domain-100-Tube.ppm",
    //                   {{eWhitePixel, 1}, {eRedPixel, 2}, {eGreenPixel, 3}, {eBluePixel, 4}});
    //mesh->save_vtk("mesh/img/Domain-318.vtk", {});
    //mesh->save_strm("mesh/img/Domain-318.strm");
    //return 1;

    mesh->read_triangle("mesh/step_.1.");

    tScalarField uc(5, mesh->NumCells());
    tScalarField up(5, mesh->NumCells());
    MhdFvSolverT<tGasPhysics> solver(mesh);
    for (uint_t cell_ind = 0; cell_ind < mesh->NumCells(); ++cell_ind) {
        /*double d = 1.0, p = 1.0, u = 1.0, v = 1.0;
        double x = mesh->get_cell(cell_ind).CenterPos().x;
        double y = mesh->Cell(cell_ind).CenterPos().y;
        d += 0.2*std::sin(M_PI*(x+y));
        std::array<real_t, 5> q{ d, p, u, v, 0.0 };
        tGasPhysics::tFluidState qq({}, nullptr, q.data());
        uc[cell_ind] = qq.cons;*/

        /*double d, p;
        if (mesh->Cell(cell_ind).CenterPos().x < 1.0) {
            d = 1.0, p = 1.0;
        } else {
            d = 0.125, p = 0.1;
        }
        std::array<real_t, 5> q{ d, p, 0.0, 0.0, 0.0 };
        MhdHydroVars v({}, nullptr, q.data());
        uc[cell_ind] = v.cons;*/

        /*double d, p, w;
        if (mesh->get_cell(cell_ind).CenterPos().x < 1.0) {
            d = 3.857134, p = 10.33333, w = 2.629369;
        } else {
            auto x = mesh->Cell(cell_ind).CenterPos().x;
            d = 1.0 + 0.2*std::sin(5.0*x);
            p = 1, w = 0;
        }
        std::array<real_t, 5> q{d, p, w, 0.0, 0.0 };
        MhdHydroVars v({}, nullptr, q.data());
        uc[cell_ind] = v.cons;*/

        /*double d, p, u, v;
        double X = mesh->get_cell(cell_ind).CenterPos().x;
        double Y = mesh->Cell(cell_ind).CenterPos().y;
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
        tGasPhysics::tFluidState qq({}, nullptr, q.data());
        uc[cell_ind] = qq.cons;*/

        /*std::array<real_t, 5> q{ 1.4, 1.0, 3.0, 0.0, 0.0 };
        MhdHydroVars v({}, nullptr, q.data());
        uc[cell_ind] = v.cons;*/

        std::array<real_t, 5> q{ 1.0, 1.0, 1.0, 0.0, 0.0 };
        MhdHydroVars v({}, nullptr, q.data());
        v.make_cons(5, uc[cell_ind].data());
    }
    real_t dt = 1e-4;
#endif

    system("rm out/*");
    const uint_t freq = 200;
    mesh->save_vtk(("out/fields-" + my_to_string(0) + ".vtk").c_str(),
                   { { "rho", 1, &uc } });
    //print_vtk<5>(0, *mesh, (std::array<real_t, 5>*)uc[0].data());
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
                //print_vtk<5>(l/freq, *mesh, (std::array<real_t, 5>*)up[0].data());
                mesh->save_vtk(("out/fields-" + my_to_string(l/freq) + ".vtk").c_str(),
                               { { "rho", 1, &uc } });
                tt = 0.0;
            }

            uc.swap(up);
        }
    }
    return 0;
}
#endif
