
#include <vector>
#include <string>

#include "SkunkBase.hh"
#include "libFeathersSchemes/SkunkFvSolver.hh"

#include <chrono>
#include <fstream>
#include <iomanip>

inline std::string my_to_string(feathers::uint_t i) {
    std::string s = std::to_string(i);
    std::string z(5-s.size(), '0');
    return z + s;
}

static void print(
        feathers::int_t nn,
        const UMesh & m,
        const std::array<real_t, 5>* u) {
    std::cout << nn << std::endl;
    std::ofstream file("out/fields-" + my_to_string(nn) + ".csv");
    file << std::setprecision(std::numeric_limits<real_t>::digits10 + 1);
    file << "x,y,z,r,p,vx,vy,vz" << std::endl;
    for (uint_t cell_ind1 = 0; cell_ind1 < m.num_marked_cells(0); ++cell_ind1) {
        const uint_t cell_ind = m.get_marked_cell_index(cell_ind1, 0);
        const auto& c = m.get_cell(cell_ind);
        MhdHydroVars v({}, u[cell_ind].data());
        auto p = m.get_cell_center_position(cell_ind);
        file << p.x << ',' << p.y << ',' << p.z-0.5 << ','
             << v.rho << ',' << v.p << ','
             << v.V.x << ',' << v.V.y << ',' << v.V.z << ','
             << std::endl;
    }
}

template<unsigned N>
static void print_vtk(feathers::uint_t nn,
                      const UMesh & m,
                      const std::array<real_t, N>* u) {
    using namespace feathers;
    std::ofstream file("out/fields-" + my_to_string(nn) + ".vtk");
    file << std::setprecision(std::numeric_limits<real_t>::digits10 + 1);
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "kek" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    file << "POINTS " << m.num_nodes() << " double" << std::endl;
    for (uint_t node_ind = 0; node_ind < m.num_nodes(); ++node_ind) {
        const Node& node = m.get_node(node_ind);
        const vec3_t& p = m.get_node_position(node_ind);
        file << p.x << " " << p.y << " " << p.z << std::endl;
    }

    file << "CELLS " << m.num_marked_cells(0) << " " << m.num_marked_cells(0)*4 << std::endl;
    for (uint_t cell_ind = m.begin_cell(0); cell_ind < m.end_cell(0); ++cell_ind) {
        const Cell& cell = m.get_cell(cell_ind);
        file << "3 ";
        std::for_each(cell.begin_node(), cell.end_node(), [&](int_t node_ind) {
            file << node_ind << " ";
        });
        file << std::endl;
    }
    file << "CELL_TYPES " << m.num_marked_cells(0) << std::endl;
    for (uint_t cell_ind = m.begin_cell(0); cell_ind < m.end_cell(0); ++cell_ind) {
        file << "5" << std::endl;
    }

    file << "CELL_DATA " << m.num_marked_cells(0) << std::endl;
    for (uint_t i = 0; i < 1; ++i) {
        const char* const names[]{"rho", "p", "vx", "vy", "vz", "bx", "by", "bz"};
        file << "SCALARS " << names[i] << " double 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        for (uint_t cell_ind = m.begin_cell(0); cell_ind < m.end_cell(0); ++cell_ind) {
            MhdHydroVars v({}, u[cell_ind].data());
            //file << m.get_cell_volume(cell_ind) << std::endl;
            file << v.prim/*u[cell_ind]*/[i] << std::endl;
        }
    }
}

#if 0
int main() {
    using namespace feathers;
    std::shared_ptr<uMesh> mesh(new feathers::structured_mesh_t(0.5, 1.0, 50,
                                                               -SKUNK_PI/125, +SKUNK_PI/125, 1,
                                                               0.0, SKUNK_PI_2, 200,
                                                               1, 2, 3, 3, 5, 4));
    //std::shared_ptr<uMesh> mesh(new feathers::structured_mesh_t(0.1, 1.0, 50,
    //                                                           -SKUNK_PI/2 +SKUNK_PI/2, 200,
    //                                                           0.0, 1.0, 1,
    //                                                           1, 2, 3, 3, 5, 4));
    /*std::shared_ptr<uMesh> mesh(new feathers::structured_mesh_t(0.0, 10.0, 200,
                                                               0.0, 2.50, 50,
                                                               -0.5, 0.5, 1,
                                                               1, 2, 3, 3, 5, 4));*/

    std::ofstream file1("../results/_nodes.csv");
    file1 << "x,y,z" << std::endl;
    for (uint_t node_ind = 0; node_ind != mesh->num_nodes(); ++node_ind) {
        const Node& face = mesh->get_node(node_ind);
        file1
            << face.get_position().x << ","
            << face.get_position().y << ","
            << face.get_position().z << std::endl;
    }
    file1.close();
    std::ofstream file("../results/_normals.csv");
    file << "x,y,z,nx,ny,nz,dx,dy,dz" << std::endl;
    for (uint_t face_ind = 0; face_ind != mesh->num_faces(); ++face_ind) {
        const Face& face = mesh->get_face(face_ind);
        vec3_t direction = mesh->get_cell_center_position(face.get_outer_cell()) -
                           mesh->get_cell_center_position(face.get_inner_cell());
        direction *= safe_inv(direction.len());
        file
                << face.get_center_position().x << ","
                << face.get_center_position().y << ","
                << face.get_center_position().z << ","
                << face.get_normal().x << ","
                << face.get_normal().y << ","
                << face.get_normal().z << ","
                << direction.x << "," << direction.y << "," << direction.z << std::endl;
    }
    file.close();

    tScalarField<5> uc(mesh->num_cells());
    tScalarField<5> up(mesh->num_cells());
    for (uint_t cell_ind = mesh->begin_cell(0); cell_ind < mesh->end_cell(0); ++cell_ind) {
        /*vec3_t n;
        auto& cell = mesh->get_cell(cell_ind);
        std::for_each(cell.begin_face(), cell.end_face(), [&](unsigned f){
            if (f != npos)
            {
                if (mesh->get_face(f).get_outer_cell() == cell_ind)
                    n += mesh->get_face_normal(f)*mesh->get_face_area(f);
                else
                    n -= mesh->get_face_normal(f)*mesh->get_face_area(f);
            }
        });
        n /= cell.get_volume();
        std::array<real_t, 5> q{ n.len(), 1.0, 0.0, 0.0, 0 };
        uc[cell_ind] = q;*/

        auto r = mesh->get_cell(cell_ind).get_center_position();
        auto v = 0.0*mhd_vec3_t({0.0, 0.0, 0.5}).cross(r);
        auto p = 0.0*Gamma*(r.x*r.x + r.y*r.y) + 1.0;
        auto d = Gamma;
        //auto v = mhd_vec3_t({0.0, 0.0, 0.0});//.cross(mesh->get_cell(cell_ind).get_center_position());

        //if (mesh->get_cell(cell_ind).get_center_position().x < 5.0) {
        //    d = 2.0, p = 5.0;
        //} else {
        //    d = 1.0, p = 1.0;
        //}

        std::array<real_t, 5> q{ d, p, v.x, v.y, v.z };
        MhdPhysicsIdealGas::tFluidState qq({}, nullptr, q.data());
        uc[cell_ind] = qq.cons;
    }

    auto dt = 1e-3;
    MhdFvSolverT<MhdPhysicsIdealGas> solver(mesh);

    print(0, *mesh, &uc[0]);
    real_t tt = 0.0;
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::nanoseconds nanoseconds;
    Clock::time_point t0, t1;
    {
        for (int_t l = 1; /*l <= 40000000*/; ++l) {

            t0 = Clock::now();
            solver.calc_step(dt, uc, up);
            t1 = Clock::now();

            tt += std::chrono::duration_cast<nanoseconds>(t1 - t0).count() * 1e-9;
            if (l % 100 == 0) {
                std::cout << "\t" << tt << "\t" << l << "\t\t" << l*dt << std::endl;
                print(l/100, *mesh, &up[0]);
                tt = 0.0;
            }

            uc.swap(up);
        }
    }
}
#endif

#if 1

#include <omp.h>

int main(int argc, char** argv) {
    omp_set_num_threads(20);

    std::shared_ptr<UMesh> mesh(new UMesh(2));

#if 1
    mesh->read_triangle("mesh/step_.1.");

    tScalarField<5> uc(mesh->num_cells());
    tScalarField<5> up(mesh->num_cells());
    MhdFvSolverT<MhdPhysicsIdealGas> solver(mesh);
    for (uint_t cell_ind = 0; cell_ind < mesh->num_cells(); ++cell_ind) {
        /*double d = 1.0, p = 1.0, u = 1.0, v = 1.0;
        double x = mesh->get_cell(cell_ind).get_center_position().x;
        double y = mesh->get_cell(cell_ind).get_center_position().y;
        d += 0.2*std::sin(M_PI*(x+y));
        std::array<real_t, 5> q{ d, p, u, v, 0.0 };
        MhdPhysicsIdealGas::tFluidState qq({}, nullptr, q.data());
        uc[cell_ind] = qq.cons;*/

        /*double d, p;
        if (mesh->get_cell(cell_ind).get_center_position().x < 1.0) {
            d = 1.0, p = 1.0;
        } else {
            d = 0.125, p = 0.1;
        }
        std::array<real_t, 5> q{ d, p, 0.0, 0.0, 0.0 };
        MhdHydroVars v({}, nullptr, q.data());
        uc[cell_ind] = v.cons;*/

        /*double d, p, w;
        if (mesh->get_cell(cell_ind).get_center_position().x < 1.0) {
            d = 3.857134, p = 10.33333, w = 2.629369;
        } else {
            auto x = mesh->get_cell(cell_ind).get_center_position().x;
            d = 1.0 + 0.2*std::sin(5.0*x);
            p = 1, w = 0;
        }
        std::array<real_t, 5> q{d, p, w, 0.0, 0.0 };
        MhdHydroVars v({}, nullptr, q.data());
        uc[cell_ind] = v.cons;*/

        /*double d, p, u, v;
        double X = mesh->get_cell(cell_ind).get_center_position().x;
        double Y = mesh->get_cell(cell_ind).get_center_position().y;
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
    //print(0, *mesh, &uc[0]);
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
