#include <iostream>
#include <iomanip>
#include "mpi.h"
#include "mesh/mesh.h"
#include "distributed_mesh.h"
#include "fluid_solver.h"
#include "mobs_solver.h"
#include "solver.h"

int main(int argc, char** argv) {
    
    MPI_Init(&argc, &argv); // 初始化MPI环境
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // 获取当前进程的rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); // 获取总进程数
    std::string filename = "example13.in";

    CYCASMesh mesh(filename, rank); // 创建网格对象

    mesh.readCells(rank); // 读取cell数据
    mesh.readFaces(rank); // 读取face数据
    mesh.readBoundaries(rank); // 读取边界数据
    
    // mesh.printInfo(); // 打印网格信息
    DistributedMesh newMesh(mesh);
    newMesh.generateLocalDomain();

    std::vector<OutputData> local_result;
    std::vector<I64> local_id;
    Solver solver(rank);    
    if (!newMesh.fluidMetadata.empty()) {
        for (auto iter = newMesh.fluidMetadata.begin(); iter != newMesh.fluidMetadata.end(); ++iter) {
            auto id = iter->first;
            FluidSolver* fluidSolver = new FluidSolver(newMesh, newMesh.fluidMetadata[id], newMesh.p[id], newMesh.ux[id], newMesh.uy[id], newMesh.uz[id]);
            solver.addFluids(fluidSolver);
        }
    } 
    if (!newMesh.solidMetadata.empty()) {
        for (auto iter = newMesh.solidMetadata.begin(); iter != newMesh.solidMetadata.end(); ++iter) {
            auto id = iter->first;
            MobSolver* mobSolver = new MobSolver(newMesh, newMesh.solidMetadata[id], newMesh.mt[id], newMesh.mx[id], newMesh.my[id], newMesh.mz[id]);
            solver.addMobs(mobSolver);
        }        
    }
    solver.solve();
    for (auto &e: solver.fluids) {
        local_result.insert(local_result.end(), e->result.begin(), e->result.end());
        local_id.insert(local_id.end(), e->idArray.begin(), e->idArray.end());
    }
    for (auto &e: solver.mobs) {
        local_result.insert(local_result.end(), e->result.begin(), e->result.end());
        local_id.insert(local_id.end(), e->idArray.begin(), e->idArray.end());
    }    
    int local_count = local_result.size();
    std::vector<int> recv_counts(size);
    MPI_Gather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> displs(size);
    std::vector<I64> all_ids;
    std::vector<OutputData> all_points;    
    if (rank == 0) {
        int total = 0;
        displs[0] = 0;
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i - 1] + recv_counts[i - 1];
        }
        total = displs[size - 1] + recv_counts[size - 1];
        all_ids.resize(total);
        all_points.resize(total);
    }

        // Step 3: Gather actual data
    MPI_Datatype MPI_DataPoint;
    // 定义结构体的每个字段信息
    int count = 2;
    int blocklengths[2] = {8, 1}; // 每个字段的个数
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_INT};
    MPI_Aint displacements[2];

    // 计算每个字段的偏移量
    MPI_Aint base_address;
    OutputData s;
    MPI_Get_address(&s, &base_address);
    MPI_Get_address(&s.ux, &displacements[0]);
    MPI_Get_address(&s.globalid, &displacements[1]);
    displacements[0] -= base_address;
    displacements[1] -= base_address;

    // 创建新数据类型
    MPI_Type_create_struct(count, blocklengths, displacements, types, &MPI_DataPoint);
    MPI_Type_commit(&MPI_DataPoint);    

    MPI_Gatherv(local_id.data(), local_count, MPI_LONG,
                all_ids.data(), recv_counts.data(), displs.data(), MPI_LONG,
                0, MPI_COMM_WORLD);

    MPI_Gatherv(local_result.data(), local_count, MPI_DataPoint,
                all_points.data(), recv_counts.data(), displs.data(), MPI_DataPoint,
                0, MPI_COMM_WORLD);
    MPI_Type_free(&MPI_DataPoint);    
   
    std::map<I64, OutputData> resultMap;
    if (rank == 0) {
        for (int i = 0; i < all_ids.size(); ++i) {
            resultMap[all_ids[i]] = all_points[i];
            // std::cout << all_ids[i] << " " << all_points[i].ux << std::endl;
        }
    }
    // Step 4: Write VTK file on rank 0
    if (rank == 0) {
        std::ofstream fout("result.vtk");
        fout << "# vtk DataFile Version 3.0\n";
        fout << "MPI gathered data with scalar p\n";
        fout << "ASCII\n";
        fout << "DATASET STRUCTURED_GRID\n";
        fout << "DIMENSIONS 50 50 50\n";
        fout << "POINTS " << all_points.size() << " double\n";
        fout << std::fixed << std::setprecision(12);
        for (auto &e: resultMap) {
            fout << e.second.locationx << " " << e.second.locationy << " " << e.second.locationz << "\n";
        }
        // for (size_t i = 0; i < all_points.size(); ++i) {
        //     auto z = i / (mesh.Ny * mesh.Nx);
        //     auto y = (i - z * mesh.Ny * mesh.Nx) / mesh.Nx;
        //     auto x = i - z * mesh.Ny * mesh.Nx - y * mesh.Nx;
        //     fout << x << " " << y << " " << z << "\n";
        // }        
        fout << "POINT_DATA " << all_points.size() << "\n";

        fout << "SCALARS pressure double 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto &e: resultMap) {
            fout << e.second.p << "\n";
        }

        fout << "VECTORS velocity double\n";
        for (auto &e: resultMap) {
            fout << e.second.ux << " " << e.second.uy << " " << e.second.uz << "\n";
        }

        fout << "SCALARS temperature double 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto &e: resultMap) {
            fout << e.second.t << "\n";
        } 
        fout << "SCALARS id int 1\n";
        fout << "LOOKUP_TABLE default\n";
        for (auto &e: resultMap) {
            fout << e.second.globalid << "\n";
        }                
        // for (const auto& pt : all_points) {
        //     fout << pt.ux << " " << pt.uy << " " << pt.uz << "\n";
        // }

        fout.close();
        std::cout << "VTK file written: result.vtk\n";
    }

    MPI_Finalize(); // 结束MPI环境

    return 0;
}

// 1. hypre 单进程 PCG AMG
// 2. 
