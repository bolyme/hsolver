#ifndef SOLVER
#define SOLVER

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include "mpi.h"
#include "distributed_mesh.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "fluid_solver.h"
#include "mobs_solver.h"

#define Re ((U0*L)/NU)

struct Solver {
    std::vector<FluidSolver*> fluids;
    std::vector<MobSolver*> mobs;
    int id;

    Solver(int id_): id(id_) {}
    void addFluids(FluidSolver*);
    void addMobs(MobSolver*);
    void solve();
    void getResult();
};

void Solver::addFluids(FluidSolver* fluid) {
    fluids.push_back(fluid);
}

void Solver::addMobs(MobSolver* mob) {
    mobs.push_back(mob);
}

void Solver::solve() {
    for (auto &e: fluids) {
        e->init(); 
    }    
    for (auto &e: mobs) {
        e->init(); 
    }     
    clock_t start;
    for (int step = 0; step < 100; ++step) {
	    if (step == 10)
		   start = clock();
        for (auto &e: fluids) {
            e->solve(); 
        }
        for (auto &e: mobs) {
            e->solve();
        }

        int requestCount = 0;
        for (auto &e: fluids) {
            for (auto &v: e->sendData) {
                v.second.resize(0);
            }
        }
        for (auto &e: fluids) {        
            for (auto &v: e->recvData) {
                v.second.resize(0);
            }   
        }  
        for (auto &e: mobs) {
            for (auto &v: e->sendData) {
                v.second.resize(0);
            }
        }
        for (auto &e: mobs) {        
            for (auto &v: e->recvData) {
                v.second.resize(0);
            }   
        }             
        for (auto &e: fluids) {
            for (auto &v: e->t_boundary) {
                requestCount++;
                // std::cout << "send fluid:";
                for (auto &a: v.second) {
                    e->sendData[v.first].push_back(e->t_temp[a]);
                    // std::cout << e->t_temp[a] << " " ;
                }  
                e->recvData[v.first].resize(e->t_halo[v.first].size());     
            }
        }
        for (auto &e: mobs) {
            for (auto &v: e->t_boundary) {
                requestCount++;
                for (auto &a: v.second) {
                    e->sendData[v.first].push_back(e->t_temp[a]);
                    
                }
                e->recvData[v.first].resize(e->sendData[v.first].size()); 
            }
        }
        MPI_Request request[requestCount*2];
        int iter = 0;
        for (auto &e: fluids) {           
            for (auto &v: e->t_boundary) {
                // std::cout << "fluid send: " << e->sendData[v.first].size() << std::endl;
                MPI_Isend(e->sendData[v.first].data(), e->sendData[v.first].size(), MPI_DOUBLE, v.first, 8,
                        MPI_COMM_WORLD, &request[iter++]);                           
            }
        }
        for (auto &e: mobs) {           
            for (auto &v: e->t_boundary) {
                MPI_Isend(e->sendData[v.first].data(), e->sendData[v.first].size(), MPI_DOUBLE, v.first, 17,
                        MPI_COMM_WORLD, &request[iter++]);                           
            }
        } 
        for (auto &e: fluids) {           
            for (auto &v: e->t_boundary) {
                MPI_Irecv(e->recvData[v.first].data(), e->recvData[v.first].size(), MPI_DOUBLE, v.first, 17,
                        MPI_COMM_WORLD, &request[iter++]);                           
            }
        }               
        for (auto &e: mobs) {           
            for (auto &v: e->t_boundary) {
                // std::cout << "mobs recv: " << e->recvData[v.first].size() << std::endl;
                MPI_Irecv(e->recvData[v.first].data(), e->recvData[v.first].size(), MPI_DOUBLE, v.first, 8,
                        MPI_COMM_WORLD, &request[iter++]);                           
            }
        }
        // std::cout << "iter: " << iter << " preCount: " << requestCount * 2 << std::endl;
        MPI_Waitall(iter, request, MPI_STATUSES_IGNORE);             
        for (auto &e: fluids) {                                  
            for (auto &v: e->t_halo) {
                for (int a = 0; a < v.second.size(); ++a) {
                    // std::cout << e->t_halo[v.first].size() << " " << e->recvData[v.first].size() << std::endl;
                    e->t_temp[e->t_halo[v.first][a]] = e->recvData[v.first][a];
                }
            } 
        }  
        for (auto &e: mobs) {                                  
            for (auto &v: e->t_halo) {
                for (int a = 0; a < v.second.size(); ++a) {
                    e->t_temp[e->t_halo[v.first][a]] = e->recvData[v.first][a];
                }
            } 
        }
        for (auto &e: fluids) {
            e->setForNext(); 
        }    
        for (auto &e: mobs) {
            e->setForNext(); 
        }                        
    }
    clock_t end = clock();
        double duration = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "耗时: " << duration << " 秒 " <<  " average: " << duration / 90.0  << std::endl;

    I64 gi, gj, gk;
    for (auto &e: fluids) {
        auto iter = e->mesh.cycasMesh.fluidRegions.begin();
        auto cellIdxArray = e->mesh.cycasMesh.fluidRegions[iter->first];
        for (decltype(cellIdxArray.size()) cell = 0; cell < cellIdxArray.size(); ++cell) {
            I64 idx = cellIdxArray[cell];
            e->mesh.indexToGlobalCoords(idx, gi, gj, gk);
            // 相对坐标
            auto i = gi - e->metadata[6], j = gj - e->metadata[7], k = gk - e->metadata[8];
            OutputData data;
            if (e->p[e->getC(i,j,k)].isEmpty)
                continue;
            data.ux = e->ux[e->getX(i,j,k)].u;
            data.uy = e->uy[e->getY(i,j,k)].u;
            data.uz = e->uz[e->getZ(i,j,k)].u;        
            data.p = e->p[e->getC(i,j,k)].data;
            data.t = e->t[e->getC(i,j,k)].data;
            data.globalid = e->id[e->getC(i,j,k)];  
            data.locationx = e->mesh.cycasMesh.xPoints[gi-1];
            data.locationy = e->mesh.cycasMesh.yPoints[gj-1];
            data.locationz = e->mesh.cycasMesh.zPoints[gk-1];              
            e->idToData[idx] = data;
        }
        e->idArray.resize(0); 
        e->result.resize(0);
        for (auto &v: e->idToData) {
            e->idArray.push_back(v.first);
            e->result.push_back(v.second); 
        }         
    }
    for (auto &e: mobs) {
        auto iter = e->mesh.cycasMesh.mobsRegions.begin();
        auto cellIdxArray = e->mesh.cycasMesh.mobsRegions[iter->first];
        for (decltype(cellIdxArray.size()) cell = 0; cell < cellIdxArray.size(); ++cell) {
            I64 idx = cellIdxArray[cell];
            e->mesh.indexToGlobalCoords(idx, gi, gj, gk);
            // 相对坐标
            auto i = gi - e->metadata[6], j = gj - e->metadata[7], k = gk - e->metadata[8];        
            OutputData data; 
            if (e->t[e->getC(i,j,k)].isEmpty) 
                continue;
            data.t = e->t[e->getC(i,j,k)].data;    
            data.locationx = e->mesh.cycasMesh.xPoints[gi-1];
            data.locationy = e->mesh.cycasMesh.yPoints[gj-1];
            data.locationz = e->mesh.cycasMesh.zPoints[gk-1];
            e->idToData[idx] = data; 
        }
        e->idArray.resize(0); 
        e->result.resize(0);
        for (auto &v: e->idToData) {
            e->idArray.push_back(v.first);
            e->result.push_back(v.second);
        } 
    }
    // std::cout << e->idArray.size() << " " << e->result.size() << std::endl;
}

#endif

// 把solver和mobssolver装到一个进程的范围里


        // for (auto &e: a) {           
        //     for (auto &p: e.pid) {
        //         MPI_Isend(sendData1[p].data(), sendData1[p].size(), MPI_DOUBLE, p, 8,
        //                 MPI_COMM_WORLD, &request[iter++]);                           
        //     }
        // }
        // for (auto &e: b) {           
        //     for (auto &p: e.pid) {
        //         MPI_Isend(sendData2[p].data(), sendData2[p].size(), MPI_DOUBLE, p, 17,
        //                 MPI_COMM_WORLD, &request[iter++]);                           
        //     }
        // } 
        // for (auto &e: a) {           
        //     for (auto &p: e.pid) {
        //         MPI_Irecv(recvData2[p].data(), recvData2[p].size(), MPI_DOUBLE, p, 17,
        //                 MPI_COMM_WORLD, &request[iter++]);                           
        //     }
        // }               
        // for (auto &e: b) {           
        //     for (auto &p: e.pid) {
        //         MPI_Irecv(recvData1[p].data(), recvData1[p].size(), MPI_DOUBLE, p, 8,
        //                 MPI_COMM_WORLD, &request[iter++]);                           
        //     }
        // }
