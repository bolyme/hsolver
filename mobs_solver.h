#ifndef TEMPERATURE
#define TEMPERATURE

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include "mpi.h"
#include "distributed_mesh.h"
#include "fluid_solver.h"

#define L 1.0                        //方腔边长，单位：米(m)
#define T0 293                     //外部温度
#define T2 313
#define T1 333                     //内部温度
#define RHO 1
                                     //时为超松弛迭代，小于1.0时为亚松弛迭代
#define NU 0.001            //粘性系数，单位:帕斯卡·秒(Pa·s)
#define CONVERGENCE2 1e-3             //全局收敛标准
#define CONVERGENCE 0.00001 
#define STEP 1              //计算域整体网格迭代上限
#define GS_STEP 100              //计算域单个网格迭代上限

#define UNDER_RELAX 0.3

#define Re ((U0*L)/NU)

struct MobSolver {
    DistributedMesh &mesh;
    std::vector<I64>& metadata;
    std::vector<LocalCell>& t;
    std::vector<LocalFace>& mx;
    std::vector<LocalFace>& my;
    std::vector<LocalFace>& mz;

    std::vector<OutputData> result; 
    std::vector<I64> idArray;
    std::map<int, OutputData> idToData;

    std::map<int, std::vector<I64>> t_boundary;
    std::map<int, std::vector<I64>> t_halo;
    std::map<int, std::vector<I64>> pr_boundary;
    std::map<int, std::vector<I64>> pr_halo;    
    std::map<int, std::vector<double>> sendData, recvData;
    std::vector<double> t_temp;

    MobSolver(DistributedMesh &m1, std::vector<I64>&m2, std::vector<LocalCell>&c,
        std::vector<LocalFace>&f1, std::vector<LocalFace>&f2, std::vector<LocalFace>&f3): 
        mesh(m1), metadata(m2), t(c), mx(f1), my(f2), mz(f3) {
            auto cx = metadata[0], cy = metadata[1], cz = metadata[2];
            auto fx = metadata[3], fy = metadata[4], fz = metadata[5];
            t_temp.resize((cx+2) * (cy+2) * (cz+2), T1);
        }
    I64 getC(I64 i, I64 j, I64 k);
    I64 getX(I64 i, I64 j, I64 k);
    I64 getY(I64 i, I64 j, I64 k);
    I64 getZ(I64 i, I64 j, I64 k);
    I64 indexT(I64 i, I64 j, I64 k);

    void applyOuterBoundaryCondition();
    void solveTemperature(std::vector<double>& t, 
                    const I64 &xStart, const I64 &xEnd, 
                    const I64 &yStart, const I64 &yEnd,
                    const I64 &zStart, const I64 &zEnd);  
    void init();
    void solve();
    void setForNext();
    void writeVTK();
    bool parseString(const std::string& token, const std::string& s, int& num1, int& num2);
    void getIJK(I64 idx, I64& i, I64& j, I64& k);
    std::vector<I64> getGlobal(I64 idx);
    bool tripletExists(const std::vector<I64>& A1, I64 e0, I64 e1, I64 e2);

};

void MobSolver::getIJK(I64 idx, I64& i, I64& j, I64& k) {
    k = idx / (metadata[0]+2) / (metadata[1]+2);
    j = (idx - k * (metadata[2]+2) * (metadata[1]+2)) / (metadata[1]+2);
    i = (idx - k * (metadata[2]+2) * (metadata[1]+2)) % (metadata[1]+2);
}
std::vector<I64> MobSolver::getGlobal(I64 idx) {
    I64 i, j, k;
    getIJK(idx, i, j, k);
    return std::vector<I64>{i-1+metadata[6], j-1+metadata[7], k-1+metadata[8]};
}

I64 MobSolver::getZ(I64 i, I64 j, I64 k) {
    return mesh.getFaceIndex(i, j, k, metadata[0], metadata[1], metadata[5]);
}
I64 MobSolver::getY(I64 i, I64 j, I64 k) {
    return mesh.getFaceIndex(i, j, k, metadata[0], metadata[4], metadata[2]);
}
I64 MobSolver::getX(I64 i, I64 j, I64 k) {
    return mesh.getFaceIndex(i, j, k, metadata[3], metadata[1], metadata[2]);
}
I64 MobSolver::getC(I64 i, I64 j, I64 k) {
    return mesh.getCellIndex(i, j, k, metadata[0], metadata[1], metadata[2]);
}
I64 MobSolver::indexT(I64 i, I64 j, I64 k) {
    return mesh.getCellIndex(i, j, k, metadata[0]+2, metadata[1]+2, metadata[2]+2);
}
void MobSolver::applyOuterBoundaryCondition() {
    //先假设是全局,后面加上if判断
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];
    for (int k = 0; k < cz; ++k) {
        for (int j = 0; j < cy; ++j) {
            for (int i = 0; i < cx; ++i) {
                t[getC(i, j, k)].data = T1;
            }
        }
    }
}
void MobSolver::init() {

    applyOuterBoundaryCondition();        
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];     
    for (int k = 0; k < cz; ++k) {
        for (int j = 0; j < cy; ++j) {
            for (int i = 0; i < fx; ++i) {
                int rank;
                auto pid = rank;
                if (mx[getX(i, j, k)].isBoundary && parseString("pr_", mx[getX(i, j, k)].info->name, rank, pid)) {
                    // std::cout << rank << " " << pid << std::endl;
                    if (pr_boundary.find(pid) == pr_boundary.end()) 
                        pr_boundary[pid] = std::vector<I64>();
                    if (pr_halo.find(pid) == pr_halo.end())
                        pr_halo[pid] = std::vector<I64>();

                    if (i + 1 >= fx) {
                        pr_boundary[pid].push_back(indexT(i-1+1, j+1, k+1));
                        pr_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (i - 1 < 0) {
                        pr_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        pr_halo[pid].push_back(indexT(i-1+1, j+1, k+1));
                    } 
                    else if (t[getC(i, j, k)].isEmpty) {
                        pr_boundary[pid].push_back(indexT(i-1+1, j+1, k+1));
                        pr_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (t[getC((i-1), (j), (k))].isEmpty) {
                        pr_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        pr_halo[pid].push_back(indexT(i-1+1, j+1, k+1));
                    }
                } else if (mx[getX(i, j, k)].isBoundary && parseString("ib_", mx[getX(i, j, k)].info->name, rank, pid)) {
                    if (t_boundary.find(pid) == t_boundary.end()) 
                        t_boundary[pid] = std::vector<I64>();
                    if (t_halo.find(pid) == t_halo.end())
                        t_halo[pid] = std::vector<I64>();

                    if (i + 1 >= fx) {
                        t_boundary[pid].push_back(indexT(i-1+1, j+1, k+1));
                        t_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (i - 1 < 0) {
                        t_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        t_halo[pid].push_back(indexT(i-1+1, j+1, k+1));
                        // std::cout << pid << " " << indexT(i+1, j+1, k+1) << " "<<  indexT(i-1+1, j+1, k+1) << std::endl;
                    } else if (t[getC(i, j, k)].isEmpty) {
                        t_boundary[pid].push_back(indexT(i-1+1, j+1, k+1));
                        t_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (t[getC((i-1), (j), (k))].isEmpty) {
                        t_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        t_halo[pid].push_back(indexT(i-1+1, j+1, k+1));
                    }                   
                }
            }
        }
    }
    
    for (int k = 0; k < cz; ++k) {
        for (int j = 0; j < fy; ++j) {
            for (int i = 0; i < cx; ++i) {
                int rank;
                auto pid = rank;
                if (my[getY(i, j, k)].isBoundary && parseString("pr_", my[getY(i, j, k)].info->name, rank, pid)) {
                    // std::cout << rank << " " << pid << std::endl;
                    if (pr_boundary.find(pid) == pr_boundary.end()) 
                        pr_boundary[pid] = std::vector<I64>();
                    if (pr_halo.find(pid) == pr_halo.end())
                        pr_halo[pid] = std::vector<I64>();

                    if (j + 1 >= fy) {
                        pr_boundary[pid].push_back(indexT(i+1, j-1+1, k+1));
                        pr_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (j - 1 < 0) {
                        pr_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        pr_halo[pid].push_back(indexT(i+1, j-1+1, k+1));
                    } else if (t[getC(i, j, k)].isEmpty) {
                        pr_boundary[pid].push_back(indexT(i+1, j-1+1, k+1));
                        pr_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (t[getC((i), (j-1), (k))].isEmpty) {
                        pr_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        pr_halo[pid].push_back(indexT(i+1, j-1+1, k+1));
                    } 
                } else if (my[getY(i, j, k)].isBoundary && parseString("ib_", my[getY(i, j, k)].info->name, rank, pid)) {
                    if (t_boundary.find(pid) == t_boundary.end()) 
                        t_boundary[pid] = std::vector<I64>();
                    if (t_halo.find(pid) == t_halo.end())
                        t_halo[pid] = std::vector<I64>();

                    if (j + 1 >= fy) {
                        t_boundary[pid].push_back(indexT(i+1, j-1+1, k+1));
                        t_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (j - 1 < 0) {
                        t_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        t_halo[pid].push_back(indexT(i+1, j-1+1, k+1));
                    } else if (t[getC(i, j, k)].isEmpty) {
                        t_boundary[pid].push_back(indexT(i+1, j-1+1, k+1));
                        t_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (t[getC((i), (j-1), (k))].isEmpty) {
                        t_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        t_halo[pid].push_back(indexT(i+1, j-1+1, k+1));
                    }                  
                }
            }
        }
    }
    
    for (int k = 0; k < fz; ++k) {
        for (int j = 0; j < cy; ++j) {
            for (int i = 0; i < cx; ++i) {
                int rank;
                auto pid = rank;
                if (mz[getZ(i, j, k)].isBoundary && parseString("pr_", mz[getZ(i, j, k)].info->name, rank, pid)) {
                    // std::cout << rank << " " << pid << std::endl;
                    if (pr_boundary.find(pid) == pr_boundary.end()) 
                        pr_boundary[pid] = std::vector<I64>();
                    if (pr_halo.find(pid) == pr_halo.end())
                        pr_halo[pid] = std::vector<I64>();

                    if (k + 1 >= fz) {
                        pr_boundary[pid].push_back(indexT(i+1, j+1, k-1+1));
                        pr_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (k - 1 < 0) {
                        pr_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        pr_halo[pid].push_back(indexT(i+1, j+1, k-1+1));
                    } else if (t[getC(i, j, k)].isEmpty) {
                        pr_boundary[pid].push_back(indexT(i+1, j+1, k-1+1));
                        pr_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (t[getC((i), (j), (k-1))].isEmpty) {
                        pr_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        pr_halo[pid].push_back(indexT(i+1, j+1, k-1+1));
                    }
                } else if (mz[getZ(i, j, k)].isBoundary && parseString("ib_", mz[getZ(i, j, k)].info->name, rank, pid)) {
                    if (t_boundary.find(pid) == t_boundary.end()) 
                        t_boundary[pid] = std::vector<I64>();
                    if (t_halo.find(pid) == t_halo.end())
                        t_halo[pid] = std::vector<I64>();

                    if (k + 1 >= fz) {
                        t_boundary[pid].push_back(indexT(i+1, j+1, k-1+1));
                        t_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (k - 1 < 0) {
                        t_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        t_halo[pid].push_back(indexT(i+1, j+1, k-1+1));
                    } else if (t[getC(i, j, k)].isEmpty) {
                        t_boundary[pid].push_back(indexT(i+1, j+1, k-1+1));
                        t_halo[pid].push_back(indexT(i+1, j+1, k+1));
                    } else if (t[getC((i), (j), (k-1))].isEmpty) {
                        t_boundary[pid].push_back(indexT(i+1, j+1, k+1));
                        t_halo[pid].push_back(indexT(i+1, j+1, k-1+1));
                    }                  
                }
            }
        }
    }

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (auto &e: t_boundary) {
        // std::sort(e.second.begin(), e.second.end());
        // auto iter = std::unique(e.second.begin(), e.second.end());
        // e.second.resize(iter - e.second.begin());
        std::cout << "mobs t boundary size" << " " << rank << " " << e.first << " " << e.second.size() << std::endl;
    }            
    for (auto &e: t_halo) {
        // std::sort(e.second.begin(), e.second.end());
        // auto iter = std::unique(e.second.begin(), e.second.end());
        // e.second.resize(iter - e.second.begin());
        std::cout << "mobs t halo size" << " " << rank << " " << e.first << " " << e.second.size() << std::endl; 
    } 
    for (auto &e: t_boundary) {
        sendData[e.first] = std::vector<double>();
        recvData[e.first] = std::vector<double>();
    }    
}

void MobSolver::setForNext() {
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];  
    for (int k = 0; k < cz; ++k) {
        for (int j = 0; j < cy; ++j) {
            for (int i = 0; i < cx; ++i) {
                if (t[getC(i,j,k)].isEmpty) continue;
                t[getC(i,j,k)].data = t_temp[indexT(i+1,j+1,k+1)];
            }
        }
    }    
}

void MobSolver::solve() {
    // set global boundary
    applyOuterBoundaryCondition();        
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];    
      
    for (int iter = 0; iter < STEP; ++iter){

        // internal cell 
        // local 0~cy-1 - 1~cy
        // absolute 0~cy+1

        I64 xStart = 1, xEnd = cx;
        I64 yStart = 1, yEnd = cy;
        I64 zStart = 1, zEnd = cz;
        // std::cout << "start mobs temp" << std::endl;
        solveTemperature(t_temp, xStart, xEnd, yStart, yEnd, zStart, zEnd);
        // 通信固体温度
        for (auto &e: sendData) {
            e.second.resize(0);
        }
        for (auto &e: recvData) {
            e.second.resize(0);
        }
        for (auto &v: pr_boundary) {
            for (auto &e: v.second) {
                sendData[v.first].push_back(t_temp[e]);
            }
            recvData[v.first].resize(sendData[v.first].size()); 
        }
        for (auto &e: pr_boundary) {          
            MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 8,
                         recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 8,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (auto &v: pr_halo) {
            for (int e = 0; e < v.second.size(); ++e) {
                t_temp[pr_halo[v.first][e]] = recvData[v.first][e];
                if (recvData[v.first][e] < 200 || recvData[v.first][e] > 400)
                    std::cout << v.first << " " << recvData[v.first][e] << std::endl ;
            }
        } 
        // 通信流体温度
        // for (auto &v: t_boundary) {
        //     for (auto &e: v.second) {
        //         sendData[v.first].push_back(t_temp[e]);
        //         // if (v.first == 2) std::cout << e << " ";
        //     }
        //     recvData[v.first].resize(sendData[v.first].size()); 
        // }
        // for (auto &e: t_boundary) {
        //     // for (auto value: sendData[e.first])
        //     // std::cout << rank << " ";
        //     // std::cout << std::endl;            
        //     MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 8,
        //                  recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 8,
        //                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // }
        // // std::cout << std::endl;
        // for (auto &v: t_halo) {
        //     for (int e = 0; e < v.second.size(); ++e) {
        //         t_temp[t_halo[v.first][e]] = recvData[v.first][e];
        //         if (recvData[v.first][e] < 200 || recvData[v.first][e] > 400)
        //             std::cout << v.first << " " << recvData[v.first][e] << std::endl ;
        //     }
        // } 
        
        // applyOuterBoundaryCondition();
    }
    // writeVTK(); 
    
    // I64 gi, gj, gk;
    // auto iter = mesh.cycasMesh.fluidRegions.begin();
    // auto cellIdxArray = mesh.cycasMesh.fluidRegions[iter->first];
    // for (decltype(cellIdxArray.size()) cell = 0; cell < cellIdxArray.size(); ++cell) {
    //     I64 idx = cellIdxArray[cell];
    //     mesh.indexToGlobalCoords(idx, gi, gj, gk);
    //     // 相对坐标
    //     auto i = gi - metadata[6], j = gj - metadata[7], k = gk - metadata[8];
    //     OutputData data;                    
    //     idToData[idx] = data;
    // }
    // iter = mesh.cycasMesh.mobsRegions.begin();
    // cellIdxArray = mesh.cycasMesh.mobsRegions[iter->first];
    // for (decltype(cellIdxArray.size()) cell = 0; cell < cellIdxArray.size(); ++cell) {
    //     I64 idx = cellIdxArray[cell];
    //     mesh.indexToGlobalCoords(idx, gi, gj, gk);
    //     // 相对坐标
    //     auto i = gi - metadata[6], j = gj - metadata[7], k = gk - metadata[8];        
    //     OutputData data; 
    //     if (t[getC(i,j,k)].isEmpty) 
    //         continue;
    //     data.t = t[getC(i,j,k)].data;    
    //     // std::cout << data.t << " " ;
    //     idToData[idx] = data; 
    // }

    // idArray.resize(0); 
    // result.resize(0);
    // for (auto &e: idToData) {
    //     idArray.push_back(e.first);
    //     result.push_back(e.second);
    //     // if (rank == 0) std::cout << e.first << " " << e.second.uz << " " << e.second.p << std::endl;
    // }   
    // std::cout << idArray.size() << " " << result.size() << std::endl;       
}

void MobSolver::writeVTK() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    std::ofstream fout("mobs" + std::to_string(rank) + ".vtk");

    fout << "# vtk DataFile Version 3.0\n";
    fout << "Vector field example\n";
    fout << "ASCII\n";
    fout << "DATASET STRUCTURED_POINTS\n";
    fout << "DIMENSIONS " << metadata[0] << " " << metadata[1] << " " <<  metadata[2] << "\n";
    fout << "ORIGIN 4 2 4\n";
    fout << "SPACING 1 1 1\n";
    fout << "POINT_DATA " << metadata[0] * metadata[1] * metadata[2] << "\n";
    fout << "VECTORS velocity float\n";

    fout << std::fixed << std::setprecision(6);
    for(int k = 0; k < metadata[2]; ++k) {
        for(int j = 0; j < metadata[1]; ++j) {
            for(int i = 0; i < metadata[0]; ++i) {
                fout << 0 << " "
                     << 0 << " "
                     << 0 << std::endl;
            }
        }
    }  
}


bool MobSolver::parseString(const std::string& token, const std::string& s, int& num1, int& num2) {
    // 检查前缀和最小长度
    if (s.substr(0, 3) != token) return false;

    std::string rest = s.substr(3); // 截取"pr_"后的部分
    size_t split_pos = rest.find('_');

    // 检查下划线是否存在且位置合法
    if (split_pos == std::string::npos || split_pos == 0 || split_pos == rest.size() - 1) {
        return false;
    }

    std::string num1_str = rest.substr(0, split_pos);
    std::string num2_str = rest.substr(split_pos + 1);

    // 检查是否为纯数字
    for (char c : num1_str) if (!isdigit(c)) return false;
    for (char c : num2_str) if (!isdigit(c)) return false;

    try {
        num1 = stoi(num1_str);
        num2 = stoi(num2_str);
    } catch (...) { // 处理转换异常（如溢出）
        return false;
    }

    return true;
}

void MobSolver::solveTemperature(std::vector<double> &t_temp,
                                     const I64 &xStart, const I64 &xEnd, 
                                     const I64 &yStart, const I64 &yEnd, 
                                     const I64 &zStart, const I64 &zEnd) {
    double residual2 = 1;
    for (int iter = 0; iter < GS_STEP; ++iter){
        // std::cout << "residual2: " << residual2 << std::endl;
        if (residual2 < CONVERGENCE2) {
            break;
        }
        for (int k = zStart; k <= zEnd; ++k) {
            for (int j = yStart; j <= yEnd; ++j) {
                for (int i = xStart; i <= xEnd; ++i) {
                    if (t[getC(i-1,j-1,k-1)].isEmpty) {
                        // std::cout << "empty cell.\n";
                        if (t[getC(i-1,j-1,k-1)].data < 1e-5 || t[getC(i-1,j-1,k-1)].data > 1000)
                            t[getC(i-1,j-1,k-1)].data = 293;
                        continue;
                    }
                    double old = t_temp[indexT(i,j,k)];
                    t_temp[indexT(i,j,k)] = 1.0 / 6.0 * (t_temp[indexT(i+1,j,k)] + t_temp[indexT(i-1,j,k)]
                                                   +t_temp[indexT(i,j-1,k)] + t_temp[indexT(i,j+1,k)]
                                                   +t_temp[indexT(i,j,k-1)] + t_temp[indexT(i,j,k+1)]);
                    t_temp[indexT(i,j,k)] = t_temp[indexT(i,j,k)] * UNDER_RELAX + (1-UNDER_RELAX) * old;
                    // std::cout << t_temp[indexT(i+1,j,k)] << " " << t_temp[indexT(i-1,j,k)] << " "
                    //           << t_temp[indexT(i,j-1,k)] << " " << t_temp[indexT(i,j+1,k)] << " "
                    //           << t_temp[indexT(i,j,k-1)] << " " << t_temp[indexT(i,j,k+1)] << std::endl;
                    residual2 += std::pow(t_temp[indexT(i,j,k)] - old, 2);
                }
            }
        }
        residual2 = std::sqrt(residual2 / (zEnd - zStart + 1) / (yEnd - yStart + 1) / (xEnd - xStart + 1));
    }
}

bool MobSolver::tripletExists(const std::vector<I64>& A1, I64 e0, I64 e1, I64 e2) {
    for (size_t i = 0; i + 2 < A1.size(); i += 3) {
        if (A1[i] == e0 && A1[i + 1] == e1 && A1[i + 2] == e2) {
            return true;
        }
    }
    return false;
}


#endif
