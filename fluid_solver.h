#ifndef DIFFERENCE
#define DIFFERENCE

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

#define L 1.0               
#define U0 0.1                    
#define RHO 1
                                 
#define NU 0.001            
#define CONVERGENCE2 1e-4             
#define CONVERGENCE 2e-3
#define STEP 1        
#define GS_STEP 100               
#define ALPHA 0.001

#define X_SPACING 0.02
#define Y_SPACING 0.02
#define Z_SPACING 0.02
#define UNDER_RELAX 0.3

#define Re ((U0*L)/NU)
struct OutputData {
    double ux=0, uy=0, uz=0, p=0;
    double t=0;
    double locationx=0.0;
    double locationy=0.0;
    double locationz=0.0;
    int globalid=-1;
};
struct KeyHash {
    std::size_t operator()(const std::tuple<int, int, int>& key) const {
        auto [i, j, k] = key;
        return std::hash<int>()(i) ^ (std::hash<int>()(j) << 1) ^ (std::hash<int>()(k) << 2);
    }
};
struct FluidSolver {
    DistributedMesh &mesh;
    std::vector<I64>& metadata;
    std::vector<LocalCell>& p;
    std::vector<LocalFace>& ux;
    std::vector<LocalFace>& uy; 
    std::vector<LocalFace>& uz;
    // MPI_Comm* hypre_comm;
    std::vector<LocalCell> t;
    std::vector<int> id;

    std::map<int, std::vector<I64>> u_boundary, v_boundary, w_boundary, p_boundary, t_boundary;
    std::map<int, std::vector<I64>> u_halo, v_halo, w_halo, p_halo, t_halo;
    std::map<int, std::vector<double>> sendData, recvData;
    std::vector<double> ulocation, vlocation, wlocation;
    std::vector<double> u_star, v_star, w_star, up, vp, wp, p_temp, p_prime, t_temp;
    std::vector<double> u_param, v_param, w_param, p_param;
    std::vector<OutputData> result;
    std::vector<I64> idArray;
    std::map<int, OutputData> idToData;
    std::unordered_map<std::tuple<int, int, int>, int, KeyHash> global_map;
    std::unordered_map<int, std::tuple<int, int, int>> result_map;  

    std::vector<int> ncols;
    std::vector<HYPRE_Int> anotherRow;
    std::vector<std::vector<int>> flag;
    std::vector<HYPRE_Int> cols1d;
    std::vector<HYPRE_Complex> values1d;
    std::vector<std::vector<double>> values;
    std::vector<std::vector<int>> cols;
    std::vector<HYPRE_Int> rows;
    std::vector<double> b_values, x_values;

    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_IJVector b;
    HYPRE_ParVector par_b;
    HYPRE_IJVector x;
    HYPRE_ParVector par_x; 
    HYPRE_Solver solver, precond;    

    FluidSolver(DistributedMesh &m1, std::vector<I64>&m2, std::vector<LocalCell>&c, 
        std::vector<LocalFace>&f1, std::vector<LocalFace>&f2, std::vector<LocalFace>&f3): 
        mesh(m1), metadata(m2), p(c), ux(f1), uy(f2), uz(f3) { 
            t = p; 
            id.resize(p.size()); 
            auto cx = metadata[0], cy = metadata[1], cz = metadata[2];
            auto fx = metadata[3], fy = metadata[4], fz = metadata[5]; 
            p_temp.resize((cx+2) * (cy+2) * (cz+2), 0);
            p_prime.resize(p_temp.size(),0);   
            u_star.resize((fx+2) * (cy+2) * (cz+2), 0);
            v_star.resize((cx+2) * (fy+2) * (cz+2), 0);
            w_star.resize((cx+2) * (cy+2) * (fz+2), 0);
            up.resize(u_star.size());
            vp.resize(v_star.size());
            wp.resize(w_star.size());
            t_temp.resize(p_temp.size(), 293);
        }
    ~FluidSolver() {    
        HYPRE_ParCSRPCGDestroy(solver);
        HYPRE_IJVectorDestroy(b);
        HYPRE_IJVectorDestroy(x);
        HYPRE_IJMatrixDestroy(A);
        HYPRE_Finalize();
    }

    I64 getZ(I64 i, I64 j, I64 k);
    I64 getY(I64 i, I64 j, I64 k);
    I64 getX(I64 i, I64 j, I64 k);
    I64 getC(I64 i, I64 j, I64 k);
    I64 indexU(I64 i, I64 j, I64 k);
    I64 indexV(I64 i, I64 j, I64 k);
    I64 indexW(I64 i, I64 j, I64 k);
    I64 indexP(I64 i, I64 j, I64 k);

    void applyXBoundaryCondition(std::vector<double>& u_star);
    void applyYBoundaryCondition(std::vector<double>& v_star);
    void applyZBoundaryCondition(std::vector<double>& w_star);
    void applyCellBoundaryCondition(std::vector<double>& p_temp);
    void solve();
    void init();
    void initParam();
    void setForNext();

    void solveXmomentum(std::vector<double> &u_star, std::vector<double> &v_star, std::vector<double> &w_star,
                   std::vector<double>& p_temp, std::vector<double>& up, const I64 &xStart, const I64 &xEnd, const I64 &yStart, const I64 &yEnd,
                   const I64 &zStart, const I64 &zEnd);

    void solveYmomentum(std::vector<double> &u_star, std::vector<double> &v_star, std::vector<double> &w_star, 
                   std::vector<double>& p_temp, std::vector<double>& vp, const I64 &xStart, const I64 &xEnd, const I64 &yStart, const I64 &yEnd,
                   const I64 &zStart, const I64 &zEnd);

    void solveZmomentum(std::vector<double> &u_star, std::vector<double> &v_star, std::vector<double> &w_star,
                    std::vector<double>& p_temp, std::vector<double>& wp, const I64 &xStart, const I64 &xEnd, const I64 &yStart, const I64 &yEnd,
                   const I64 &zStart, const I64 &zEnd);

    void solveCorrectionPressure(const std::vector<double> &u_star, const std::vector<double> &v_star,
                            const std::vector<double> &w_star, std::vector<double> &p_prime,
                            std::vector<double> &up, std::vector<double> &vp, std::vector<double> &wp, 
                            const I64 &xStart, const I64 &xEnd, const I64 &yStart, const I64 &yEnd, const I64 &zStart,
                            const I64 &zEnd);

    void HYPREIJMatrix_SolveCorrectionPressure(const std::vector<double> &u_star, const std::vector<double> &v_star,
                                     const std::vector<double> &w_star, std::vector<double> &p_prime,
                                     std::vector<double> & up, std::vector<double> & vp, std::vector<double> & wp,
                                     I64 cx, I64 cy, I64 cz, I64 offsetX, I64 offsetY, I64 offsetZ, I64 globalStart, I64 local_size); 
    int global_idx(int i, int j, int k, int nx, int ny, int nz, std::unordered_map<std::tuple<int, int, int>, int, KeyHash>&);
    void solveTemp(const std::vector<double> &u_star, const std::vector<double> &v_star,
                                     const std::vector<double> &w_star, std::vector<double> &t_temp,
                                     const I64 &xStart, const I64 &xEnd, 
                                     const I64 &yStart, const I64 &yEnd, 
                                     const I64 &zStart, const I64 &zEnd);                        

    bool parseString(const std::string& token, const std::string& s, int& num1, int& num2);

    void initMap(int size, int xstart, int xend, int ystart, int yend, int zstart, int zend,
                 std::unordered_map<std::tuple<int, int, int>, int, KeyHash> &local_map, int ilower, int iter);

    void initMatrix(I64 local_size, int ilower, int iupper);

    void
    initColsAndSetValues(I64 cx, I64 cy, I64 cz, I64 local_size, int rank, int nx, int ny, int nz, int xstart, int xend,
                         int ystart, int yend, int zstart, int zend, int ilower, int nnz,
                         const std::vector<double> &u_star, const std::vector<double> &v_star, const std::vector<double> &w_star, 
                         std::vector<double> &up, std::vector<double> &vp, std::vector<double> &wp);
    void setCoffValues(const std::vector<double> &u_star, const std::vector<double> &v_star,
                                    const std::vector<double> &w_star, const std::vector<double> &up,
                                    const std::vector<double> &vp, const std::vector<double> &wp, I64 local_size,
                                    int xstart, int xend, int ystart, int yend, int zstart, int zend, int ilower);

};

I64 FluidSolver::getZ(I64 i, I64 j, I64 k) {
    return mesh.getFaceIndex(i, j, k, metadata[0], metadata[1], metadata[5]);
}
I64 FluidSolver::indexW(I64 i, I64 j, I64 k) {
    return mesh.getFaceIndex(i, j, k, metadata[0]+2, metadata[1]+2, metadata[5]+2);
}
I64 FluidSolver::getY(I64 i, I64 j, I64 k) {
    return mesh.getFaceIndex(i, j, k, metadata[0], metadata[4], metadata[2]);
}
I64 FluidSolver::indexV(I64 i, I64 j, I64 k) {
    return mesh.getFaceIndex(i, j, k, metadata[0]+2, metadata[4]+2, metadata[2]+2);
}
I64 FluidSolver::getX(I64 i, I64 j, I64 k) {
    return mesh.getFaceIndex(i, j, k, metadata[3], metadata[1], metadata[2]);
}
I64 FluidSolver::indexU(I64 i, I64 j, I64 k) {
    return mesh.getFaceIndex(i, j, k, metadata[3]+2, metadata[1]+2, metadata[2]+2);
}
I64 FluidSolver::getC(I64 i, I64 j, I64 k) {
    return mesh.getCellIndex(i, j, k, metadata[0], metadata[1], metadata[2]);
}
I64 FluidSolver::indexP(I64 i, I64 j, I64 k) {
    return mesh.getCellIndex(i, j, k, metadata[0]+2, metadata[1]+2, metadata[2]+2);
}

void FluidSolver::applyXBoundaryCondition(std::vector<double>& u_star) {
// 1. i == fx for ux and u_star
// 2. i == 0 and i == fx + 1 for u_star
// 3. v_star and w_star for i == 0 and i = fx + 1
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];
    // x区域最小halo  
    if (metadata[6] == 1) {
        for (int k = 0; k < cz + 2; ++k) { // k: 0-51 for u_star 1-50 for ux
            for (int j = 0; j < cy + 2; ++j) { // j: 0-51 for u_star 1-50 for ux
                u_star[indexU(0,j, k)] = 0;
                v_star[indexV(0,j, k)] = 0;
                w_star[indexW(0,j, k)] = 0;
            }
        }
    // x区域最大halo    
    } else if (metadata[6]-1+fx == mesh.cycasMesh.Nx+1) {
        for (int k = 0; k < cz + 2; ++k) { // k: 0-51 for u_star 1-50 for ux
            for (int j = 0; j < cy + 2; ++j) { // j: 0-51 for u_star 1-50 for ux
                u_star[indexU(fx+1, j, k)] = 0;
                v_star[indexV(cx+1,j, k)] = 0;
                w_star[indexW(cx+1,j, k)] = 0;
            }
        }
    }       
    // x区域最大boundary        
    for (int k = 0; k < cz; ++k) { // 0-51 for u_star 1-50 for ux
        for (int j = 0; j < cy; ++j) { // 0-51 for u_star 1-50 for ux
            if (metadata[6]-1+fx-1 == mesh.cycasMesh.Nx) {
                ux[getX(fx-1, j, k)].u = 0;
                u_star[indexU(fx,j+1,k+1)] = 0;
            } 
        }
    }          
    /*
      for (int k = 0; k < cz; ++k) { // 0-51 for u_star 1-50 for ux
        for (int j = 0; j < cy; ++j) { // 0-51 for u_star 1-50 for ux
            for (int i = 0; i < fx; ++i) {
                if (metadata[7]-1+j == mesh.cycasMesh.Ny-1) {
                    ux[getX(i, j, k)].u = U0;
                    u_star[indexU(i+1,j+1,k+1)] = U0;
                }            
            }
        }
    } */ 
}

void FluidSolver::applyZBoundaryCondition(std::vector<double>& w_star) {
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];  
    // z区域最小halo
    if (metadata[8] == 1) {
        for (int j = 0; j < cy + 2; ++j) { // k: 0-51 for u_star 1-50 for ux
            for (int i = 0; i < cx + 2; ++i) { // j: 0-51 for u_star 1-50 for ux
                u_star[indexU(i,j, 0)] = 0;
                v_star[indexV(i,j, 0)] = 0;
                w_star[indexW(i,j, 0)] = 0;
            }
        }
    // z区域最大halo    
    } else if (metadata[8]-1+fz == mesh.cycasMesh.Nz+1) {
        for (int j = 0; j < cy + 2; ++j) { // k: 0-51 for u_star 1-50 for ux
            for (int i = 0; i < cx + 2; ++i) { // j: 0-51 for u_star 1-50 for ux
                u_star[indexU(i,j, cz+1)] = 0;
                v_star[indexV(i,j, cz+1)] = 0;
                w_star[indexW(i,j, fz+1)] = 0;
            }
        }
    }       
    // z区域最大boundary        
    for (int j = 0; j < cy; ++j) { // k: 0-51 for u_star 1-50 for ux
        for (int i = 0; i < cx; ++i) { // j: 0-51 for u_star 1-50 for ux
            if (metadata[8]-1+fz-1 == mesh.cycasMesh.Nz) {
                uz[getZ(i, j, fz-1)].u = 0;
                w_star[indexW(i+1,j+1,fz)] = 0;
            } 
        }
    }
}

void FluidSolver::applyYBoundaryCondition(std::vector<double>& v_star) {
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];  
    // y区域最小halo
    if (metadata[7] == 1) {
        for (int k = 0; k < cz + 2; ++k) { // 0-51 for u_star 1-50 for ux
            for (int i = 0; i < cx + 2; ++i) { // 0-52 for u_star 1-51 for ux        
                u_star[indexU(i,0, k)] = 0;
                v_star[indexV(i,0, k)] = 0;
                w_star[indexW(i,0, k)] = 0;    
            }        
        }
    // y区域最大halo    
    } else if (metadata[7]-1+fy == mesh.cycasMesh.Ny+1) {
        for (int k = 0; k < cz + 2; ++k) { // 0-51 for u_star 1-50 for ux
            for (int i = 0; i < cx + 2; ++i) { // 0-52 for u_star 1-51 for ux        
                u_star[indexU(i,cy+1, k)] = U0;
                v_star[indexV(i,fy+1, k)] = 0;
                w_star[indexW(i,cy+1, k)] = 0;    
            }        
        }        
    }
    // y区域最大boundary
    for (int k = 0; k < cz; ++k) { // 0-51 for u_star 1-50 for ux
        for (int i = 0; i < cx; ++i) { // 0-52 for u_star 1-51 for ux
            if (metadata[7]-1+fy-1 == mesh.cycasMesh.Ny) {
                uy[getY(i, fy-1, k)].u = 0;
                v_star[indexV(i+1,fy,k+1)] = 0;
                ux[getX(i, cy-1, k)].u = U0;
                u_star[indexU(i+1,cy,k+1)] = U0;                     
            }
        }
    }
      
}
void FluidSolver::applyCellBoundaryCondition(std::vector<double>& p_temp) {
    //先假设是全局,后面加上if判断
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];    
    for (int k = 0; k < cz; ++k) {
        for (int j = 0; j < cy; ++j) {
            for (int i = 0; i < cx; ++i) {
                auto gx = i+metadata[6]-1;
                auto gy = j+metadata[7]-1; 
                auto gz = k+metadata[8]-1;
                // if (gx == 0 || gy == 0 || gz == 0 || gx == mesh.cycasMesh.Nx - 1 || gz == mesh.cycasMesh.Nz - 1) {
                //     p[getC(i, j, k)].data = 0.0;
                // } else if (gy == mesh.cycasMesh.Ny - 1) {
                //     p[getC(i, j, k)].data = 0.01;
                // }
                if (gx == 0)
                    p[getC(i, j, k)].data = p[getC(i+1, j, k)].data; 
                else if (gx == mesh.cycasMesh.Nx - 1)
                    p[getC(i, j, k)].data = p[getC(i-1, j, k)].data;  
                else if (gy == 0) 
                    p[getC(i, j, k)].data = p[getC(i, j+1, k)].data;
                else if (gy == mesh.cycasMesh.Ny - 1)
                    p[getC(i, j, k)].data = p[getC(i, j-1, k)].data;
                else if (gz == 0) 
                    p[getC(i, j, k)].data = p[getC(i, j, k+1)].data;
                else if (gz == mesh.cycasMesh.Nz - 1)
                    p[getC(i, j, k)].data = p[getC(i, j, k-1)].data;  

                if (gx == 0) {
                    p_temp[indexP(i+1,j+1,k+1)] = p[getC(i, j, k)].data;
                    p_temp[indexP(i, j+1,k+1)] = p[getC(i, j, k)].data;
                } 
                else if (gx == mesh.cycasMesh.Nx - 1) {
                    p_temp[indexP(i+1,j+1,k+1)] = p[getC(i, j, k)].data; 
                    p_temp[indexP(i+2,j+1,k+1)] = p[getC(i, j, k)].data;
                }else if (gy == 0) {
                    p_temp[indexP(i+1,j+1,k+1)] = p[getC(i, j, k)].data;
                    p_temp[indexP(i+1, j,k+1)] = p[getC(i, j, k)].data;
                }else if (gy == mesh.cycasMesh.Ny - 1){
                    p_temp[indexP(i+1,j+1,k+1)] = p[getC(i, j, k)].data;
                    p_temp[indexP(i+1,j+2,k+1)] = p[getC(i, j, k)].data;
                    if (i+1 == cx) {
                        p_temp[indexP(i+2,j+2,k+1)] = p[getC(i, j, k)].data;
                    }
                    if (k+1 == cz) {
                        p_temp[indexP(i+1,j+2,k+2)] = p[getC(i, j, k)].data;
                    }
                    if (i+1 == 1) {
                        p_temp[indexP(i,j+2,k+1)] = p[getC(i, j, k)].data;
                    }
                    if (k+1 == 1) {
                        p_temp[indexP(i+1,j+2,k)] = p[getC(i, j, k)].data;
                    }
                    if (i+1 == cx && k+1 == cz) {
                        p_temp[indexP(i+2,j+2,k+2)] = p[getC(i, j, k)].data;
                    }
                    if (i+1 == 1 && k+1 ==1) {
                        p_temp[indexP(i,j+2,k)] = p[getC(i, j, k)].data;
                    }
                }else if (gz == 0) {
                    p_temp[indexP(i+1,j+1,k+1)] = p[getC(i, j, k)].data;
                    p_temp[indexP(i+1,j+1,k)] = p[getC(i, j, k)].data;
                }else if (gz == mesh.cycasMesh.Nz - 1){
                    p_temp[indexP(i+1,j+1,k+1)] = p[getC(i, j, k)].data;
                    p_temp[indexP(i+1,j+1,k+2)] = p[getC(i, j, k)].data;   
                }                 
            }
        }
    }    
    p[getC(0,0,0)].data = 0;
    for (int k = 0; k < cz; ++k) {
        for (int j = 0; j < cy; ++j) {
            for (int i = 0; i < cx; ++i) {
                auto gx = i+metadata[6]-1;
                auto gy = j+metadata[7]-1; 
                auto gz = k+metadata[8]-1;
                if (gx == 0 || gy == 0 || gz == 0 || gx == mesh.cycasMesh.Nx - 1 || gz == mesh.cycasMesh.Nz - 1 || gy == mesh.cycasMesh.Ny - 1)
                    t[getC(i, j, k)].data = 293;
            }
        }
    }                   
}

void FluidSolver::solveZmomentum(std::vector<double> &u_star, std::vector<double> &v_star, std::vector<double> &w_star,
    std::vector<double>& p_temp, std::vector<double>& wp, const I64 &xStart, const I64 &xEnd, const I64 &yStart,
                            const I64 &yEnd, const I64 &zStart, const I64 &zEnd) {
    double residual = 1, global_residual = residual;
    int step = 0, left_step = step;
    int iter = 0;
    while (step < GS_STEP) {
        if (global_residual < CONVERGENCE) {
            // left_step = GS_STEP - step;
            break;
        } 
        iter = 0;    
        for (int k = zStart; k <= zEnd-1; ++k) {
            for (int j = yStart; j <= yEnd; ++j) {
                for (int i = xStart; i <= xEnd; ++i) {
                    if (p[getC(i-1, j-1, k-1)].isEmpty) {
                        continue;
                    }     

                    double facexy = w_param[iter*9];
                    double facexz = w_param[iter*9+1];
                    double faceyz = w_param[iter*9+2];

                    double fb = 0.5 * RHO * facexy * (w_star[indexW(i, j, k)] + w_star[indexW(i, j, k - 1)]);
                    double ff = 0.5 * RHO * facexy * (w_star[indexW(i, j, k)] + w_star[indexW(i, j, k + 1)]);
                    double fw = 0.5 * RHO * faceyz * (u_star[indexU(i, j, k)] + u_star[indexU(i, j, k - 1)]);
                    double fe = 0.5 * RHO * faceyz * (u_star[indexU(i + 1, j, k)] + u_star[indexU(i + 1, j, k - 1)]);
                    double fs = 0.5 * RHO * facexz * (v_star[indexV(i, j, k)] + v_star[indexV(i, j, k - 1)]);
                    double fn = 0.5 * RHO * facexz * (v_star[indexV(i, j + 1, k)] + v_star[indexV(i, j + 1, k - 1)]);

                    double Dw = w_param[iter*9+3], De = w_param[iter*9+4];
                    double Ds = w_param[iter*9+5], Dn = w_param[iter*9+6];
                    double Db = w_param[iter*9+7], Df = w_param[iter*9+8];

                    double ae = De + std::max(-fe, 0.0);
                    double aw = Dw + std::max( fw, 0.0);
                    double an = Dn + std::max(-fn, 0.0);
                    double as = Ds + std::max( fs, 0.0);
                    double ab = Db + std::max(-fb, 0.0);
                    double af = Df + std::max( ff, 0.0);
                    double ap = ae + aw + an + as + ab + af;
                    double b = facexy * (p_temp[indexP(i, j, k-1)] - p_temp[indexP(i, j, k)]);
                    double old = w_star[indexW(i, j, k)];                   
                    w_star[indexW(i, j, k)] = (ae * w_star[indexW(i+1, j, k)] + aw * w_star[indexW(i-1, j, k)] +
                                            an * w_star[indexW(i, j+1, k)] + as * w_star[indexW(i, j-1, k)] +
                                            af * w_star[indexW(i, j, k-1)] + ab * w_star[indexW(i, j, k+1)] + b) / ap;
                    w_star[indexW(i, j, k)] = UNDER_RELAX * w_star[indexW(i, j, k)] + (1 - UNDER_RELAX) * old;
                    wp[indexW(i, j, k)] = ap;
                    residual += std::pow(w_star[indexW(i, j, k)] - old, 2);
                    iter++;
                }
            }
        }
        
        for (auto &e: sendData) {
            e.second.resize(0);
        }
        for (auto &e: recvData) {
            e.second.resize(0);
        }               
        for (auto &v: w_boundary) {
            for (auto &e: v.second) {
                sendData[v.first].push_back(w_star[e]);
            }
            recvData[v.first].resize(recvData[v.first].size() + sendData[v.first].size()); 
        }    
        for (auto &e: w_boundary) {
            MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 6, 
                         recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 6, 
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (auto &v: w_halo) {
            for (int e = 0; e < v.second.size(); ++e) {
                w_star[w_halo[v.first][e]] = recvData[v.first][e];
            }
        }        
        step++;
        residual = residual / (zEnd - zStart + 1) / (yEnd - yStart + 1) / (xEnd - xStart + 1);
        MPI_Allreduce(&residual, &global_residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_residual = std::sqrt(global_residual);
        applyZBoundaryCondition(w_star);
    }   
}

void
FluidSolver::solveYmomentum(std::vector<double> &u_star, std::vector<double> &v_star, std::vector<double> &w_star,
    std::vector<double>& p_temp, std::vector<double>& vp, const I64 &xStart, const I64 &xEnd, const I64 &yStart,
                            const I64 &yEnd, const I64 &zStart, const I64 &zEnd) {
    
    double residual = 1, global_residual = residual;
    int step = 0, left_step = step;
    int iter = 0;
    while (step < GS_STEP) {
        if (global_residual < CONVERGENCE) {
            // left_step = GS_STEP - step;
            break;
        } 
        iter = 0;
        for (int k = zStart; k <= zEnd; ++k) {
            for (int j = yStart; j <= yEnd-1; ++j) {
                for (int i = xStart; i <= xEnd; ++i) {
                    if (p[getC(i-1, j-1, k-1)].isEmpty) {
                        continue;
                    }
                    double facexy = v_param[iter*9];
                    double facexz = v_param[iter*9+1];
                    double faceyz = v_param[iter*9+2];

                    double fn = 0.5 * RHO * facexz * (v_star[indexV(i, j, k)] + v_star[indexV(i, j - 1, k)]);
                    double fs = 0.5 * RHO * facexz * (v_star[indexV(i, j, k)] + v_star[indexV(i, j + 1, k)]);
                    double fe = 0.5 * RHO * faceyz * (u_star[indexU(i, j, k)] + u_star[indexU(i, j - 1, k)]);
                    double fw = 0.5 * RHO * faceyz * (u_star[indexU(i + 1, j, k)] + u_star[indexU(i + 1, j - 1, k)]);
                    double fb = 0.5 * RHO * facexy * (w_star[indexW(i, j, k + 1)] + w_star[indexW(i, j - 1, k + 1)]);
                    double ff = 0.5 * RHO * facexy * (w_star[indexW(i, j, k)] + w_star[indexW(i, j - 1, k)]);

                    double Dw = v_param[iter*9+3], De = v_param[iter*9+4];
                    double Ds = v_param[iter*9+5], Dn = v_param[iter*9+6];
                    double Db = v_param[iter*9+7], Df = v_param[iter*9+8];

                    double ae = De + std::max(-fe, 0.0);
                    double aw = Dw + std::max( fw, 0.0);
                    double an = Dn + std::max(-fn, 0.0);
                    double as = Ds + std::max( fs, 0.0);
                    double ab = Db + std::max(-fb, 0.0);
                    double af = Df + std::max( ff, 0.0);
                    double ap = ae + aw + an + as + ab + af;
                    double b = facexz * (p_temp[indexP(i, j-1, k)] - p_temp[indexP(i, j, k)]);
                    double old = v_star[indexV(i, j, k)];                    
                    v_star[indexV(i, j, k)] = (ae * v_star[indexV(i+1, j, k)] + aw * v_star[indexV(i-1, j, k)] +
                                            an * v_star[indexV(i, j+1, k)] + as * v_star[indexV(i, j-1, k)] +
                                            af * v_star[indexV(i, j, k-1)] + ab * v_star[indexV(i, j, k+1)] + b) / ap;
                    v_star[indexV(i, j, k)] = UNDER_RELAX * v_star[indexV(i, j, k)] + (1 - UNDER_RELAX) * old;
                    vp[indexV(i, j, k)] = ap;
                    residual += std::pow(v_star[indexV(i, j, k)] - old, 2);
                    iter++;
                }
            }
        }

        for (auto &e: sendData) {
            e.second.resize(0);
        }
        for (auto &e: recvData) {
            e.second.resize(0);
        }

        for (auto &v: v_boundary) {
            for (auto &e: v.second) {
                sendData[v.first].push_back(v_star[e]);
            }
            recvData[v.first].resize(recvData[v.first].size() + sendData[v.first].size()); 
        } 
        for (auto &e: v_boundary) {
            MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 5, 
                         recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 5, 
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for (auto &v: v_halo) {
            for (int e = 0; e < v.second.size(); ++e) {
                v_star[v_halo[v.first][e]] = recvData[v.first][e];
            }
        }         
        step++;
        residual = residual / (zEnd - zStart + 1) / (yEnd - yStart + 1) / (xEnd - xStart + 1);
        MPI_Allreduce(&residual, &global_residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_residual = std::sqrt(global_residual);
        applyYBoundaryCondition(v_star);
    } 
}

void FluidSolver::solveXmomentum(std::vector<double> &u_star, std::vector<double> &v_star, std::vector<double> &w_star,
    std::vector<double>& p_temp, std::vector<double>& up, const I64 &xStart, const I64 &xEnd, const I64 &yStart,
                            const I64 &yEnd, const I64 &zStart, const I64 &zEnd) {
    double residual = 1, global_residual = residual;
    int step = 0, left_step = step;
    int iter = 0;
    while (step < GS_STEP) {
        if (global_residual < CONVERGENCE) {
            break;
        } 
        iter = 0;
        for (int k = zStart; k <= zEnd; ++k) {
            for (int j = yStart; j <= yEnd; ++j) {
                for (int i = xStart; i <= xEnd-1; ++i) {
                    if (p[getC(i-1, j-1, k-1)].isEmpty) {
                        continue;
                    }
                    double facexy = u_param[iter*9];
                    double facexz = u_param[iter*9+1];
                    double faceyz = u_param[iter*9+2];

                    double fw = 0.5 * RHO * faceyz * (u_star[indexU(i, j, k)] + u_star[indexU(i - 1, j, k)]);
                    double fe = 0.5 * RHO * faceyz * (u_star[indexU(i, j, k)] + u_star[indexU(i + 1, j, k)]);
                    double fs = 0.5 * RHO * facexz * (v_star[indexV(i, j, k)] + v_star[indexV(i - 1, j, k)]);
                    double fn = 0.5 * RHO * facexz * (v_star[indexV(i, j + 1, k)] + v_star[indexV(i - 1, j + 1, k)]);
                    double fb = 0.5 * RHO * facexy * (w_star[indexW(i, j, k + 1)] + w_star[indexW(i - 1, j, k + 1)]);
                    double ff = 0.5 * RHO * facexy * (w_star[indexW(i, j, k)] + w_star[indexW(i - 1, j, k)]);

                    double Dw = u_param[iter*9+3], De = u_param[iter*9+4];
                    double Ds = u_param[iter*9+5], Dn = u_param[iter*9+6];
                    double Db = u_param[iter*9+7], Df = u_param[iter*9+8];

                    double ae = De + std::max(-fe, 0.0);
                    double aw = Dw + std::max( fw, 0.0);
                    double an = Dn + std::max(-fn, 0.0);
                    double as = Ds + std::max( fs, 0.0);
                    double ab = Db + std::max(-fb, 0.0);
                    double af = Df + std::max( ff, 0.0);
                    double ap = ae + aw + an + as + ab + af;
                    double b = faceyz * (p_temp[indexP((i-1), j, k)] - p_temp[indexP(i, j, k)]);
                    double old = u_star[indexU(i, j, k)];
                    u_star[indexU(i, j, k)] = (ae * u_star[indexU(i+1, j, k)] + aw * u_star[indexU(i-1, j, k)] +
                                            an * u_star[indexU(i, j+1, k)] + as * u_star[indexU(i, j-1, k)] +
                                            af * u_star[indexU(i, j, k-1)] + ab * u_star[indexU(i, j, k+1)] + b) / ap;
                    u_star[indexU(i, j, k)] = UNDER_RELAX * u_star[indexU(i, j, k)] + (1 - UNDER_RELAX) * old;              
                    up[indexU(i, j, k)] = ap;
                    residual += std::pow(u_star[indexU(i, j, k)] - old, 2);
                    iter++;
                }
            }
        }
        for (auto &e: sendData) {
            e.second.resize(0);
        }
        for (auto &e: recvData) {
            e.second.resize(0);
        }        
        for (auto &v: u_boundary) {
            for (auto &e: v.second) {
                sendData[v.first].push_back(u_star[e]);
            }
            recvData[v.first].resize(recvData[v.first].size() + sendData[v.first].size());        
        } 
                         
        for (auto &e: u_boundary) {
            MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 0, 
                         recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 0, 
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (auto &v: u_halo) {
            for (int e = 0; e < v.second.size(); ++e) {
                u_star[u_halo[v.first][e]] = recvData[v.first][e];
            }
        }
        step++;
        residual = residual / (zEnd - zStart + 1) / (yEnd - yStart + 1) / (xEnd - xStart + 1);
        MPI_Allreduce(&residual, &global_residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_residual = std::sqrt(global_residual);
        applyXBoundaryCondition(u_star);

    }
}

void FluidSolver::solveCorrectionPressure(const std::vector<double> &u_star, const std::vector<double> &v_star,
                                     const std::vector<double> &w_star, std::vector<double> &p_prime,
                                     std::vector<double> & up, std::vector<double> & vp, std::vector<double> & wp,
                                     const I64 &xStart, const I64 &xEnd, 
                                     const I64 &yStart, const I64 &yEnd, 
                                     const I64 &zStart, const I64 &zEnd) {
    double residual2 = 1;
    clock_t start = clock();
    int iter = 0, left_step = iter;   
    for (; iter < GS_STEP; ++iter) {
        if (residual2 < CONVERGENCE2) {
            left_step = GS_STEP - iter;
            break;
        }
        for (int k = zStart; k <= zEnd-1; ++k) {
            for (int j = yStart; j <= yEnd-1; ++j) {
                for (int i = xStart; i <= xEnd-1; ++i) {
                    if (p[getC(i-1, j-1, k-1)].isEmpty) {
                        continue;
                    }
                    double ap_u_approx = up[indexU(i, j, k)];
                    double ap_v_approx = vp[indexV(i, j, k)];
                    double ap_w_approx = wp[indexW(i, j, k)];

                    // 防止除以零
                    if (ap_u_approx == 0) ap_u_approx = 1e-10;
                    if (ap_v_approx == 0) ap_v_approx = 1e-10;
                    if (ap_w_approx == 0) ap_w_approx = 1e-10;
                    double xp = (i+1-1+metadata[6]-1) < ulocation.size()? (ulocation[i+1-1+metadata[6]-1]-ulocation[i-1+metadata[6]-1])
                                                                        : (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]);
                    double xw = (i-1-1+metadata[6]-1) >= 0? (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]): xp;
                    double xe = (i+2-1+metadata[6]-1) < ulocation.size()? (ulocation[i+2-1+metadata[6]-1]-ulocation[i+1-1+metadata[6]-1]): xp;

                    double yp = (j+1-1+metadata[7]-1) < vlocation.size()? (vlocation[j+1-1+metadata[7]-1]-vlocation[j-1+metadata[7]-1])
                                                                        : (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]);
                    double ys = (j-1-1+metadata[7]-1) >= 0? (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]): yp;
                    double yn = (j+2-1+metadata[7]-1) < vlocation.size()? (vlocation[j+2-1+metadata[7]-1]-vlocation[j+1-1+metadata[7]-1]): yp;

                    double zp = (k+1-1+metadata[8]-1) < wlocation.size()? (wlocation[k+1-1+metadata[8]-1]-wlocation[k-1+metadata[8]-1])
                                                                        : (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]);
                    double zf = (k-1-1+metadata[8]-1) >= 0? (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]): zp;
                    double zb = (k+2-1+metadata[8]-1) < wlocation.size()? (wlocation[k+2-1+metadata[8]-1]-wlocation[k+1-1+metadata[8]-1]): zp;

                    double facexy = xp * yp;
                    double facexz = xp * zp;
                    double faceyz = yp * zp;

                    double a_e_prime = RHO * (faceyz) / ap_u_approx;
                    double a_w_prime = RHO * (faceyz) / ap_u_approx;
                    double a_n_prime = RHO * (facexz) / ap_v_approx;
                    double a_s_prime = RHO * (facexz) / ap_v_approx;
                    double a_f_prime = RHO * (facexy) / ap_w_approx;
                    double a_b_prime = RHO * (facexy) / ap_w_approx;

                    // 连续性方程的源项 (质量不平衡)
                    double b_cont = RHO * ((-u_star[indexU(i + 1, j, k)] + u_star[indexU(i, j, k)]) +
                                           (-v_star[indexV(i, j + 1, k)] + v_star[indexV(i, j, k)]) +
                                           (-w_star[indexW(i, j, k + 1)] + w_star[indexW(i, j, k)]));

                    double ap_prime_val = 0;
                    double p_prime_neighbors_sum = 0.0;
                    p_prime_neighbors_sum = a_e_prime * p_prime[indexP(i+1,j,k)] + a_w_prime * p_prime[indexP(i-1,j,k)]
                                        + a_n_prime * p_prime[indexP(i,j+1,k)] + a_s_prime * p_prime[indexP(i,j-1,k)]
                                        + a_f_prime * p_prime[indexP(i,j,k-1)] + a_b_prime * p_prime[indexP(i,j,k+1)];
                    ap_prime_val = a_e_prime + a_w_prime + a_n_prime + a_s_prime + a_f_prime + a_b_prime;
                    double old_p_prime_val = p_prime[indexP(i, j, k)];
                    if (ap_prime_val == 0) ap_prime_val = 1e-10; 

                    p_prime[indexP(i, j, k)] = (p_prime_neighbors_sum + b_cont) / ap_prime_val;
                    residual2 += std::pow(p_prime[indexP(i, j, k)] - old_p_prime_val, 2);
                }
            }
        }
        for (auto &e: sendData) {
            e.second.resize(0);
        }
        for (auto &e: recvData) {
            e.second.resize(0);
        }

        for (auto &v: p_boundary) {
            for (auto &e: v.second) {
                sendData[v.first].push_back(p_prime[e]);
            }
            recvData[v.first].resize(recvData[v.first].size() + sendData[v.first].size()); 
        }
        for (auto &e: p_boundary) {
            MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 1,
                         recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 1,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (auto &v: p_halo) {
            for (int e = 0; e < v.second.size(); ++e) {
                p_prime[p_halo[v.first][e]] = recvData[v.first][e];
            }
        }                      
        // p_prime[indexP(0, 0, 0)] = 0;
        residual2 = std::sqrt(residual2 / (zEnd - zStart + 1) / (yEnd - yStart + 1) / (xEnd - xStart + 1));
    }
    while (left_step > 0) {
        for (auto &e: sendData) {
            e.second.resize(0);
        }
        for (auto &e: recvData) {
            e.second.resize(0);
        }

        for (auto &v: p_boundary) {
            for (auto &e: v.second) {
                sendData[v.first].push_back(p_prime[e]);
            }
            recvData[v.first].resize(recvData[v.first].size() + sendData[v.first].size()); 
        }
        for (auto &e: p_boundary) {
            MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 1,
                         recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 1,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (auto &v: p_halo) {
            for (int e = 0; e < v.second.size(); ++e) {
                p_prime[p_halo[v.first][e]] = recvData[v.first][e];
            }
        }         
        left_step--;
    }

    clock_t end = clock();  
    double duration = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "CPU耗时: " << duration << " 秒 " << iter << std::endl;        
}

void FluidSolver::solveTemp(const std::vector<double> &u_star, const std::vector<double> &v_star,
                                     const std::vector<double> &w_star, std::vector<double> &t_temp,
                                     const I64 &xStart, const I64 &xEnd, 
                                     const I64 &yStart, const I64 &yEnd, 
                                     const I64 &zStart, const I64 &zEnd) {
    double residual = 1, global_residual = residual;
    int step = 0, left_step = step;
    int iter = 0;   
    while (step < GS_STEP) {
        if (global_residual < CONVERGENCE) {
            break;
        } 
        iter = 0;
        for (int k = zStart; k <= zEnd; ++k) {
            for (int j = yStart; j <= yEnd; ++j) {
                for (int i = xStart; i <= xEnd; ++i) {
                    if (p[getC(i-1, j-1, k-1)].isEmpty) {
                        continue;
                    }
                    double xp = (i+1-1+metadata[6]-1) < ulocation.size()? (ulocation[i+1-1+metadata[6]-1]-ulocation[i-1+metadata[6]-1])
                                                                        : (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]);
                    double xw = (i-1-1+metadata[6]-1) >= 0? (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]): xp;
                    double xe = (i+2-1+metadata[6]-1) < ulocation.size()? (ulocation[i+2-1+metadata[6]-1]-ulocation[i+1-1+metadata[6]-1]): xp;

                    double yp = (j+1-1+metadata[7]-1) < vlocation.size()? (vlocation[j+1-1+metadata[7]-1]-vlocation[j-1+metadata[7]-1])
                                                                        : (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]);
                    double ys = (j-1-1+metadata[7]-1) >= 0? (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]): yp;
                    double yn = (j+2-1+metadata[7]-1) < vlocation.size()? (vlocation[j+2-1+metadata[7]-1]-vlocation[j+1-1+metadata[7]-1]): yp;

                    double zp = (k+1-1+metadata[8]-1) < wlocation.size()? (wlocation[k+1-1+metadata[8]-1]-wlocation[k-1+metadata[8]-1])
                                                                        : (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]);
                    double zf = (k-1-1+metadata[8]-1) >= 0? (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]): zp;
                    double zb = (k+2-1+metadata[8]-1) < wlocation.size()? (wlocation[k+2-1+metadata[8]-1]-wlocation[k+1-1+metadata[8]-1]): zp;

                    double facexy = xp * yp;
                    double facexz = xp * zp;
                    double faceyz = yp * zp;
                    double spacingx = (xp+xw)/2.0;
                    double spacingy = (yp+ys)/2.0;
                    double spacingz = (zp+zf)/2.0;
                    double conv = -u_star[indexU(i,j,k)]*(t_temp[indexP(i+1,j,k)]-t_temp[indexP(i-1,j,k)])
                                   -v_star[indexV(i,j,k)]*(t_temp[indexP(i,j+1,k)]-t_temp[indexP(i,j-1,k)])
                                   -w_star[indexW(i,j,k)]*(t_temp[indexP(i,j,k+1)]-t_temp[indexP(i,j,k-1)]);

                    double diff = ALPHA * (t_temp[indexP(i+1,j,k)]+t_temp[indexP(i-1,j,k)]) / spacingx / spacingx
                                + ALPHA * (t_temp[indexP(i,j+1,k)]+t_temp[indexP(i,j-1,k)]) / spacingy / spacingy
                                + ALPHA * (t_temp[indexP(i,j,k+1)]+t_temp[indexP(i,j,k-1)]) / spacingz / spacingz;

                    double denom = 2 * ALPHA * (1/spacingx/spacingx+1/spacingy/spacingy+1/spacingz/spacingz);
                    double old = t_temp[indexP(i,j,k)];
                    t_temp[indexP(i,j,k)] = (conv + diff) / denom * (UNDER_RELAX) + (1-UNDER_RELAX) * old;
                    residual += std::pow(t_temp[indexP(i,j,k)] - old, 2);
                    iter++;
                }
            }
        }

        for (auto &e: sendData) {
            e.second.resize(0);
        }
        for (auto &e: recvData) {
            e.second.resize(0);
        }   
        // 先根据p_boundary通信流体温度
        for (auto &v: p_boundary) {
            for (auto &e: v.second) {
                sendData[v.first].push_back(t_temp[e]);
            }
            // std::cout << std::endl;
            recvData[v.first].resize(recvData[v.first].size() + sendData[v.first].size()); 
        }
        for (auto &e: p_boundary) {
            MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 2,
                         recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 2,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (auto &v: p_halo) {
            for (int e = 0; e < v.second.size(); ++e) {
                t_temp[p_halo[v.first][e]] = recvData[v.first][e];               
            }
        }
        step++;
        residual = residual / (zEnd - zStart + 1) / (yEnd - yStart + 1) / (xEnd - xStart + 1);
        MPI_Allreduce(&residual, &global_residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_residual = std::sqrt(global_residual);                 
    }
}
void FluidSolver::init() {
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];  
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];        

    // applyOuterBoundaryCondition(u_star, p_temp);
    applyXBoundaryCondition(u_star);
    applyYBoundaryCondition(v_star);
    applyZBoundaryCondition(w_star);
    applyCellBoundaryCondition(p_temp);
    // std::cout << metadata[6] - 1 << " " << metadata[6] - 1 + fx - 1 << " " << metadata[7] - 1 << " " << metadata[7] - 1 + fy - 1
    //             << " " << metadata[8] - 1 << " " <<  metadata[8] - 1 + fz - 1 << std::endl;

    ulocation = mesh.cycasMesh.xPoints;
    vlocation = mesh.cycasMesh.yPoints;
    wlocation = mesh.cycasMesh.zPoints;

    // std::cout << std::endl;    
    // std::cout << "init done. \n";    
    for (int k = 0; k < cz; ++k) {
        for (int j = 0; j < cy; ++j) {
            for (int i = 0; i < fx; ++i) {
                int rank;
                auto pid = rank;
                if (ux[getX(i, j, k)].isBoundary 
                    && parseString("pr_", ux[getX(i, j, k)].info->name, rank, pid)) {
                    // std::cout << rank << " " << pid << std::endl;
                    if (u_boundary.find(pid) == u_boundary.end()) 
                        u_boundary[pid] = std::vector<I64>();
                    if (v_boundary.find(pid) == v_boundary.end())
                        v_boundary[pid] = std::vector<I64>(); 
                    if (w_boundary.find(pid) == w_boundary.end())
                        w_boundary[pid] = std::vector<I64>();
                    if (p_boundary.find(pid) == p_boundary.end())
                        p_boundary[pid] = std::vector<I64>(); 
                    if (u_halo.find(pid) == u_halo.end())
                        u_halo[pid] = std::vector<I64>();
                    if (v_halo.find(pid) == v_halo.end())
                        v_halo[pid] = std::vector<I64>();
                    if (w_halo.find(pid) == w_halo.end())
                        w_halo[pid] = std::vector<I64>();
                    if (p_halo.find(pid) == p_halo.end())
                        p_halo[pid] = std::vector<I64>();

                    if (i + 1 >= fx) {
                        u_boundary[pid].push_back(indexU((i-1)+1, j+1, k+1));
                        p_boundary[pid].push_back(indexP((i-1)+1, j+1, k+1));
                        v_boundary[pid].push_back(indexV((i-1)+1, j+1, k+1));
                        v_boundary[pid].push_back(indexV((i-1)+1, j+1+1, k+1));
                        w_boundary[pid].push_back(indexW((i-1)+1, j+1, k+1));
                        w_boundary[pid].push_back(indexW((i-1)+1, j+1, k+1+1));

                        u_halo[pid].push_back(indexU((i)+1, j+1, k+1));
                        p_halo[pid].push_back(indexP((i)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, j+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, j+1, k+1+1));

                    } else if (i - 1 < 0) {
                        u_boundary[pid].push_back(indexU((i)+1, j+1, k+1));
                        p_boundary[pid].push_back(indexP((i)+1, j+1, k+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1, k+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1+1, k+1));
                        w_boundary[pid].push_back(indexW((i)+1, j+1, k+1));
                        w_boundary[pid].push_back(indexW((i)+1, j+1, k+1+1)); 

                        u_halo[pid].push_back(indexU((i-1)+1, j+1, k+1));
                        p_halo[pid].push_back(indexP((i-1)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i-1)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i-1)+1, j+1+1, k+1));
                        w_halo[pid].push_back(indexW((i-1)+1, j+1, k+1));
                        w_halo[pid].push_back(indexW((i-1)+1, j+1, k+1+1));                                                                       
                    } else if (i <= cx - 1 && p[getC(i, j, k)].isEmpty) {
                        u_boundary[pid].push_back(indexU((i-1)+1, j+1, k+1));
                        p_boundary[pid].push_back(indexP((i-1)+1, j+1, k+1));
                        v_boundary[pid].push_back(indexV((i-1)+1, j+1, k+1));
                        v_boundary[pid].push_back(indexV((i-1)+1, j+1+1, k+1));
                        w_boundary[pid].push_back(indexW((i-1)+1, j+1, k+1));
                        w_boundary[pid].push_back(indexW((i-1)+1, j+1, k+1+1));

                        u_halo[pid].push_back(indexU((i)+1, j+1, k+1));
                        p_halo[pid].push_back(indexP((i)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, j+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, j+1, k+1+1));
                    } else if (i >= 1 && p[getC((i-1), (j), (k))].isEmpty) {
                        u_boundary[pid].push_back(indexU((i)+1, j+1, k+1));
                        p_boundary[pid].push_back(indexP((i)+1, j+1, k+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1, k+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1+1, k+1));
                        w_boundary[pid].push_back(indexW((i)+1, j+1, k+1));
                        w_boundary[pid].push_back(indexW((i)+1, j+1, k+1+1)); 

                        u_halo[pid].push_back(indexU((i-1)+1, j+1, k+1));
                        p_halo[pid].push_back(indexP((i-1)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i-1)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i-1)+1, j+1+1, k+1));
                        w_halo[pid].push_back(indexW((i-1)+1, j+1, k+1));
                        w_halo[pid].push_back(indexW((i-1)+1, j+1, k+1+1));                          
                    }
                } else if (ux[getX(i, j, k)].isBoundary 
                    && ux[getX(i, j, k)].info != nullptr && parseString("ib_", ux[getX(i, j, k)].info->name, rank, pid)) {
                    // if (rank != pid) { 
                        if (t_boundary.find(pid) == t_boundary.end()) 
                            t_boundary[pid] = std::vector<I64>();
                        if (t_halo.find(pid) == t_halo.end())
                            t_halo[pid] = std::vector<I64>();                        
                        if (i + 1 >= fx) {
                            t_boundary[pid].push_back(indexP(i-1+1, j+1, k+1));
                            t_halo[pid].push_back(indexP(i+1, j+1, k+1));
                        } else if (i - 1 < 0) {
                            t_boundary[pid].push_back(indexP(i+1, j+1, k+1));
                            t_halo[pid].push_back(indexP(i-1+1, j+1, k+1));
                        } 
                        else if (i <= cx - 1 && p[getC(i, j, k)].isEmpty) {
                            t_boundary[pid].push_back(indexP(i-1+1, j+1, k+1));
                            t_halo[pid].push_back(indexP(i+1, j+1, k+1));
                        } else if (i >= 1 && p[getC((i-1), (j), (k))].isEmpty) {
                            t_boundary[pid].push_back(indexP(i+1, j+1, k+1));
                            t_halo[pid].push_back(indexP(i-1+1, j+1, k+1));
                        }
                    // }
                    ux[getX(i, j, k)].u = u_star[indexU(i+1,j+1,k+1)] = 0;
                }
            }
        }
    }    

    for (int k = 0; k < cz; ++k) {
        for (int j = 0; j < fy; ++j) {
            for (int i = 0; i < cx; ++i) {
                int rank;
                auto pid = rank;                
                if (uy[getY(i, j, k)].isBoundary && uy[getY(i, j, k)].info != nullptr 
                    && parseString("pr_", uy[getY(i, j, k)].info->name, rank, pid)) {
                    if (u_boundary.find(pid) == u_boundary.end())
                        u_boundary[pid] = std::vector<I64>();
                    if (v_boundary.find(pid) == v_boundary.end())
                        v_boundary[pid] = std::vector<I64>(); 
                    if (w_boundary.find(pid) == w_boundary.end())
                        w_boundary[pid] = std::vector<I64>();
                    if (p_boundary.find(pid) == p_boundary.end())
                        p_boundary[pid] = std::vector<I64>(); 
                    if (u_halo.find(pid) == u_halo.end())
                        u_halo[pid] = std::vector<I64>();
                    if (v_halo.find(pid) == v_halo.end())
                        v_halo[pid] = std::vector<I64>();
                    if (w_halo.find(pid) == w_halo.end())
                        w_halo[pid] = std::vector<I64>();
                    if (p_halo.find(pid) == p_halo.end())
                        p_halo[pid] = std::vector<I64>();

                    if (j + 1 >= fy) {
                        v_boundary[pid].push_back(indexV(i+1, (j-1)+1, k+1));
                        p_boundary[pid].push_back(indexP(i+1, (j-1)+1, k+1));
                        u_boundary[pid].push_back(indexU(i+1, (j-1)+1, k+1));
                        u_boundary[pid].push_back(indexU(i+1+1, (j-1)+1, k+1));
                        w_boundary[pid].push_back(indexW(i+1, (j-1)+1, k+1));
                        w_boundary[pid].push_back(indexW(i+1, (j-1)+1, k+1+1));

                        v_halo[pid].push_back(indexV((i)+1, j+1, k+1));
                        p_halo[pid].push_back(indexP((i)+1, j+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1, j+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1+1, j+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, j+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, j+1, k+1+1));

                    } else if (j - 1 < 0) {
                        v_boundary[pid].push_back(indexV(i+1, (j)+1, k+1));
                        p_boundary[pid].push_back(indexP(i+1, (j)+1, k+1));
                        u_boundary[pid].push_back(indexU(i+1, (j)+1, k+1));
                        u_boundary[pid].push_back(indexU(i+1+1, (j)+1, k+1));
                        w_boundary[pid].push_back(indexW(i+1, (j)+1, k+1));
                        w_boundary[pid].push_back(indexW(i+1, (j)+1, k+1+1));

                        v_halo[pid].push_back(indexV((i)+1, (j-1)+1, k+1));
                        p_halo[pid].push_back(indexP((i)+1, (j-1)+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1, (j-1)+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1+1, (j-1)+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, (j-1)+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, (j-1)+1, k+1+1));  

                    } else if (j <= cy - 1 && p[getC(i, j, k)].isEmpty) {
                        v_boundary[pid].push_back(indexV(i+1, (j-1)+1, k+1));
                        p_boundary[pid].push_back(indexP(i+1, (j-1)+1, k+1));
                        u_boundary[pid].push_back(indexU(i+1, (j-1)+1, k+1));
                        u_boundary[pid].push_back(indexU(i+1+1, (j-1)+1, k+1));
                        w_boundary[pid].push_back(indexW(i+1, (j-1)+1, k+1));
                        w_boundary[pid].push_back(indexW(i+1, (j-1)+1, k+1+1));

                        v_halo[pid].push_back(indexV((i)+1, j+1, k+1));
                        p_halo[pid].push_back(indexP((i)+1, j+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1, j+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1+1, j+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, j+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, j+1, k+1+1));

                    } else if (j >= 1 && p[getC((i), (j-1), (k))].isEmpty) {
                        v_boundary[pid].push_back(indexV(i+1, (j)+1, k+1));
                        p_boundary[pid].push_back(indexP(i+1, (j)+1, k+1));
                        u_boundary[pid].push_back(indexU(i+1, (j)+1, k+1));
                        u_boundary[pid].push_back(indexU(i+1+1, (j)+1, k+1));
                        w_boundary[pid].push_back(indexW(i+1, (j)+1, k+1));
                        w_boundary[pid].push_back(indexW(i+1, (j)+1, k+1+1));

                        v_halo[pid].push_back(indexV((i)+1, (j-1)+1, k+1));
                        p_halo[pid].push_back(indexP((i)+1, (j-1)+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1, (j-1)+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1+1, (j-1)+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, (j-1)+1, k+1));
                        w_halo[pid].push_back(indexW((i)+1, (j-1)+1, k+1+1));
                    }
                } else if (uy[getY(i, j, k)].isBoundary && uy[getY(i, j, k)].info != nullptr 
                        && parseString("ib_", uy[getY(i, j, k)].info->name, rank, pid)) {
                    // if (rank != pid) { 
                        if (t_boundary.find(pid) == t_boundary.end()) 
                            t_boundary[pid] = std::vector<I64>();
                        if (t_halo.find(pid) == t_halo.end())
                            t_halo[pid] = std::vector<I64>();                         
                        if (j + 1 >= fy) {
                            t_boundary[pid].push_back(indexP(i+1, j-1+1, k+1));
                            t_halo[pid].push_back(indexP(i+1, j+1, k+1));
                        } else if (j - 1 < 0) {
                            t_boundary[pid].push_back(indexP(i+1, j+1, k+1));
                            t_halo[pid].push_back(indexP(i+1, j-1+1, k+1));
                        } else if (j <= cy - 1 && p[getC(i, j, k)].isEmpty) {
                            t_boundary[pid].push_back(indexP(i+1, j-1+1, k+1));
                            t_halo[pid].push_back(indexP(i+1, j+1, k+1));
                        } else if (j >= 1 && p[getC((i), (j-1), (k))].isEmpty) {
                            t_boundary[pid].push_back(indexP(i+1, j+1, k+1));
                            t_halo[pid].push_back(indexP(i+1, j-1+1, k+1));
                        }
                    // }
                    uy[getY(i, j, k)].u = v_star[indexV(i+1,j+1,k+1)] = 0;                    
                }
            }
        }
    } 

    for (int k = 0; k < fz; ++k) {
        for (int j = 0; j < cy; ++j) {
            for (int i = 0; i < cx; ++i) {
                int rank;
                auto pid = rank;                 
                if (uz[getZ(i, j, k)].isBoundary && uz[getZ(i, j, k)].info != nullptr && parseString("pr_", uz[getZ(i, j, k)].info->name, rank, pid)) {
                    if (u_boundary.find(pid) == u_boundary.end())
                        u_boundary[pid] = std::vector<I64>();
                    if (v_boundary.find(pid) == v_boundary.end())
                        v_boundary[pid] = std::vector<I64>(); 
                    if (w_boundary.find(pid) == w_boundary.end())
                        w_boundary[pid] = std::vector<I64>();
                    if (p_boundary.find(pid) == p_boundary.end())
                        p_boundary[pid] = std::vector<I64>(); 
                    if (u_halo.find(pid) == u_halo.end())
                        u_halo[pid] = std::vector<I64>();
                    if (v_halo.find(pid) == v_halo.end())
                        v_halo[pid] = std::vector<I64>();
                    if (w_halo.find(pid) == w_halo.end())
                        w_halo[pid] = std::vector<I64>();
                    if (p_halo.find(pid) == p_halo.end())
                        p_halo[pid] = std::vector<I64>();

                    if (k + 1 >= fz) {
                        w_boundary[pid].push_back(indexW((i)+1, j+1, (k-1)+1));
                        p_boundary[pid].push_back(indexP((i)+1, j+1, (k-1)+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1, (k-1)+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1+1, (k-1)+1));
                        u_boundary[pid].push_back(indexU((i)+1, j+1, (k-1)+1));
                        u_boundary[pid].push_back(indexU((i)+1+1, j+1, (k-1)+1));

                        w_halo[pid].push_back(indexW((i)+1, j+1, k+1));
                        p_halo[pid].push_back(indexP((i)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1, j+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1+1, j+1, k+1));

                    } else if (k < 1) {
                        w_boundary[pid].push_back(indexW((i)+1, j+1, (k)+1));
                        p_boundary[pid].push_back(indexP((i)+1, j+1, (k)+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1, (k)+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1+1, (k)+1));
                        u_boundary[pid].push_back(indexU((i)+1, j+1, (k)+1));
                        u_boundary[pid].push_back(indexU((i)+1+1, j+1, (k)+1));   

                        w_halo[pid].push_back(indexW((i)+1, j+1, (k-1)+1));
                        p_halo[pid].push_back(indexP((i)+1, j+1, (k-1)+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1, (k-1)+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1+1, (k-1)+1));
                        u_halo[pid].push_back(indexU((i)+1, j+1, (k-1)+1));
                        u_halo[pid].push_back(indexU((i)+1+1, j+1, (k-1)+1));                                                                    
                    } else if (k <= cz - 1 && p[getC(i, j, k)].isEmpty) {

                        w_boundary[pid].push_back(indexW((i)+1, j+1, (k-1)+1));
                        p_boundary[pid].push_back(indexP((i)+1, j+1, (k-1)+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1, (k-1)+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1+1, (k-1)+1));
                        u_boundary[pid].push_back(indexU((i)+1, j+1, (k-1)+1));
                        u_boundary[pid].push_back(indexU((i)+1+1, j+1, (k-1)+1));

                        w_halo[pid].push_back(indexW((i)+1, j+1, k+1));
                        p_halo[pid].push_back(indexP((i)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1, k+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1, j+1, k+1));
                        u_halo[pid].push_back(indexU((i)+1+1, j+1, k+1));

                    } else if (k >= 1 && p[getC((i), (j), (k-1))].isEmpty) {
                        w_boundary[pid].push_back(indexW((i)+1, j+1, (k)+1));
                        p_boundary[pid].push_back(indexP((i)+1, j+1, (k)+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1, (k)+1));
                        v_boundary[pid].push_back(indexV((i)+1, j+1+1, (k)+1));
                        u_boundary[pid].push_back(indexU((i)+1, j+1, (k)+1));
                        u_boundary[pid].push_back(indexU((i)+1+1, j+1, (k)+1));   

                        w_halo[pid].push_back(indexW((i)+1, j+1, (k-1)+1));
                        p_halo[pid].push_back(indexP((i)+1, j+1, (k-1)+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1, (k-1)+1));
                        v_halo[pid].push_back(indexV((i)+1, j+1+1, (k-1)+1));
                        u_halo[pid].push_back(indexU((i)+1, j+1, (k-1)+1));
                        u_halo[pid].push_back(indexU((i)+1+1, j+1, (k-1)+1)); 
                    }                        
                } else if (uz[getZ(i, j, k)].isBoundary && uz[getZ(i, j, k)].info != nullptr 
                        && parseString("ib_", uz[getZ(i, j, k)].info->name, rank, pid)) {

                        if (t_boundary.find(pid) == t_boundary.end()) 
                            t_boundary[pid] = std::vector<I64>();
                        if (t_halo.find(pid) == t_halo.end())
                            t_halo[pid] = std::vector<I64>();                         
                        if (k + 1 >= fz) {
                            t_boundary[pid].push_back(indexP(i+1, j+1, k-1+1));
                            t_halo[pid].push_back(indexP(i+1, j+1, k+1));
                        } else if (k - 1 < 0) {
                            t_boundary[pid].push_back(indexP(i+1, j+1, k+1));
                            t_halo[pid].push_back(indexP(i+1, j+1, k-1+1));
                        } else if (k <= cz - 1 && p[getC(i, j, k)].isEmpty) {
                            t_boundary[pid].push_back(indexP(i+1, j+1, k-1+1));
                            t_halo[pid].push_back(indexP(i+1, j+1, k+1));
                        } else if (k >= 1 && p[getC((i), (j), (k-1))].isEmpty) {
                            t_boundary[pid].push_back(indexP(i+1, j+1, k+1));
                            t_halo[pid].push_back(indexP(i+1, j+1, k-1+1));
                        }
                    uz[getZ(i, j, k)].u = w_star[indexW(i+1,j+1,k+1)] = 0;     
                }
            }
        }
    } 
     
    for (auto &e: t_boundary) {
        sendData[e.first] = std::vector<double>();
        recvData[e.first] = std::vector<double>();
    }  
    
    initParam();
}

void FluidSolver::initParam() {
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];  
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];      
    for (int k = 1; k <= cz; ++k) {
        for (int j = 1; j <= cy; ++j) {
            for (int i = 1; i <= fx-1; ++i) {
                if (p[getC(i-1, j-1, k-1)].isEmpty) {
                    continue;
                } 
                double xp = (i+1-1+metadata[6]-1) < ulocation.size()? (ulocation[i+1-1+metadata[6]-1]-ulocation[i-1+metadata[6]-1])
                                                                    : (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]);
                double xw = (i-1-1+metadata[6]-1) >= 0? (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]): xp;
                double xe = (i+2-1+metadata[6]-1) < ulocation.size()? (ulocation[i+2-1+metadata[6]-1]-ulocation[i+1-1+metadata[6]-1]): xp;

                double yp = (j+1-1+metadata[7]-1) < vlocation.size()? (vlocation[j+1-1+metadata[7]-1]-vlocation[j-1+metadata[7]-1])
                                                                    : (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]);
                double ys = (j-1-1+metadata[7]-1) >= 0? (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]): yp;
                double yn = (j+2-1+metadata[7]-1) < vlocation.size()? (vlocation[j+2-1+metadata[7]-1]-vlocation[j+1-1+metadata[7]-1]): yp;

                double zp = (k+1-1+metadata[8]-1) < wlocation.size()? (wlocation[k+1-1+metadata[8]-1]-wlocation[k-1+metadata[8]-1])
                                                                    : (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]);
                double zf = (k-1-1+metadata[8]-1) >= 0? (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]): zp;
                double zb = (k+2-1+metadata[8]-1) < wlocation.size()? (wlocation[k+2-1+metadata[8]-1]-wlocation[k+1-1+metadata[8]-1]): zp;

                double facexy = xp * yp;
                double facexz = xp * zp;
                double faceyz = yp * zp;
                double Dw = 2 * NU * faceyz / (xp+xw), De = 2 * NU * faceyz / (xp+xe);
                double Ds = 2 * NU * facexz / (ys+yp), Dn = 2 * NU * facexz / (yn+yp);
                double Db = 2 * NU * facexy / (zb+zp), Df = 2 * NU * facexy / (zf+zp);

                u_param.push_back(facexy);
                u_param.push_back(facexz);
                u_param.push_back(faceyz);
                u_param.push_back(Dw);
                u_param.push_back(De);
                u_param.push_back(Ds);
                u_param.push_back(Dn);
                u_param.push_back(Db);
                u_param.push_back(Df);
            }
        }
    }
    for (int k = 1; k <= cz; ++k) {
        for (int j = 1; j <= fy-1; ++j) {
            for (int i = 1; i <= cx; ++i) {
                if (p[getC(i-1, j-1, k-1)].isEmpty) {
                    continue;
                }
                double xp = (i+1-1+metadata[6]-1) < ulocation.size()? (ulocation[i+1-1+metadata[6]-1]-ulocation[i-1+metadata[6]-1])
                                                                    : (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]);
                double xw = (i-1-1+metadata[6]-1) >= 0? (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]): xp;
                double xe = (i+2-1+metadata[6]-1) < ulocation.size()? (ulocation[i+2-1+metadata[6]-1]-ulocation[i+1-1+metadata[6]-1]): xp;

                double yp = (j+1-1+metadata[7]-1) < vlocation.size()? (vlocation[j+1-1+metadata[7]-1]-vlocation[j-1+metadata[7]-1])
                                                                    : (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]);
                double ys = (j-1-1+metadata[7]-1) >= 0? (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]): yp;
                double yn = (j+2-1+metadata[7]-1) < vlocation.size()? (vlocation[j+2-1+metadata[7]-1]-vlocation[j+1-1+metadata[7]-1]): yp;

                double zp = (k+1-1+metadata[8]-1) < wlocation.size()? (wlocation[k+1-1+metadata[8]-1]-wlocation[k-1+metadata[8]-1])
                                                                    : (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]);
                double zf = (k-1-1+metadata[8]-1) >= 0? (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]): zp;
                double zb = (k+2-1+metadata[8]-1) < wlocation.size()? (wlocation[k+2-1+metadata[8]-1]-wlocation[k+1-1+metadata[8]-1]): zp;

                double facexy = xp * yp;
                double facexz = xp * zp;
                double faceyz = yp * zp;

                double Dw = 2 * NU * faceyz / (xp+xw), De = 2 * NU * faceyz / (xp+xe);
                double Ds = 2 * NU * facexz / (ys+yp), Dn = 2 * NU * facexz / (yn+yp);
                double Db = 2 * NU * facexy / (zb+zp), Df = 2 * NU * facexy / (zf+zp); 
                
                v_param.push_back(facexy);
                v_param.push_back(facexz);
                v_param.push_back(faceyz);
                v_param.push_back(Dw);
                v_param.push_back(De);
                v_param.push_back(Ds);
                v_param.push_back(Dn);
                v_param.push_back(Db);
                v_param.push_back(Df);                
            }
        }
    }   
    for (int k = 1; k <= fz-1; ++k) {
        for (int j = 1; j <= cy; ++j) {
            for (int i = 1; i <= cx; ++i) {
                if (p[getC(i-1, j-1, k-1)].isEmpty) {
                    continue;
                }
                double xp = (i+1-1+metadata[6]-1) < ulocation.size()? (ulocation[i+1-1+metadata[6]-1]-ulocation[i-1+metadata[6]-1])
                                                                    : (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]);
                double xw = (i-1-1+metadata[6]-1) >= 0? (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]): xp;
                double xe = (i+2-1+metadata[6]-1) < ulocation.size()? (ulocation[i+2-1+metadata[6]-1]-ulocation[i+1-1+metadata[6]-1]): xp;

                double yp = (j+1-1+metadata[7]-1) < vlocation.size()? (vlocation[j+1-1+metadata[7]-1]-vlocation[j-1+metadata[7]-1])
                                                                    : (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]);
                double ys = (j-1-1+metadata[7]-1) >= 0? (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]): yp;
                double yn = (j+2-1+metadata[7]-1) < vlocation.size()? (vlocation[j+2-1+metadata[7]-1]-vlocation[j+1-1+metadata[7]-1]): yp;

                double zp = (k+1-1+metadata[8]-1) < wlocation.size()? (wlocation[k+1-1+metadata[8]-1]-wlocation[k-1+metadata[8]-1])
                                                                    : (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]);
                double zf = (k-1-1+metadata[8]-1) >= 0? (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]): zp;
                double zb = (k+2-1+metadata[8]-1) < wlocation.size()? (wlocation[k+2-1+metadata[8]-1]-wlocation[k+1-1+metadata[8]-1]): zp;

                double facexy = xp * yp;
                double facexz = xp * zp;
                double faceyz = yp * zp;

                double Dw = 2 * NU * faceyz / (xp+xw), De = 2 * NU * faceyz / (xp+xe);
                double Ds = 2 * NU * facexz / (ys+yp), Dn = 2 * NU * facexz / (yn+yp);
                double Db = 2 * NU * facexy / (zb+zp), Df = 2 * NU * facexy / (zf+zp); 

                w_param.push_back(facexy);
                w_param.push_back(facexz);
                w_param.push_back(faceyz);
                w_param.push_back(Dw);
                w_param.push_back(De);
                w_param.push_back(Ds);
                w_param.push_back(Dn);
                w_param.push_back(Db);
                w_param.push_back(Df);
            }
        }
    } 
    for (int k = 1; k <= cz; ++k) {
        for (int j = 1; j <= cy; ++j) {
            for (int i = 1; i <= cx; ++i) {
                if (p[getC(i-1, j-1, k-1)].isEmpty) {
                    continue;
                }
                double xp = (i+1-1+metadata[6]-1) < ulocation.size()? (ulocation[i+1-1+metadata[6]-1]-ulocation[i-1+metadata[6]-1])
                                                                    : (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]);
                double xw = (i-1-1+metadata[6]-1) >= 0? (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]): xp;
                double xe = (i+2-1+metadata[6]-1) < ulocation.size()? (ulocation[i+2-1+metadata[6]-1]-ulocation[i+1-1+metadata[6]-1]): xp;

                double yp = (j+1-1+metadata[7]-1) < vlocation.size()? (vlocation[j+1-1+metadata[7]-1]-vlocation[j-1+metadata[7]-1])
                                                                    : (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]);
                double ys = (j-1-1+metadata[7]-1) >= 0? (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]): yp;
                double yn = (j+2-1+metadata[7]-1) < vlocation.size()? (vlocation[j+2-1+metadata[7]-1]-vlocation[j+1-1+metadata[7]-1]): yp;

                double zp = (k+1-1+metadata[8]-1) < wlocation.size()? (wlocation[k+1-1+metadata[8]-1]-wlocation[k-1+metadata[8]-1])
                                                                    : (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]);
                double zf = (k-1-1+metadata[8]-1) >= 0? (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]): zp;
                double zb = (k+2-1+metadata[8]-1) < wlocation.size()? (wlocation[k+2-1+metadata[8]-1]-wlocation[k+1-1+metadata[8]-1]): zp;

                double facexy = xp * yp;
                double facexz = xp * zp;
                double faceyz = yp * zp;
                double spacingx = (xp+xw)/2.0;
                double spacingy = (yp+ys)/2.0;
                double spacingz = (zp+zf)/2.0;

                p_param.push_back(spacingx);
                p_param.push_back(spacingy);
                p_param.push_back(spacingz);

            }
        }
    }
}

void FluidSolver::solve() {
    // set global boundary
    // std::cout << "start applying. \n"; 
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];  
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];        

    for (int iter = 0; iter < STEP; ++iter){

        // internal cell 
        // local 0~cy - 1 1~cy
        // absolute 0~cy+1
         I64 xStart = 1, xEnd = fx;
         I64 yStart = 1, yEnd = cy;
         I64 zStart = 1, zEnd = cz;
        //  std::cout << "start solve x momentum" << std::endl;        
         solveXmomentum(u_star, v_star, w_star, p_temp, up, xStart, xEnd, yStart, yEnd, zStart, zEnd);

         xStart = 1, xEnd = cx;
         yStart = 1, yEnd = fy;
         zStart = 1, zEnd = cz;
        //  std::cout << "start solve y momentum" << std::endl;        
         solveYmomentum(u_star, v_star, w_star, p_temp, vp, xStart, xEnd, yStart, yEnd, zStart, zEnd);

         xStart = 1, xEnd = cx;
         yStart = 1, yEnd = cy;
         zStart = 1, zEnd = fz;
        //  std::cout << "start solve z momentum" << std::endl;        
         solveZmomentum(u_star, v_star, w_star, p_temp, wp, xStart, xEnd, yStart, yEnd, zStart, zEnd);

         xStart = 1, xEnd = cx;
         yStart = 1, yEnd = cy;
         zStart = 1, zEnd = cz;
        //  std::cout << "start solve correction" << std::endl;
        //solveCorrectionPressure(u_star, v_star, w_star, p_prime,
//                               up, vp, wp, xStart, xEnd, yStart, yEnd, zStart, zEnd);
  
        // HYPRE_SolveCorrectionPressure(u_star, v_star, w_star, p_prime,
        //                          up, vp, wp, cx, cy, cz, metadata[6], metadata[7], metadata[8]);
        // reallocation 
        auto num = cx * cy * cz;
        for (int i = 0; i < cx; ++i)
            for (int j = 0; j < cy; ++j)
                for (int k = 0; k < cz; ++k)
                    num -= p[getC(i,j,k)].isEmpty? 1: 0;

        int global_start;
        MPI_Exscan(&num, &global_start, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0) {
            global_start = 0;
        }

        HYPREIJMatrix_SolveCorrectionPressure(u_star, v_star, w_star, p_prime,
                                 up, vp, wp, cx, cy, cz, metadata[6], metadata[7], metadata[8], global_start, num);

        //std::cout << "start correction" << std::endl;
        // 修正U速度，直接用ux的坐标，u_star是该坐标+1
        xStart = 1, xEnd = fx;
        yStart = 1, yEnd = cy;
        zStart = 1, zEnd = cz;        
        for (int k = zStart; k <= zEnd; ++k) {
            for (int j = yStart; j <= yEnd; ++j) {
                for (int i = xStart; i <= xEnd-1; ++i) {
                    if (p[getC(i-1, j-1, k-1)].isEmpty) {
                        // p[getC(i-1, j-1, k-1)].data = 0;
                        continue;
                    } 
                    double xp = (i+1-1+metadata[6]-1) < ulocation.size()? (ulocation[i+1-1+metadata[6]-1]-ulocation[i-1+metadata[6]-1])
                                                                        : (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]);
                    double xw = (i-1-1+metadata[6]-1) >= 0? (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]): xp;
                    double xe = (i+2-1+metadata[6]-1) < ulocation.size()? (ulocation[i+2-1+metadata[6]-1]-ulocation[i+1-1+metadata[6]-1]): xp;

                    double yp = (j+1-1+metadata[7]-1) < vlocation.size()? (vlocation[j+1-1+metadata[7]-1]-vlocation[j-1+metadata[7]-1])
                                                                        : (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]);
                    double ys = (j-1-1+metadata[7]-1) >= 0? (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]): yp;
                    double yn = (j+2-1+metadata[7]-1) < vlocation.size()? (vlocation[j+2-1+metadata[7]-1]-vlocation[j+1-1+metadata[7]-1]): yp;

                    double zp = (k+1-1+metadata[8]-1) < wlocation.size()? (wlocation[k+1-1+metadata[8]-1]-wlocation[k-1+metadata[8]-1])
                                                                        : (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]);
                    double zf = (k-1-1+metadata[8]-1) >= 0? (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]): zp;
                    double zb = (k+2-1+metadata[8]-1) < wlocation.size()? (wlocation[k+2-1+metadata[8]-1]-wlocation[k+1-1+metadata[8]-1]): zp;

                    double facexy = xp * yp;
                    double facexz = xp * zp;
                    double faceyz = yp * zp;                  
                    double ap_u_approx = up[indexU(i, j, k)];
                    if (ap_u_approx == 0) ap_u_approx = 1e-10;
                    double correction = (p_prime[indexP(i-1,j,k)] - p_prime[indexP(i,j,k)]) * faceyz / ap_u_approx;
                    ux[getX(i-1,j-1,k-1)].u = u_star[indexU(i,j,k)] + correction;
                    u_star[indexU(i,j,k)] = ux[getX(i-1,j-1,k-1)].u;
                }
            }
        }

        // 修正V速度
        xStart = 1, xEnd = cx;
        yStart = 1, yEnd = fy;
        zStart = 1, zEnd = cz;        
        for (int k = zStart; k <= zEnd; ++k) {
            for (int j = yStart; j <= yEnd-1; ++j) {
                for (int i = xStart; i <= xEnd; ++i) {
                    if (p[getC(i-1, j-1, k-1)].isEmpty) {
                        // p[getC(i-1, j-1, k-1)].data = 0;
                        continue;
                    } 
                    double xp = (i+1-1+metadata[6]-1) < ulocation.size()? (ulocation[i+1-1+metadata[6]-1]-ulocation[i-1+metadata[6]-1])
                                                                        : (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]);
                    double xw = (i-1-1+metadata[6]-1) >= 0? (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]): xp;
                    double xe = (i+2-1+metadata[6]-1) < ulocation.size()? (ulocation[i+2-1+metadata[6]-1]-ulocation[i+1-1+metadata[6]-1]): xp;

                    double yp = (j+1-1+metadata[7]-1) < vlocation.size()? (vlocation[j+1-1+metadata[7]-1]-vlocation[j-1+metadata[7]-1])
                                                                        : (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]);
                    double ys = (j-1-1+metadata[7]-1) >= 0? (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]): yp;
                    double yn = (j+2-1+metadata[7]-1) < vlocation.size()? (vlocation[j+2-1+metadata[7]-1]-vlocation[j+1-1+metadata[7]-1]): yp;

                    double zp = (k+1-1+metadata[8]-1) < wlocation.size()? (wlocation[k+1-1+metadata[8]-1]-wlocation[k-1+metadata[8]-1])
                                                                        : (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]);
                    double zf = (k-1-1+metadata[8]-1) >= 0? (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]): zp;
                    double zb = (k+2-1+metadata[8]-1) < wlocation.size()? (wlocation[k+2-1+metadata[8]-1]-wlocation[k+1-1+metadata[8]-1]): zp;

                    double facexy = xp * yp;
                    double facexz = xp * zp;
                    double faceyz = yp * zp;                   
                    double ap_v_approx = vp[indexV(i, j, k)];
                    if (ap_v_approx == 0) ap_v_approx = 1e-10;
                    double correction = (p_prime[indexP(i,j-1,k)] - p_prime[indexP(i,j,k)]) * facexz / ap_v_approx;
                    uy[getY(i-1,j-1,k-1)].u = v_star[indexV(i,j,k)] + correction;
                    v_star[indexV(i,j,k)] = uy[getY(i-1,j-1,k-1)].u;
                }
            }
        }

        // 修正W速度
        xStart = 1, xEnd = cx;
        yStart = 1, yEnd = cy;
        zStart = 1, zEnd = fz;        
        for (int k = zStart; k <= zEnd-1; ++k) {
            for (int j = yStart; j <= yEnd; ++j) {
                for (int i = xStart; i <= xEnd; ++i) {
                    if (p[getC(i-1, j-1, k-1)].isEmpty) {
                        // p[getC(i-1, j-1, k-1)].data = 0;
                        continue;
                    } 
                    double xp = (i+1-1+metadata[6]-1) < ulocation.size()? (ulocation[i+1-1+metadata[6]-1]-ulocation[i-1+metadata[6]-1])
                                                                        : (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]);
                    double xw = (i-1-1+metadata[6]-1) >= 0? (ulocation[i-1+metadata[6]-1]-ulocation[i-1-1+metadata[6]-1]): xp;
                    double xe = (i+2-1+metadata[6]-1) < ulocation.size()? (ulocation[i+2-1+metadata[6]-1]-ulocation[i+1-1+metadata[6]-1]): xp;

                    double yp = (j+1-1+metadata[7]-1) < vlocation.size()? (vlocation[j+1-1+metadata[7]-1]-vlocation[j-1+metadata[7]-1])
                                                                        : (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]);
                    double ys = (j-1-1+metadata[7]-1) >= 0? (vlocation[j-1+metadata[7]-1]-vlocation[j-1-1+metadata[7]-1]): yp;
                    double yn = (j+2-1+metadata[7]-1) < vlocation.size()? (vlocation[j+2-1+metadata[7]-1]-vlocation[j+1-1+metadata[7]-1]): yp;

                    double zp = (k+1-1+metadata[8]-1) < wlocation.size()? (wlocation[k+1-1+metadata[8]-1]-wlocation[k-1+metadata[8]-1])
                                                                        : (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]);
                    double zf = (k-1-1+metadata[8]-1) >= 0? (wlocation[k-1+metadata[8]-1]-wlocation[k-1-1+metadata[8]-1]): zp;
                    double zb = (k+2-1+metadata[8]-1) < wlocation.size()? (wlocation[k+2-1+metadata[8]-1]-wlocation[k+1-1+metadata[8]-1]): zp;

                    double facexy = xp * yp;
                    double facexz = xp * zp;
                    double faceyz = yp * zp;                                 
                    double ap_w_approx = wp[indexW(i, j, k)];
                    if (ap_w_approx == 0) ap_w_approx = 1e-10;
                    double correction = (p_prime[indexP(i,j,k-1)] - p_prime[indexP(i,j,k)]) * facexy / ap_w_approx;
                    uz[getZ(i-1,j-1,k-1)].u = w_star[indexW(i,j,k)] + correction;
                    w_star[indexW(i,j,k)] =  uz[getZ(i-1,j-1,k-1)].u;
                }
            }
        }

              
            // 修正压力
        for (int k = 0; k < cz; ++k) {
            for (int j = 0; j < cy; ++j) {
                for (int i = 0; i < cx; ++i) {
                    if (p[getC(i, j, k)].isEmpty) {
                        // p[getC(i-1, j-1, k-1)].data = 0;
                        continue;
                    }                    
                    p[getC(i,j,k)].data += UNDER_RELAX * p_prime[indexP(i+1,j+1,k+1)];
                }
            }
        }

        for (int k = 0; k < cz; ++k) {
            for (int j = 0; j < cy; ++j) {
                for (int i = 0; i < cx ; ++i) {
                    if (p[getC(i, j, k)].isEmpty) {
                        // p[getC(i-1, j-1, k-1)].data = 0;
                        continue;
                    }                    
                    p_temp[indexP(i+1,j+1,k+1)] = p[getC(i,j,k)].data;
                }
            }
        }
        std::fill(p_prime.begin(), p_prime.end(), 0);
        for (auto &e: sendData) {
            e.second.resize(0);
        }
        for (auto &e: recvData) {
            e.second.resize(0);
        }

        for (auto &v: p_boundary) {
            for (auto &e: v.second) {
                sendData[v.first].push_back(p_temp[e]);
            }
            recvData[v.first].resize(recvData[v.first].size() + sendData[v.first].size()); 
        }
        for (auto &e: p_boundary) {
            MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 11,
                         recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 11,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // 为下次循环
        for (auto &v: p_halo) {
            for (int e = 0; e < v.second.size(); ++e) {
                p_temp[p_halo[v.first][e]] = recvData[v.first][e];
            }
        }                 
        // std::fill(p_prime.begin(), p_prime.end(), 0);
        // solve temperature equation
        xStart = 1, xEnd = cx;
        yStart = 1, yEnd = cy;
        zStart = 1, zEnd = cz;
        // std::cout << "solve fluid temp" << std::endl;    
        solveTemp(u_star, v_star, w_star, t_temp,
                   xStart, xEnd, yStart, yEnd, zStart, zEnd);          
                   
    }   
}
void FluidSolver::setForNext() {

    auto cx = metadata[0], cy = metadata[1], cz = metadata[2];  
    auto fx = metadata[3], fy = metadata[4], fz = metadata[5];      
    for (int k = 0; k < cz; ++k) {
        for (int j = 0; j < cy; ++j) {
            for (int i = 0; i < cx; ++i) {
                auto tmp  = t_temp[indexP(i+1,j+1,k+1)];
                t[getC(i,j,k)].data = tmp;
            }
        }
    }
    //std::cout << "compensate boundary conditions\n";
    applyXBoundaryCondition(u_star);
    applyYBoundaryCondition(v_star);
    applyZBoundaryCondition(w_star);
    applyCellBoundaryCondition(p_temp);
    for (int k = 0; k < cz; ++k) {
        for (int j = 0; j < cy; ++j) {
            for (int i = 0; i < cx; ++i) {
                p_temp[indexP(i+1, j+1, k+1)] = p[getC(i, j, k)].data;
            }
        }
    } 
    for (auto &e: sendData) {
        e.second.resize(0);
    }
    for (auto &e: recvData) {
        e.second.resize(0);
    }

    for (auto &v: p_boundary) {
        for (auto &e: v.second) {
            sendData[v.first].push_back(p_temp[e]);
        }
        recvData[v.first].resize(recvData[v.first].size() + sendData[v.first].size()); 
    }
    for (auto &e: p_boundary) {
        MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 12,
                        recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 12,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // 为下次循环
    for (auto &v: p_halo) {
        for (int e = 0; e < v.second.size(); ++e) {
            p_temp[p_halo[v.first][e]] = recvData[v.first][e];
        }
    } 
}

bool FluidSolver::parseString(const std::string& token, const std::string& s, int& num1, int& num2) {
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

void FluidSolver::HYPREIJMatrix_SolveCorrectionPressure(const std::vector<double> &u_star, const std::vector<double> &v_star,
                                     const std::vector<double> &w_star, std::vector<double> &p_prime,
                                     std::vector<double> & up, std::vector<double> & vp, std::vector<double> & wp,
    I64 cx, I64 cy, I64 cz, I64 offsetX, I64 offsetY, I64 offsetZ, I64 globalStart, I64 local_size) {

    // HYPRE_Initialize();
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // std::cout << "Rank "<< rank<< " " << std::endl;
    int num_procs = size;
    int myid = rank;

    int nx = mesh.cycasMesh.Nx;
    int ny = mesh.cycasMesh.Ny;
    int nz = mesh.cycasMesh.Nz;

    offsetX -= 1;
    offsetY -= 1;
    offsetZ -= 1;

    int xstart = offsetX;
    int xend = offsetX + cx - 1;

    int ystart = offsetY;
    int yend = offsetY + cy - 1;

    int zstart = offsetZ;
    int zend = offsetZ + cz - 1;
    // renumbering
    std::unordered_map<std::tuple<int, int, int>, int, KeyHash> local_map;
    int ilower = globalStart;
    int iupper = local_size - 1 + globalStart;
    int iter = 0;
    int nnz = 0;   
                             
    if (global_map.size() == 0) {
        initMap(size, xstart, xend, ystart, yend, zstart, zend, local_map, ilower, iter);
        initMatrix(local_size, ilower, iupper);
        initColsAndSetValues(cx, cy, cz, local_size, rank, nx, ny, nz, xstart, xend, ystart, yend, zstart,
                             zend, ilower, nnz, u_star, v_star, w_star, up, vp, wp);
    } else {
        setCoffValues(u_star, v_star, w_star, up, vp, wp, local_size, xstart, xend, ystart, yend, zstart, zend, ilower);         
    }

    HYPRE_IJMatrixAssemble(A);
    HYPRE_IJMatrixGetObject(A, (void **)&parcsr_A);
    /* Create the rhs and solution */

    HYPRE_IJVectorSetValues(b, local_size, rows.data(), b_values.data());
    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void **)&par_b);

    HYPRE_IJVectorSetValues(x, local_size, rows.data(), x_values.data());
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **)&par_x);

    // std::cout << "has already assembled\n";  
    int num_iterations;
    double final_res_norm;
    // std::cout << "start setup and solver" << std::endl;
/* Now set up the AMG preconditioner and specify any parameters */
/*
    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetPrintLevel(precond, 1);
    HYPRE_BoomerAMGSetCoarsenType(precond, 6);
    HYPRE_BoomerAMGSetOldDefault(precond);
    HYPRE_BoomerAMGSetRelaxType(precond, 6);
    HYPRE_BoomerAMGSetNumSweeps(precond, 1);
    HYPRE_BoomerAMGSetTol(precond, 0.0);
    HYPRE_BoomerAMGSetMaxIter(precond, 1);
    HYPRE_BoomerAMGSetStrongThreshold(solver, 0.6);
    HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                        (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);
*/
    // std::cout << "start solve" << std::endl;
    HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
    HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);
    /* Run info - needed logging turned on */
    HYPRE_PCGGetNumIterations(solver, &num_iterations);
    HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

    HYPRE_Int n = iupper - ilower + 1;
    std::vector<double> result(n);
    // 获取值
  
    HYPRE_IJVectorGetValues(x, n, rows.data(), result.data());
    for (int iter = ilower; iter <= iupper; ++iter) {
        if (result_map.find(iter) == result_map.end())
            continue;
        auto ix = std::get<0>(result_map[iter]); 
        auto iy = std::get<1>(result_map[iter]); 
        auto iz = std::get<2>(result_map[iter]); 
        auto i = ix - xstart;
        auto j = iy - ystart;
        auto k = iz - zstart;
        // std::cout << rank << " " << iter << " " << ix << " " << iy << " " << iz << std::endl;
        p_prime[indexP(i+1,j+1,k+1)] = result[iter - ilower];
    }

    for (auto &e: sendData) {
        e.second.resize(0);
    }
    for (auto &e: recvData) {
        e.second.resize(0);
    }

    for (auto &v: p_boundary) {
        for (auto &e: v.second) {
            sendData[v.first].push_back(p_prime[e]);
        }
        recvData[v.first].resize(recvData[v.first].size() + sendData[v.first].size()); 
    }
    for (auto &e: p_boundary) {
        MPI_Sendrecv(sendData[e.first].data(), sendData[e.first].size(), MPI_DOUBLE, e.first, 1,
                        recvData[e.first].data(), recvData[e.first].size(), MPI_DOUBLE, e.first, 1,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for (auto &v: p_halo) {
        for (int e = 0; e < v.second.size(); ++e) {
            p_prime[p_halo[v.first][e]] = recvData[v.first][e];
        }
    } 
                                      
}

void FluidSolver::initMatrix(I64 local_size, int ilower, int iupper) {
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);
    rows.resize(local_size);
    for (int i = 0; i < local_size; i++)
    {
        rows[i] = ilower + i; // global row number
    }
    flag.resize(local_size);
    for (auto &e: flag) {
        e.resize(7, 0);
    }
    values.resize(local_size);
    for (auto &e: values) {
        e = std::vector<double>();
    }
    cols.resize(local_size);
    for (auto &e: cols) {
        e = std::vector<int>();
    }
    b_values.resize(local_size, 0);
    x_values.resize(local_size, 0);
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);   
    
    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

    HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
    HYPRE_PCGSetTol(solver, 1e-5);     /* conv. tolerance */
    HYPRE_PCGSetTwoNorm(solver, 1);    /* use the two norm as the stopping criteria */
/*    HYPRE_PCGSetPrintLevel(solver, 2); 
    HYPRE_PCGSetLogging(solver, 1);  
        HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetPrintLevel(precond, 1);
    HYPRE_BoomerAMGSetCoarsenType(precond, 6);
    HYPRE_BoomerAMGSetOldDefault(precond);
    HYPRE_BoomerAMGSetRelaxType(precond, 6);
    HYPRE_BoomerAMGSetNumSweeps(precond, 1);
    HYPRE_BoomerAMGSetTol(precond, 0.0);
    HYPRE_BoomerAMGSetMaxIter(precond, 1);
    HYPRE_BoomerAMGSetStrongThreshold(solver, 0.6);
    HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                        (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);    */   
}

void FluidSolver::initMap(int size, int xstart, int xend, int ystart, int yend, int zstart, int zend,
                          std::unordered_map<std::tuple<int, int, int>, int, KeyHash> &local_map, int ilower,
                          int iter) {
    for (int k = zstart; k <= zend; k++) {
        for (int j = ystart; j <= yend; j++) {
            for (int i = xstart; i <= xend; i++) {
                int ix = i - xstart;
                int iy = j - ystart;
                int iz = k - zstart;
                if (p[getC(ix, iy, iz)].isEmpty)
                    continue;
                local_map[{i, j, k}] = ilower + iter;
                ++iter;
            }
        }
    }
    std::vector<int> local_data;
    for (const auto& [key, val] : local_map) {
        auto [i, j, k] = key;
        local_data.push_back(i);
        local_data.push_back(j);
        local_data.push_back(k);
        local_data.push_back(val);
    }

    int local_N = local_data.size();

    int *recv_counts = (int *)malloc(size * sizeof(int));
    MPI_Allgather(&local_N, 1, MPI_INT,  // 发送自己的 local_N
                recv_counts, 1, MPI_INT,  // 接收所有进程的 local_N
                MPI_COMM_WORLD);

    // 步骤2：计算位移 displs 和总长度 total_N
    int *displs = (int *)malloc(size * sizeof(int));
    int total_N = 0;
    for (int p = 0; p < size; p++) {
        displs[p] = total_N;  // 当前进程的位移
        total_N += recv_counts[p];  // 累加总长度
    }
    std::vector<int> global_data(total_N);
    MPI_Allgatherv(local_data.data(), local_N, MPI_INT,  // 发送本地数据
            global_data.data(), recv_counts, displs, MPI_INT, // 接收缓冲区（仅主进程有效）
            MPI_COMM_WORLD);

    for (int i = 0; i < total_N; i += 4) {
        int ii = global_data[i];
        int jj = global_data[i + 1];
        int kk = global_data[i + 2];
        int id = global_data[i + 3];
        global_map[{ii, jj, kk}] = id;
    }
    free(recv_counts);
    free(displs);
    for (auto &e: local_map) {
        result_map[e.second] = e.first;
    }    
}

void
FluidSolver::initColsAndSetValues(I64 cx, I64 cy, I64 cz, I64 local_size, int rank, int nx, int ny, int nz, int xstart,
                                  int xend, int ystart, int yend, int zstart, int zend, int ilower, int nnz,
                                  const std::vector<double> &u_star, const std::vector<double> &v_star,
                                  const std::vector<double> &w_star, std::vector<double> &up, std::vector<double> &vp,
                                  std::vector<double> &wp) {
    for (int ix = xstart; ix <= xend; ++ix) {
        for (int iy = ystart; iy <= yend; ++iy) {
            for (int iz = zstart; iz <= zend; ++iz) {
                if (global_map.find({ix, iy, iz}) == global_map.end())
                    continue;
                auto row = global_map[{ix, iy, iz}];
                auto i = ix - xstart;
                auto j = iy - ystart;
                auto k = iz - zstart;
                if (p[getC(i, j, k)].isEmpty)
                    continue;
                id[getC(i, j, k)] = row; // output
                nnz = 0;

                double ap_u_approx = up[indexU(i + 1, j + 1, k + 1)];
                double ap_v_approx = vp[indexV(i + 1, j + 1, k + 1)];
                double ap_w_approx = wp[indexW(i + 1, j + 1, k + 1)];
                double xp = (i + 1 + metadata[6] - 1) < ulocation.size() ? (ulocation[i + 1 + metadata[6] - 1] -
                                                                            ulocation[i + metadata[6] - 1])
                                                                            : (ulocation[i + metadata[6] - 1] -
                                                                            ulocation[i - 1 + metadata[6] - 1]);

                double yp = (j + 1 + metadata[7] - 1) < vlocation.size() ? (vlocation[j + 1 + metadata[7] - 1] -
                                                                            vlocation[j + metadata[7] - 1])
                                                                            : (vlocation[j + metadata[7] - 1] -
                                                                            vlocation[j - 1 + metadata[7] - 1]);

                double zp = (k + 1 + metadata[8] - 1) < wlocation.size() ? (wlocation[k + 1 + metadata[8] - 1] -
                                                                            wlocation[k + metadata[8] - 1])
                                                                            : (wlocation[k + metadata[8] - 1] -
                                                                            wlocation[k - 1 + metadata[8] - 1]);

                double facexy = xp * yp;
                double facexz = xp * zp;
                double faceyz = yp * zp;
                double a_e_prime = (faceyz) / ap_u_approx;
                double a_w_prime = (faceyz) / ap_u_approx;
                double a_n_prime = (facexz) / ap_v_approx;
                double a_s_prime = (facexz) / ap_v_approx;
                double a_f_prime = (facexy) / ap_w_approx;
                double a_b_prime = (facexy) / ap_w_approx;
                double a_p_prime = (a_e_prime + a_w_prime + a_n_prime + a_s_prime + a_f_prime + a_b_prime);
                
                // std::cout << getC(i,j,k) << " " ;
                b_values[row - ilower]
                        = -(u_star[indexU(i + 1 + 1, j + 1, k + 1)] - u_star[indexU(i + 1, j + 1, k + 1)])
                          -(v_star[indexV(i + 1, j + 1 + 1, k + 1)] - v_star[indexV(i + 1, j + 1, k + 1)])
                          -(w_star[indexW(i + 1, j + 1, k + 1 + 1)] - w_star[indexW(i + 1, j + 1, k + 1)]);
                int pid = -1;
                anotherRow.push_back(row);
                int idx1  = global_idx(ix, iy, iz - 1, nx, ny, nz,  global_map);
                if (i >= 0 && i <= cx - 1 && j >= 0 && j <= cy - 1 && k - 1 >= 0 && k-1 <= cz - 1
                    && p[getC(i, j, k - 1)].isEmpty
                    && uz[getZ(i, j, k)].isBoundary
                    && parseString("ib_", uz[getZ(i, j, k)].info->name, rank, pid)) idx1 = -1;
                if (idx1 >= 0)
                {
                    nnz++;
                    cols[row - ilower].push_back(idx1);
                    values[row - ilower].push_back(-a_f_prime);
                    flag[row - ilower][0] = 1;
                }
                int idx2  = global_idx(ix, iy - 1, iz, nx, ny, nz,  global_map);
                if (i >= 0 && i <= cx - 1 && j-1 >= 0 && j-1 <= cy - 1 && k >= 0 && k <= cz - 1
                    && p[getC(i, j - 1, k)].isEmpty
                    && uy[getY(i, j, k)].isBoundary
                    && parseString("ib_", uy[getY(i, j, k)].info->name, rank, pid)) idx2 = -1;
                if (idx2 >= 0)
                {
                    nnz++;
                    cols[row - ilower].push_back(idx2);
                    values[row - ilower].push_back(-a_s_prime);
                    flag[row - ilower][1] = 1;
                }
                int idx3  = global_idx(ix - 1, iy, iz, nx, ny, nz,  global_map);
                if (i-1 >= 0 && i-1 <= cx - 1 && j >= 0 && j <= cy - 1 && k >= 0 && k <= cz - 1
                    && p[getC(i - 1, j, k)].isEmpty
                    && ux[getX(i, j, k)].isBoundary
                    && parseString("ib_", ux[getX(i, j, k)].info->name, rank, pid)) idx3 = -1;
                if (idx3 >= 0)
                {
                    nnz++;
                    cols[row - ilower].push_back(idx3);
                    values[row - ilower].push_back(-a_w_prime);
                    flag[row - ilower][2] = 1;
                }

                nnz++;
                cols[row - ilower].push_back(row);
                values[row - ilower].push_back(a_p_prime);
                flag[row - ilower][3] = 1;

                int idx4  = global_idx(ix + 1, iy, iz, nx, ny, nz,  global_map);
                if (i+1 >= 0 && i+1 <= cx - 1 && j >= 0 && j <= cy - 1 && k >= 0 && k <= cz - 1
                    && p[getC(i + 1, j, k)].isEmpty
                    && ux[getX(i + 1, j, k)].isBoundary
                    && parseString("ib_", ux[getX(i + 1, j, k)].info->name, rank, pid)) idx4 = -1;
                if (idx4 >= 0)
                {
                    nnz++;
                    cols[row - ilower].push_back(idx4);
                    values[row - ilower].push_back(-a_e_prime);
                    flag[row - ilower][4] = 1;
                }

                int idx5 = global_idx(ix, iy + 1, iz, nx, ny, nz,  global_map);
                if (i >= 0 && i <= cx - 1 && j+1 >= 0 && j+1 <= cy - 1 && k >= 0 && k <= cz - 1
                    && p[getC(i, j + 1, k)].isEmpty
                    && uy[getY(i, j + 1, k)].isBoundary
                    && parseString("ib_", uy[getY(i, j + 1, k)].info->name, rank, pid)) idx5 = -1;
                if (idx5 >= 0 )
                {
                    nnz++;
                    cols[row - ilower].push_back(idx5);
                    values[row - ilower].push_back(-a_n_prime);
                    flag[row - ilower][5] = 1;
                }
                int idx6 = global_idx(ix, iy, iz + 1, nx, ny, nz,  global_map);
                if (i >= 0 && i <= cx - 1 && j >= 0 && j <= cy - 1 && k+1 >= 0 && k+1 <= cz - 1
                    && p[getC(i, j, k + 1)].isEmpty
                    && ux[getX(i, j, k + 1)].isBoundary
                    && parseString("ib_", ux[getX(i, j, k + 1)].info->name, rank, pid)) idx6 = -1;
                if (idx6 >= 0)
                {
                    nnz++;
                    cols[row - ilower].push_back(idx6);
                    values[row - ilower].push_back(-a_b_prime);
                    flag[row - ilower][6] = 1;
                }
                ncols.push_back(nnz);
            }
        }
    }
    for (int p = 0; p < anotherRow.size(); ++p) {
        auto row = anotherRow[p];
        for (auto &e: cols[row - ilower])
            cols1d.push_back(e);
        for (auto &e: values[row - ilower])
            values1d.push_back(e);
    }
    HYPRE_IJMatrixSetValues(A, local_size, ncols.data(), anotherRow.data(), cols1d.data(), values1d.data());
}

int FluidSolver::global_idx(int i, int j, int k, int nx, int ny, int nz, std::unordered_map<std::tuple<int, int, int>, int, KeyHash>& global_map) {
    if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) 
        return -1; // 边界外标记为-1
    if (global_map.find({i, j, k}) == global_map.end())
        return -1;
    else 
        return global_map[{i, j, k}];
}

void FluidSolver::setCoffValues(const std::vector<double> &u_star, const std::vector<double> &v_star,
                                const std::vector<double> &w_star, const std::vector<double> &up,
                                const std::vector<double> &vp, const std::vector<double> &wp, I64 local_size,
                                int xstart, int xend, int ystart, int yend, int zstart, int zend, int ilower) {
    int nnz;
    for (auto& e: result_map) {
        auto row = e.first;
        auto ix = std::get<0>(e.second);
        auto i = ix - xstart;
        auto iy = std::get<1>(e.second);
        auto j = iy - ystart;
        auto iz = std::get<2>(e.second);
        auto k = iz - zstart;
        this->id[this->getC(i, j, k)] = row;
        nnz = 0;

        double ap_u_approx = up[this->indexU(i + 1, j + 1, k + 1)];
        double ap_v_approx = vp[this->indexV(i + 1, j + 1, k + 1)];
        double ap_w_approx = wp[this->indexW(i + 1, j + 1, k + 1)];

        double xp = (i + 1 + this->metadata[6] - 1) < this->ulocation.size() ? (
                this->ulocation[i + 1 + this->metadata[6] - 1] -
                this->ulocation[i + this->metadata[6] - 1])
                                                                                : (
                            this->ulocation[i + this->metadata[6] - 1] -
                            this->ulocation[i - 1 + this->metadata[6] - 1]);

        double yp = (j + 1 + this->metadata[7] - 1) < this->vlocation.size() ? (
                this->vlocation[j + 1 + this->metadata[7] - 1] -
                this->vlocation[j + this->metadata[7] - 1])
                                                                                : (
                            this->vlocation[j + this->metadata[7] - 1] -
                            this->vlocation[j - 1 + this->metadata[7] - 1]);

        double zp = (k + 1 + this->metadata[8] - 1) < this->wlocation.size() ? (
                this->wlocation[k + 1 + this->metadata[8] - 1] -
                this->wlocation[k + this->metadata[8] - 1])
                                                                                : (
                            this->wlocation[k + this->metadata[8] - 1] -
                            this->wlocation[k - 1 + this->metadata[8] - 1]);

        double facexy = xp * yp;
        double facexz = xp * zp;
        double faceyz = yp * zp;
        double a_e_prime = (faceyz) / ap_u_approx;
        double a_w_prime = (faceyz) / ap_u_approx;
        double a_n_prime = (facexz) / ap_v_approx;
        double a_s_prime = (facexz) / ap_v_approx;
        double a_f_prime = (facexy) / ap_w_approx;
        double a_b_prime = (facexy) / ap_w_approx;
        double a_p_prime = (a_e_prime + a_w_prime + a_n_prime + a_s_prime + a_f_prime + a_b_prime);
        
        this->b_values[row - ilower]
                    = -(u_star[this->indexU(i + 1 + 1, j + 1, k + 1)] - u_star[this->indexU(i + 1, j + 1, k + 1)])
                        -(v_star[this->indexV(i + 1, j + 1 + 1, k + 1)] - v_star[this->indexV(i + 1, j + 1, k + 1)])
                        -(w_star[this->indexW(i + 1, j + 1, k + 1 + 1)] - w_star[this->indexW(i + 1, j + 1, k + 1)]);
        int it = 0;

        if (this->flag[row - ilower][0] == 1){
            this->values[row - ilower][it++] = (-a_f_prime);
        }
        if (this->flag[row - ilower][1] == 1) {
            this->values[row - ilower][it++] = (-a_s_prime);
        }
        if (this->flag[row - ilower][2] == 1) {
            this->values[row - ilower][it++] = (-a_w_prime);
        }
        if (this->flag[row - ilower][3] == 1) {
            this->values[row - ilower][it++] = (a_p_prime);
        }
        if (this->flag[row - ilower][4] == 1) {
            this->values[row - ilower][it++] = (-a_e_prime);
        }
        if (this->flag[row - ilower][5] == 1) {
            this->values[row - ilower][it++] = (-a_n_prime);
        }
        if (this->flag[row - ilower][6] == 1) {
            this->values[row - ilower][it++] = (-a_b_prime);
        }        
    }

    for (int p = 0; p < this->anotherRow.size(); ++p) {
        auto row = this->anotherRow[p];
        for (auto &e: this->values[row - ilower])
            this->values1d.push_back(e);
    }
    HYPRE_IJMatrixSetValues(this->A, local_size, this->ncols.data(), this->anotherRow.data(), this->cols1d.data(), this->values1d.data());
}

#endif
