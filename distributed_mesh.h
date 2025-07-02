#ifndef DISTRIBUTED_MESH
#define DISTRIBUTED_MESH

#include <iostream>
#include <algorithm>
#include "mpi.h"
#include "mesh/mesh.h"

struct BoundaryInfo {
    std::string name;
    std::string type;
    int nFaces;
    int myProcNo;
    int neighbProcNo;
};

struct LocalCell {
    double data = 0.0;
    bool isEmpty = true;
};

struct LocalFace {
    bool isBoundary = false;
    double u = 0.0;
    BoundaryInfo* info = nullptr;
    double location;
};


struct DistributedMesh {
    DistributedMesh(const CYCASMesh& mesh): cycasMesh(mesh) {
        idx2GlobalIJK();
    }
    CYCASMesh cycasMesh;
    std::map<int, std::vector<std::vector<I64>>> fluidGlobalIJK;
    std::map<int, std::vector<std::vector<I64>>> mobsGlobalIJK;

    /* 1. store values by vector, indexed by local index */
    std::map<int, std::vector<LocalCell>> p;
    std::map<int, std::vector<LocalFace>> ux, uy, uz; 
    std::map<int, std::vector<I64>> fluidMetadata;

    std::map<int, std::vector<LocalCell>> mt;
    //仅用于标记
    std::map<int, std::vector<LocalFace>> mx, my, mz;
    std::map<int, std::vector<I64>> solidMetadata;
    // std::map<int, std::vector<LocalFace>> sx, sy, sz; 

    /* 2. store and allocate according to global idx and global i j k */
    void indexToGlobalCoords(I64 index, I64& i, I64& j, I64& k) const;
    // void indexToGlobalCoordsFrom0(I64 index, I64& i, I64& j, I64& k) const;
    void generateLocalDomain();
    void initializeDomain();
    void idx2GlobalIJK();

    /* 3. access by local i j k*/
    I64 getCellIndex(I64 i, I64 j, I64 k, I64 nx, I64 ny, I64 nz) const;
    I64 getFaceIndex(I64 i, I64 j, I64 k, I64 nx, I64 ny, I64 nz) const;

    void getGlobalFaceIndexX(I64 globalId, I64 &i, I64 &j, I64 &k) const;
    void getGlobalFaceIndexY(I64 globalId, I64 &i, I64 &j, I64 &k) const;
    void getGlobalFaceIndexZ(I64 globalId, I64 &i, I64 &j, I64 &k) const;

};

void DistributedMesh::idx2GlobalIJK() {
    for (auto iter = cycasMesh.fluidRegions.begin(); iter != cycasMesh.fluidRegions.end(); ++iter) {
        auto regionId = iter->first;
        auto cellIdxArray = iter->second;
        for (decltype(cellIdxArray.size()) cell = 0; cell < cellIdxArray.size(); ++cell) {
            I64 idx = cellIdxArray[cell];
            I64 i, j, k;
            indexToGlobalCoords(idx, i, j, k);
            if (fluidGlobalIJK.find(regionId) == fluidGlobalIJK.end()) {
                fluidGlobalIJK[regionId] = std::vector<std::vector<I64>>();
            }
            fluidGlobalIJK[regionId].push_back({i, j, k});
        }
    } // 流体区域
    for (auto iter = cycasMesh.mobsRegions.begin(); iter != cycasMesh.mobsRegions.end(); ++iter) {
        auto regionId = iter->first;
        auto cellIdxArray = iter->second;
        for (decltype(cellIdxArray.size()) cell = 0; cell < cellIdxArray.size(); ++cell) {
            I64 idx = cellIdxArray[cell];
            I64 i, j, k;
            indexToGlobalCoords(idx, i, j, k);
            if (mobsGlobalIJK.find(regionId) == mobsGlobalIJK.end()) {
                mobsGlobalIJK[regionId] = std::vector<std::vector<I64>>();
            }
            mobsGlobalIJK[regionId].push_back({i, j, k});
        }
        
    } // 固体区域

}

void DistributedMesh::generateLocalDomain() {
    // resized by the number of fluid regions 
    for (auto &e: fluidGlobalIJK) {
        p[e.first] = std::vector<LocalCell>();
        ux[e.first]= std::vector<LocalFace>();
        uy[e.first]= std::vector<LocalFace>();
        uz[e.first]= std::vector<LocalFace>();       
    }
    for (auto &e: mobsGlobalIJK) {
        mt[e.first] = std::vector<LocalCell>(); 
        mx[e.first]= std::vector<LocalFace>();
        my[e.first]= std::vector<LocalFace>();
        mz[e.first]= std::vector<LocalFace>();   
    }

    auto Nx = cycasMesh.Nx, Ny = cycasMesh.Ny, Nz = cycasMesh.Nz;
    auto NxP = cycasMesh.NxP, NyP = cycasMesh.NyP;   
    auto NxF = NxP * Ny * Nz;
    auto NyF = Nx * NyP * Nz; // y方向面网格数量  
    // std::vector<std::vector<std::vector<I64>>> bd;
    for (auto iter = fluidGlobalIJK.begin(); iter != fluidGlobalIJK.end(); ++iter) {
        auto regionId = iter->first;
        auto globalCoords = iter->second;
        auto minXIter = min_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) { 
            return a[0] < b[0];
        });
        auto minYIter = min_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) {
            return a[1] < b[1];
        });
        auto minZIter = min_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) {
            return a[2] < b[2];
        }); 
        auto maxXIter = max_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) {
            return a[0] < b[0];
        });
        auto maxYIter = max_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) {
            return a[1] < b[1];
        });
        auto maxZIter = max_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) {
            return a[2] < b[2];
        });  

        I64 absoluteOffsetX = (*minXIter)[0], absoluteOffsetY = (*minYIter)[1], absoluteOffsetZ = (*minZIter)[2];
        I64 cellSizeX = (*maxXIter)[0] - (*minXIter)[0] + 1;
        I64 cellSizeY = (*maxYIter)[1] - (*minYIter)[1] + 1;
        I64 cellSizeZ = (*maxZIter)[2] - (*minZIter)[2] + 1;
        I64 faceSizeX = cellSizeX + 1;
        I64 faceSizeY = cellSizeY + 1;
        I64 faceSizeZ = cellSizeZ + 1;
        fluidMetadata[regionId] = std::vector<I64>();
        fluidMetadata[regionId].insert(fluidMetadata[regionId].end(), 
            {cellSizeX, cellSizeY, cellSizeZ, faceSizeX, faceSizeY, faceSizeZ, absoluteOffsetX, absoluteOffsetY, absoluteOffsetZ});

        // std::cout << "maxx: " << (*maxXIter)[0] << " maxy: " << (*maxYIter)[1] << " maxz: " << (*maxZIter)[2] << std::endl;
        // std::cout << "minx: " << (*minXIter)[0] << " miny: " << (*minYIter)[1] << " minz: " << (*minZIter)[2] << std::endl; 
        // std::cout << "cx: " << cellSizeX << " cy: " << cellSizeY << " cz: " << cellSizeZ << std::endl;
        // std::cout << "fx: " << faceSizeX << " fy: " << faceSizeY << " fz: " << faceSizeZ << std::endl; 

        p[regionId].resize(cellSizeX * cellSizeY * cellSizeZ);
        ux[regionId].resize(faceSizeX * cellSizeY * cellSizeZ);
        uy[regionId].resize(cellSizeX * faceSizeY * cellSizeZ);
        uz[regionId].resize(cellSizeX * cellSizeY * faceSizeZ);

        int count = 0;
        for(auto globalijk = globalCoords.begin(); globalijk != globalCoords.end(); ++globalijk) {
            auto &v = *globalijk;
            count++;
            p[regionId][getCellIndex(v[0] - absoluteOffsetX, v[1] - absoluteOffsetY, v[2] - absoluteOffsetZ,              
                cellSizeX, cellSizeY, cellSizeZ)].isEmpty = false;
        }
        std::cout << "count: " << count << std::endl;

        for(auto bdIter = cycasMesh.fluidBoundaries[regionId].begin(); bdIter != cycasMesh.fluidBoundaries[regionId].end(); ++bdIter) {
            auto name = bdIter->first;
            auto bdFaces = bdIter->second;
            auto startFace = bdFaces.startFace;
            BoundaryInfo* info = new BoundaryInfo();
            info->myProcNo = bdIter->second.myProcNo;
            info->neighbProcNo = bdIter->second.neighbProcNo;
            info->type = bdIter->second.type;
            info->name = bdIter->second.name;
            info->nFaces =  bdIter->second.nFaces;
            // bd.push_back(std::vector<std::vector<I64>>());
            for (auto iter = cycasMesh.fluidFaces[regionId].begin() + startFace; 
                    iter != cycasMesh.fluidFaces[regionId].begin() + startFace + bdFaces.nFaces; ++iter) {
                // I64 i = globalCoords[iter->owner][0], j = globalCoords[iter->owner][1], k = globalCoords[iter->owner][2]; 
                I64 i, j, k;
                if (iter->globalId > (NyF + NxF) - 1) {
                    // 0-4, 0-9, 0-5
                    getGlobalFaceIndexZ(iter->globalId, i, j, k);
                    // std::cout << i - absoluteOffsetX << " " << j - absoluteOffsetY << " " << k - absoluteOffsetZ << " z " << iter->globalId << std::endl; 
                    auto &f = uz[regionId][getFaceIndex(i - absoluteOffsetX, j - absoluteOffsetY, k - absoluteOffsetZ,
                                                        cellSizeX, cellSizeY, cellSizeZ)];
                    f.isBoundary = true;
                    f.info = info; 
                    // bd.back().push_back({bd.size() - 1, 2, i - absoluteOffsetX, j - absoluteOffsetY, k - absoluteOffsetZ});                  
                } else if (iter->globalId > NxF - 1) {
                    // NxF~NyF+NxF-1
                    // 0-4, 0-10, 0-4
                    getGlobalFaceIndexY(iter->globalId, i, j, k);
                    // std::cout << i - absoluteOffsetX << " " << j - absoluteOffsetY << " " << k - absoluteOffsetZ << " y " << iter->globalId << std::endl;
                    auto &f = uy[regionId][getFaceIndex(i - absoluteOffsetX, j - absoluteOffsetY, k - absoluteOffsetZ,
                                                        cellSizeX, faceSizeY, cellSizeZ)];
                    f.isBoundary = true; 
                    f.info = info; 
                    // bd.back().push_back({bd.size() - 1, 1, i - absoluteOffsetX, j - absoluteOffsetY, k - absoluteOffsetZ});                   
                } else {
                    // 0~NxF-1
                    // 0-5, 0-9, 0-4
                    getGlobalFaceIndexX(iter->globalId, i, j, k);
                    // std::cout << i - absoluteOffsetX << " " << j - absoluteOffsetY << " " << k - absoluteOffsetZ << " x " << iter->globalId << std::endl;
                    auto &f = ux[regionId][getFaceIndex(i - absoluteOffsetX, j - absoluteOffsetY, k - absoluteOffsetZ,
                                                        faceSizeX, cellSizeY, cellSizeZ)];
                    f.isBoundary = true;
                    f.info = info;
                    // bd.back().push_back({bd.size() - 1, 0, i - absoluteOffsetX, j - absoluteOffsetY, k - absoluteOffsetZ});
                }
            } 
            
        }
            int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);     
    std::cout << "rank: " << rank << std::endl;

    } // 流体区域   
    // for (auto &v: bd) {
    //     std::sort(v.begin(), v.end());
    //     for (auto &e: v) {
    //         std::cout << e[0] << " " << e[1] << " " << e[2] << " " << e[3] << " " << e[4] << std::endl;
    //     }
    // }

    for (auto iter = mobsGlobalIJK.begin(); iter != mobsGlobalIJK.end(); ++iter) {
        auto regionId = iter->first;
        auto globalCoords = iter->second;
        auto minXIter = min_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) { 
            return a[0] < b[0];
        });
        auto minYIter = min_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) {
            return a[1] < b[1];
        });
        auto minZIter = min_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) {
            return a[2] < b[2];
        }); 
        auto maxXIter = max_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) {
            return a[0] < b[0];
        });
        auto maxYIter = max_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) {
            return a[1] < b[1];
        });
        auto maxZIter = max_element(globalCoords.begin(), globalCoords.end(),
            [](const std::vector<I64>& a, const std::vector<I64>& b) {
            return a[2] < b[2];
        });  

        I64 absoluteOffsetX = (*minXIter)[0], absoluteOffsetY = (*minYIter)[1], absoluteOffsetZ = (*minZIter)[2];
        I64 cellSizeX = (*maxXIter)[0] - (*minXIter)[0] + 1;
        I64 cellSizeY = (*maxYIter)[1] - (*minYIter)[1] + 1;
        I64 cellSizeZ = (*maxZIter)[2] - (*minZIter)[2] + 1;
        I64 faceSizeX = cellSizeX + 1;
        I64 faceSizeY = cellSizeY + 1;
        I64 faceSizeZ = cellSizeZ + 1;

        solidMetadata[regionId] = std::vector<I64>();
        solidMetadata[regionId].insert(solidMetadata[regionId].end(), 
            {cellSizeX, cellSizeY, cellSizeZ, faceSizeX, faceSizeY, faceSizeZ, absoluteOffsetX, absoluteOffsetY, absoluteOffsetZ});
        // std::cout << "maxx: " << (*maxXIter)[0] << " maxy: " << (*maxYIter)[1] << " maxz: " << (*maxZIter)[2] << std::endl;
        // std::cout << "minx: " << (*minXIter)[0] << " miny: " << (*minYIter)[1] << " minz: " << (*minZIter)[2] << std::endl; 
        // std::cout << "cx: " << cellSizeX << " cy: " << cellSizeY << " cz: " << cellSizeZ << std::endl;
        // std::cout << "fx: " << faceSizeX << " fy: " << faceSizeY << " fz: " << faceSizeZ << std::endl; 
        mt[regionId].resize(cellSizeX * cellSizeY * cellSizeZ);
        mx[regionId].resize(faceSizeX * cellSizeY * cellSizeZ);
        my[regionId].resize(cellSizeX * faceSizeY * cellSizeZ);
        mz[regionId].resize(cellSizeX * cellSizeY * faceSizeZ);
        // access by solidT[regionId][getCellIndex(i, j, k)] 
        int count = 0;
        for(auto globalijk = globalCoords.begin(); globalijk != globalCoords.end(); ++globalijk) {
            auto &v = *globalijk;
            count++;
            mt[regionId][getCellIndex(v[0] - absoluteOffsetX, v[1] - absoluteOffsetY, v[2] - absoluteOffsetZ, 
                         cellSizeX, cellSizeY, cellSizeZ)].isEmpty = false;
        }
        std::cout << "count: " << count << std::endl;
        for(auto bdIter = cycasMesh.mobsBoundaries[regionId].begin(); bdIter !=cycasMesh.mobsBoundaries[regionId].end(); ++bdIter) {
            auto name = bdIter->first;
            auto bdFaces = bdIter->second;
            auto startFace = bdFaces.startFace;
            BoundaryInfo* info = new BoundaryInfo();
            info->myProcNo = bdIter->second.myProcNo;
            info->neighbProcNo = bdIter->second.neighbProcNo;
            info->type = bdIter->second.type;
            info->name = bdIter->second.name;
            info->nFaces =  bdIter->second.nFaces;
            for (auto iter = cycasMesh.mobsFaces[regionId].begin() + startFace; 
                    iter != cycasMesh.mobsFaces[regionId].begin() + startFace + bdFaces.nFaces; ++iter) {
                I64 i, j, k;
                // TODO 事实上只有边界有用，但其他内存也分配了，后续应该只分配边界空间
                if (iter->globalId > (NyF + NxF) - 1) {
                    // 0-4, 0-9, 0-5
                    getGlobalFaceIndexZ(iter->globalId, i, j, k);
                    auto &f = mz[regionId][getFaceIndex(i - absoluteOffsetX, j - absoluteOffsetY, k - absoluteOffsetZ, 
                                                        cellSizeX, cellSizeY, cellSizeZ)];
                    f.isBoundary = true;
                    f.info = info;                    
                } else if (iter->globalId > NxF - 1) {
                    getGlobalFaceIndexY(iter->globalId, i, j, k);
                    auto &f = my[regionId][getFaceIndex(i - absoluteOffsetX, j - absoluteOffsetY, k - absoluteOffsetZ, 
                                                        cellSizeX, faceSizeY, cellSizeZ)];
                    f.isBoundary = true; 
                    f.info = info;                    
                } else {
                    getGlobalFaceIndexX(iter->globalId, i, j, k);
                    auto &f = mx[regionId][getFaceIndex(i - absoluteOffsetX, j - absoluteOffsetY, k - absoluteOffsetZ, 
                                                        faceSizeX, cellSizeY, cellSizeZ)];
                    f.isBoundary = true;
                    f.info = info;
                }
            } 
            
        }

     } // 固体区域    
}

void DistributedMesh::indexToGlobalCoords(I64 index, I64& i, I64& j, I64& k) const {
    auto Nx = cycasMesh.Nx, Ny = cycasMesh.Ny;
    index -= ONE; // index from 1
    // global i j k
    k = index / (Nx * Ny);
    j = (index - k * Nx * Ny) / Nx;
    i = index - k * Nx * Ny - j * Nx;
    i += ONE;
    j += ONE;
    k += ONE;
}

inline I64 DistributedMesh::getCellIndex(I64 i, I64 j, I64 k, I64 nx, I64 ny, I64 nz) const {
    return (i) + (j)*nx + (k)*nx*ny;
}

inline I64 DistributedMesh::getFaceIndex(I64 i, I64 j, I64 k, I64 nx, I64 ny, I64 nz) const {
    return (i) + (j)*nx + (k)*nx*ny;
}

void DistributedMesh::getGlobalFaceIndexX(I64 globalId, I64 &i, I64 &j, I64 &k) const {
    auto Ny = cycasMesh.Ny;
    auto NxP = cycasMesh.NxP;     
    k = globalId / (NxP*Ny);
    j = (globalId - k * NxP * Ny) / NxP;
    i = globalId - k * NxP * Ny - j * NxP;
    i += ONE;
    j += ONE;
    k += ONE;
}

void DistributedMesh::getGlobalFaceIndexY(I64 globalId, I64 &i, I64 &j, I64 &k) const {
    auto Nx = cycasMesh.Nx, Ny = cycasMesh.Ny, Nz = cycasMesh.Nz;
    auto NxP = cycasMesh.NxP, NyP = cycasMesh.NyP;   
    auto NxF = NxP * Ny * Nz; // x方向面网格数量
    globalId -= NxF;
    k = globalId / (NyP*Nx);
    j = (globalId - k * NyP * Nx) / Nx;
    i = globalId - k * Nx * NyP - j * Nx;
    i += ONE;
    j += ONE;
    k += ONE;
}

void DistributedMesh::getGlobalFaceIndexZ(I64 globalId, I64 &i, I64 &j, I64 &k) const {
    auto Nx = cycasMesh.Nx, Ny = cycasMesh.Ny, Nz = cycasMesh.Nz;
    auto NxP = cycasMesh.NxP, NyP = cycasMesh.NyP;   
    auto NxF = NxP * Ny * Nz; // x方向面网格数量
    auto NyF = Nx * NyP * Nz; // y方向面网格数量   
    globalId -= (NxF + NyF);
    k = globalId / (Ny * Nx);
    j = (globalId - k * Ny * Nx) / Nx;
    i = globalId - k * Nx * Ny - j * Nx;
    i += ONE;
    j += ONE;
    k += ONE;
}

#endif