#include <iostream>
#include <cmath>
#include <queue>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <algorithm>
#include <cctype>
// #include <filesystem>

#include <sstream>
#include <string>
#include <cstdlib>

#include "mesh/mesh.h"
//#define METIS
#ifdef METIS
#include <metis.h>
#else
#include <scotch.h>
#endif

void sliceInputString(const std::string& input_str, std::string& output_str) {
    bool skip = false;

    // Remove parts between ';' and '\n'
    for (char c : input_str) {
        if (c == ';') {
            skip = true;
        } else if (c == '\n') {
            skip = false;
        } else if (!skip) {
            output_str += c;
        }
    }
}

StructuredMesh::StructuredMesh(std::string filename) {

    std::cout << "Info: initializing mesh ..." << std::endl;

    readMesh(filename);

    setDims(); // 设置网格维度

    std::cout << "Info: initializing cell type." << std::endl;
    initCellsType(); // 根据输入设置单元类型

    std::cout << "Info: initializing face type." << std::endl;
    initFacesType(); // 根据输入设置面类型

    std::cout << "Info: mesh initialized." << std::endl;

}

void StructuredMesh::setDims()
{
    NxP = xPoints.size();
    NyP = yPoints.size();
    NzP = zPoints.size();
    Nx = NxP - 1;
    if (perioFlag == 360)
        Ny = NyP;
    else
        Ny = NyP - 1;
    Nz = NzP - 1;
    NCells = Nx * Ny * Nz;
    NPoints = NxP * NyP * NzP;
    NxF = NxP * Ny * Nz; // x方向面网格数量
    NyF = Nx * NyP * Nz; // y方向面网格数量
    NzF = Nx * Ny * NzP; // z方向面网格数量
    NFaces = NxF + NyF + NzF;
            /* (Nx * Ny + Ny * Nz + Nz * Nx) // 外部面
             + ((Nx-1)*Ny*Nz + (Ny-1)*Nz*Nx + (Nz-1)*Nx*Ny) // 内部面 */
    NxNeigh = 1;
    NyNeigh = Nx;
    NzNeigh = Nx*Ny;
}

void StructuredMesh::initCellsType()
{
    domain.resize(NCells);
    //cellTypes.resize(NCells);

    for(I64 k=1; k<=Nz; ++k) {
        for(I64 j=1; j<=Ny; ++j) {
            for(I64 i=1; i<=Nx; ++i) {
                I64 index = getCellIndex(i, j, k);
                StCell& cell = getCell(index);
                cell.globalId = index; // 设置全局ID
                cell.type = CellType::Fluid; // 默认流体单元
                //cellTypes[index-1] = CellType::Fluid; // 默认流体单元
                cell.groupId = -1; // 默认组ID为-1
                cell.regionId = -1; // 默认区域ID为-1

                // 检查是否为障碍物单元
                // 这里假设障碍物单元是通过给定的坐标范围来定义的
                for(size_t m=0; m<mobsPoints.size(); ++m) {
                    if(i >= mobsPoints[m][0] && i <= mobsPoints[m][1]
                        && j >= mobsPoints[m][2] && j <= mobsPoints[m][3]
                        && k >= mobsPoints[m][4] && k <= mobsPoints[m][5]) {
                        cell.type = CellType::Obstacle; // 障碍物单元
                        //cellTypes[index-1] = CellType::Obstacle; // 障碍物单元
                        break;
                    }
                }
            }
        }
    }
}

void StructuredMesh::initFacesType()
{
    // 初始化面类型
    //domainFaces.resize(NFaces);
    faceTypes.resize(NFaces);

    // 内部面
    for(I64 k=2; k<Nz; ++k) {
        for(I64 j=2; j<Ny; ++j) {
            for(I64 i=2; i<Nx; ++i) {
                I64 faceIndex = getFaceIndexX(i, j, k);
                //domainFaces[faceIndex].type = FaceType::Internal;
                faceTypes[faceIndex] = FaceType::Internal;
            }
        }
    }

    // x方向边界面
    for(I64 k=1; k<=Nz; ++k) {
        for(I64 j=1; j<=Ny; ++j) {
            I64 faceIndex = getFaceIndexX(1, j, k);
            //domainFaces[faceIndex].type = FaceType::Wall; // 边界面
            faceTypes[faceIndex] = FaceType::Wall; // 边界面
            faceIndex = getFaceIndexX(NxP, j, k);
            //domainFaces[faceIndex].type = FaceType::Wall; // 边界面
            faceTypes[faceIndex] = FaceType::Wall; // 边界面
        }
    }

    // y方向边界面, perioFlag != 360视为圆周方向，不设置边界
    if (perioFlag != 360) {
        for(I64 k=1; k<=Nz; ++k) {
            for(I64 i=1; i<=Nx; ++i) {
                I64 faceIndex = getFaceIndexY(i, 1, k);
                //domainFaces[faceIndex].type = FaceType::Wall; // 边界面
                faceTypes[faceIndex] = FaceType::Wall; // 边界面
                faceIndex = getFaceIndexY(i, NyP, k);
                //domainFaces[faceIndex].type = FaceType::Wall; // 边界面
                faceTypes[faceIndex] = FaceType::Wall; // 边界面
            }
        }
    } else {
        for(I64 k=1; k<=Nz; ++k) {
            for(I64 i=1; i<=Nx; ++i) {
                I64 faceIndex = getFaceIndexY(i, 1, k);
                //domainFaces[faceIndex].type = FaceType::Wall; // 边界面
                faceTypes[faceIndex] = FaceType::Internal; // 边界面
                faceIndex = getFaceIndexY(i, NyP, k);
                //domainFaces[faceIndex].type = FaceType::Wall; // 边界面
                faceTypes[faceIndex] = FaceType::Internal; // 边界面
            }
        }        
    }

    // z方向边界面
    for(I64 j=1; j<=Ny; ++j) {
        for(I64 i=1; i<=Nx; ++i) {
            I64 faceIndex = getFaceIndexZ(i, j, 1);
            //domainFaces[faceIndex].type = FaceType::Wall; // 边界面
            faceTypes[faceIndex] = FaceType::Wall; // 边界面
            faceIndex = getFaceIndexZ(i, j, NzP);
            //domainFaces[faceIndex].type = FaceType::Wall; // 边界面
            faceTypes[faceIndex] = FaceType::Wall; // 边界面
        }
    }

    // 内部wall边界面
    for(size_t m=0; m<wallPoints.size(); ++m) {
        int i0 = wallPoints[m][0];
        int i1 = wallPoints[m][1];
        int j0 = wallPoints[m][2];
        int j1 = wallPoints[m][3];
        int k0 = wallPoints[m][4];
        int k1 = wallPoints[m][5];
        if(i0==i1) {
            for(I64 i=1; i<=NxP; ++i) {
                if(i==i0) {
                    for(I64 j=1; j<=Ny; ++j) {
                        for(I64 k=1; k<=Nz; ++k) {
                            if(j>=j0 && j<=j1 && k>=k0 && k<=k1) {
                                I64 faceIndex = getFaceIndexX(i, j, k);
                                //domainFaces[faceIndex].type = FaceType::Wall; // 墙面
                                faceTypes[faceIndex] = FaceType::Wall; // 边界面
                            }
                        }
                    }
                }
            }
        }

        if(j0==j1) {
            for(I64 j=1; j<=NyP; j++) {
                if(j==j0) {
                    for(I64 i=1; i<=Nx; ++i) {
                        for(I64 k=1; k<=Nz; ++k) {
                            if(i>=i0 && i<=i1 && k>=k0 && k<=k1) {
                                I64 faceIndex = getFaceIndexY(i, j, k);
                                //domainFaces[faceIndex].type = FaceType::Wall; // 墙面
                                faceTypes[faceIndex] = FaceType::Wall; // 边界面
                            }
                        }
                    }
                }
            }
        }

        if(k0==k1) {
            for(I64 k=1; k<=NzP; ++k) {
                if(k==k0) {
                    for(I64 i=1; i<=Nx; ++i) {
                        for(I64 j=1; j<=Ny; ++j) {
                            if(i>=i0 && i<=i1 && j>=j0 && j<=j1) {
                                I64 faceIndex = getFaceIndexZ(i, j, k);
                                //domainFaces[faceIndex].type = FaceType::Wall; // 墙面
                                faceTypes[faceIndex] = FaceType::Wall; // 边界面
                            }
                        }
                    }
                }
            }
        }
    }

    // 内部holes
    for(size_t m=0; m<holePoints.size(); ++m) {
        int i0 = holePoints[m][0];
        int i1 = holePoints[m][1];
        int j0 = holePoints[m][2];
        int j1 = holePoints[m][3];
        int k0 = holePoints[m][4];
        int k1 = holePoints[m][5];
        // if(i0==i1) {
        //     for(I64 i=1; i<=NxP; ++i) {
        //         if(i==i0) {
        //             for(I64 j=1; j<=Ny; ++j) {
        //                 for(I64 k=1; k<=Nz; ++k) {
        //                     if(j>=j0 && j<=j1 && k>=k0 && k<=k1 && (holePoints[m][7] == 1 || holePoints[m][8] == 1)) {
        //                         I64 faceIndex = getFaceIndexX(i, j, k);
        //                         //domainFaces[faceIndex].type = FaceType::Wall; // 墙面
        //                         faceTypes[faceIndex] = FaceType::Hole; // 边界面
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }

        // if(j0==j1) {
        //     for(I64 j=1; j<=NyP; j++) {
        //         if(j==j0) {
        //             for(I64 i=1; i<=Nx; ++i) {
        //                 for(I64 k=1; k<=Nz; ++k) {
        //                     if(i>=i0 && i<=i1 && k>=k0 && k<=k1 && (holePoints[m][9] == 1 || holePoints[m][10] == 1)) {
        //                         I64 faceIndex = getFaceIndexY(i, j, k);
        //                         //domainFaces[faceIndex].type = FaceType::Wall; // 墙面
        //                         faceTypes[faceIndex] = FaceType::Hole; // 边界面
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }

        // if(k0==k1) {
        //     for(I64 k=1; k<=NzP; ++k) {
        //         if(k==k0) {
        //             for(I64 i=1; i<=Nx; ++i) {
        //                 for(I64 j=1; j<=Ny; ++j) {
        //                     if(i>=i0 && i<=i1 && j>=j0 && j<=j1 && (holePoints[m][11] == 1 || holePoints[m][12] == 1)) {
        //                         I64 faceIndex = getFaceIndexZ(i, j, k);
        //                         //domainFaces[faceIndex].type = FaceType::Wall; // 墙面
        //                         faceTypes[faceIndex] = FaceType::Hole; // 边界面
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }
    }    
    
}

void StructuredMesh::getConnectedRegions(std::vector<std::vector<I64>>& fRegions,
                                    std::vector<std::vector<I64>>& mRegions)
{
    // 获取cell之间的连接关系
    std::vector<std::vector<I64>> cellAdjs(NCells);
    for(I64 k=1; k<=Nz; ++k) {
        for(I64 j=1; j<=Ny; ++j) {
            for(I64 i=1; i<=Nx; ++i) {
                I64 cellIndex = getCellIndex(i, j, k);
                StCell& cell = getCell(cellIndex);
                int numNeigh = 6;
                std::vector<I64> neighbourCellIndexes(numNeigh, -1);
                if(i>1) {
                    neighbourCellIndexes[0] = cellIndex - NxNeigh; // x-1
                }
                if(i<Nx) {
                    neighbourCellIndexes[1] = cellIndex + NxNeigh; // x+1
                }
                if(j>1) {
                    neighbourCellIndexes[2] = cellIndex - NyNeigh; // y-1
                }
                if(j<Ny) {
                    neighbourCellIndexes[3] = cellIndex + NyNeigh; // y+1
                }
                if (perioFlag == 360 && j == 1) {
                    neighbourCellIndexes[2] = cellIndex + NyNeigh * (Ny - 1);
                } 
                if (perioFlag == 360 && j == Ny) {
                    neighbourCellIndexes[3] = cellIndex - NyNeigh * (Ny - 1);
                }
                if(k>1) {
                    neighbourCellIndexes[4] = cellIndex - NzNeigh; // z-1
                }
                if(k<Nz) {
                    neighbourCellIndexes[5] = cellIndex + NzNeigh; // z+1
                }

                
                for(int n=0; n<numNeigh; ++n) {
                    if(neighbourCellIndexes[n] <=0 || neighbourCellIndexes[n] > NCells) {
                        // 邻接单元不存在，跳过
                        continue;
                    }

                    if(n==0) {
                        I64 faceIndex = getFaceIndexX(i, j, k);
                        StCell& neighCell = getCell(neighbourCellIndexes[n]);
                        if(cell.type==neighCell.type && faceTypes[faceIndex] != FaceType::Wall) {
                            cellAdjs[cellIndex-1].push_back(neighbourCellIndexes[n]-1);
                        }
                    } else if(n==1) {
                        I64 faceIndex = getFaceIndexX(i+1, j, k);
                        StCell& neighCell = getCell(neighbourCellIndexes[n]);
                        if(cell.type==neighCell.type && faceTypes[faceIndex] != FaceType::Wall) {
                            cellAdjs[cellIndex-1].push_back(neighbourCellIndexes[n]-1);
                        }
                    } else if(n==2) {
                        I64 faceIndex = getFaceIndexY(i, j, k);
                        StCell& neighCell = getCell(neighbourCellIndexes[n]);
                        if(cell.type==neighCell.type && faceTypes[faceIndex] != FaceType::Wall) {
                            cellAdjs[cellIndex-1].push_back(neighbourCellIndexes[n]-1);
                        }
                    } else if(n==3) {
                        I64 faceIndex = getFaceIndexY(i, j+1, k);
                        StCell& neighCell = getCell(neighbourCellIndexes[n]);
                        if(cell.type==neighCell.type && faceTypes[faceIndex] != FaceType::Wall) {
                            cellAdjs[cellIndex-1].push_back(neighbourCellIndexes[n]-1);
                        }
                    } else if(n==4) {
                        I64 faceIndex = getFaceIndexZ(i, j, k);
                        StCell& neighCell = getCell(neighbourCellIndexes[n]);
                        if(cell.type==neighCell.type && faceTypes[faceIndex] != FaceType::Wall) {
                            cellAdjs[cellIndex-1].push_back(neighbourCellIndexes[n]-1);
                        }
                    } else if(n==5) {
                        I64 faceIndex = getFaceIndexZ(i, j, k+1);
                        StCell& neighCell = getCell(neighbourCellIndexes[n]);
                        if(cell.type==neighCell.type && faceTypes[faceIndex] != FaceType::Wall) {
                            cellAdjs[cellIndex-1].push_back(neighbourCellIndexes[n]-1);
                        }
                    }
                }
            }
        }
    }

    // BFS遍历cellAdjs，获取连接的区域
    std::cout << "Info: BFS to get connected regions." << std::endl;
    std::vector<bool> visited(NCells, false);
    for(I64 i=0; i<NCells; i++)
    {
        if(!visited[i]) {
            std::vector<I64> region;
            BFS(i, cellAdjs, visited, region);
            if(region.size() > 0) {
                StCell& cell = getCell(region[0]+1);
                if(cell.type == CellType::Fluid) {
                    fRegions.push_back(region);  // cellId inside start from zero
                } else if(cell.type == CellType::Obstacle) {
                    mRegions.push_back(region); // cellId inside start from zero
                } else {
                    std::cerr << "TODO: Hole region." << std::endl;
                }
            }
        }
    }

    std::cout << "Info: number of fluid regions: " << fRegions.size() << std::endl;
    std::cout << "Info: number of obstacle regions: " << mRegions.size() << std::endl;

    // 排序
    for(size_t i=0; i<fRegions.size(); ++i) {
        for(size_t j=0; j<fRegions[i].size(); ++j) {
            I64 cellIndex = fRegions[i][j] + 1;
            StCell& cell = getCell(cellIndex);
            cell.regionId = i; // 设置区域ID
            if(cell.type != CellType::Fluid) {
                std::cerr << "Error: wrong type after BFS." << std::endl;
            }
        }
        fRegions[i].clear();
    }

    for(size_t i=0; i<mRegions.size(); ++i) {
        for(size_t j=0; j<mRegions[i].size(); ++j) {
            I64 cellIndex = mRegions[i][j] + 1;
            StCell& cell = getCell(cellIndex);
            cell.regionId = i; // 设置区域ID
            if(cell.type != CellType::Obstacle) {
                std::cerr << "Error: wrong type after BFS." << std::endl;
            }
        }
        mRegions[i].clear();
    }

    // 保证region里cell的顺序跟全局顺序一致
    for(size_t i=1; i<=NCells; ++i) {

       StCell& cell = getCell(i);
        if(cell.type == CellType::Fluid) {
            if(cell.regionId >= 0) {
                fRegions[cell.regionId].push_back(i);
            } else {
                std::cerr << "Fluid cell without region ID." << std::endl;
            }
        } else if(cell.type == CellType::Obstacle) {
            if(cell.regionId >= 0) {
                mRegions[cell.regionId].push_back(i);
            } else {
                std::cerr << "Obstacle cell without region ID." << std::endl;
            }
        } else {
            std::cerr << "TODO: Hole cell." << std::endl;
        }
    }

}

void StructuredMesh::flatten() {
    for (int i = 0 ;i < fluidRegions.size(); ++i) {

        std::vector<std::vector<I64>> procMetadata(fluidRegions[i].size());
        std::map<I64, I64> procSize;
        std::map<I64, I64> preSize; 
        std::vector<I64> sortedVector;
        std::map<std::pair<I64, I64>, std::vector<I64>> interaction; 
        sortedVector.resize(fluidRegions[i].size());
        for (int j = 0 ; j < fluidRegions[i].size(); ++j) {
            I64 minX = NCells, minY = NCells, minZ = NCells;
            I64 maxX = 0, maxY = 0, maxZ = 0;            
            procMetadata[j] = std::vector<I64>();
            for (int k = 0; k < fluidRegions[i][j].size(); ++k) {
                I64 cellIndex = fluidRegions[i][j][k]; // global
                I64 ix, iy, iz;
                indexToCoords(cellIndex, ix, iy, iz);
                minX = std::min(ix, minX);
                minY = std::min(iy, minY);
                minZ = std::min(iz, minZ);
                maxX = std::max(ix, maxX);
                maxY = std::max(iy, maxY);
                maxZ = std::max(iz, maxZ);
                preSize[j]++;
            }
            procMetadata[j].insert(procMetadata[j].end(), {minX, maxX, minY, maxY, minZ, maxZ});
            procSize[j] = ((maxZ - minZ)+1)*((maxY - minY)+1)*((maxX - minX)+1);
            sortedVector[j] = j;
        }
        // auto sorted by sorted map

        for (auto iter = procSize.begin(); iter != procSize.end(); ++iter) { //计数
            std::cout << "i: " << i << " " << iter->first << " " << iter->second;
            std::cout << " " << procMetadata[iter->first][0] << " " << procMetadata[iter->first][1] << " " <<
                      procMetadata[iter->first][2] << " " << procMetadata[iter->first][3] << " " << 
                      procMetadata[iter->first][4] << " " << procMetadata[iter->first][5] << std::endl;
            
        }
        for (auto iter = preSize.begin(); iter != preSize.end(); ++iter) {
            std::cout << "i: " << iter->first << " " << iter->second << std::endl;
        }
        std::sort(sortedVector.begin(), sortedVector.end(),
             [&](const I64& i, const I64& j) {return procSize[i] < procSize[j];});
        for (int k = 0; k < sortedVector.size(); ++k) {
            std::cout << sortedVector[k] << " " ;
        }
        std::cout << std::endl;

        for (int k1 = 0; k1 < sortedVector.size(); ++k1) {
            for (int k2 = k1 + 1; k2 < sortedVector.size(); ++k2) {
                auto j1 = sortedVector[k1];
                auto j2 = sortedVector[k2];
                auto p = std::make_pair(j1, j2);
                interaction.insert({p, std::vector<I64>()});
                auto left1 = std::max(procMetadata[j1][0], procMetadata[j2][0]);
                auto right1 = std::min(procMetadata[j1][1], procMetadata[j2][1]);
                if (right1 >= left1) {
                    interaction[p].push_back(left1);
                    interaction[p].push_back(right1);
                }
                auto left2 = std::max(procMetadata[j1][2], procMetadata[j2][2]);
                auto right2 = std::min(procMetadata[j1][3], procMetadata[j2][3]);
                if (right2 >= left2) {
                    interaction[p].push_back(left2);                    
                    interaction[p].push_back(right2);
                }
                auto left3 = std::max(procMetadata[j1][4], procMetadata[j2][4]);
                auto right3 = std::min(procMetadata[j1][5], procMetadata[j2][5]);                
                if (right3 >= left3) {
                    interaction[p].push_back(left3);                    
                    interaction[p].push_back(right3);
                }                                             
            }
        }
        for (auto &e: interaction) {
            if (e.second.size() < 6){
                e.second.clear();
            }
            std::cout << e.first.first << " " << e.first.second << " : value ";
            for (auto &v: e.second) std::cout << v << " " ;
            std::cout << std::endl;
        }
        // set baseline 
        for (int k1 = 0; k1 < sortedVector.size(); ++k1) {
            for (int k2 = k1 + 1; k2 < sortedVector.size(); ++k2) {
                auto j1 = sortedVector[k1];
                auto j2 = sortedVector[k2];    
                auto p = std::make_pair(j1, j2);
                if (interaction[p].size() == 0) 
                    continue;
                auto d1 = interaction[p][1] - interaction[p][0];
                auto d2 = interaction[p][3] - interaction[p][2];
                auto d3 = interaction[p][5] - interaction[p][4];
                auto ciend = interaction[p][1], cjend = interaction[p][3], ckend = interaction[p][5];
                // if (d1 < d2 && d1 < d3) {
                //     ciend = interaction[p][1]/2;
                // } else if (d2 < d1 && d2 < d3) {
                //     cjend = interaction[p][3]/2;
                // } else if (d3 < d1 && d3 < d2) {
                //     ckend = interaction[p][5]/2;
                // }
                for (I64 ci = interaction[p][0]; ci <= ciend; ci++) {
                    for (I64 cj = interaction[p][2]; cj <= cjend; cj++) {
                        for (I64 ck = interaction[p][4]; ck <= ckend; ck++) {
                            auto c = getCellIndex(ci, cj, ck);
                            StCell& cell = getCell(c);
                            if (cell.cellFlag == 0 && cell.groupId == j2) {
                                cell.groupId = j1;  
                                cell.procId = j1; 
                                cell.cellFlag = 1;      
                            }                  
                        }
                    }
                }
            }
        }    
    }
    // for(size_t i=0; i<fluidRegions.size(); ++i) {
    //     for (int j = 0; j < fluidRegions[i].size(); ++j)
    //         fluidRegions[i][j].clear();
    // }     
    // for(I64 i=1; i<=NCells; ++i) {
    //     StCell& cell = getCell(i);
    //     if(cell.type == CellType::Fluid) {
    //         if(cell.regionId >= 0) {
    //             if (fluidRegions[cell.regionId].size() <= cell.groupId) {
    //                 fluidRegions[cell.regionId].resize(cell.groupId + 1);
    //                 fluidRegions[cell.regionId][cell.groupId] = std::vector<I64>();
    //             }
    //             fluidRegions[cell.regionId][cell.groupId].push_back(i);
    //         } else {
    //             std::cerr << "Fluid cell without region ID." << std::endl;
    //         }
    //     } else if(cell.type == CellType::Obstacle) {
    //         if(cell.regionId >= 0) {
    //             mobsRegions[cell.regionId][0].push_back(i);
    //         } else {
    //             std::cerr << "Obstacle cell without region ID." << std::endl;
    //         }
    //     } else {
    //         std::cerr << "TODO: Hole cell." << std::endl;
    //     }
    // }
    // for (int i = 0 ;i < fluidRegions.size(); ++i) {

    //     std::vector<std::vector<I64>> procMetadata(fluidRegions[i].size());
    //     std::map<I64, I64> procSize;
    //     std::map<I64, I64> preSize; 
    //     std::vector<I64> sortedVector;
    //     std::map<std::pair<I64, I64>, std::vector<I64>> interaction; 
    //     sortedVector.resize(fluidRegions[i].size());
    //     for (int j = 0 ; j < fluidRegions[i].size(); ++j) {
    //         I64 minX = NCells, minY = NCells, minZ = NCells;
    //         I64 maxX = 0, maxY = 0, maxZ = 0;            
    //         procMetadata[j] = std::vector<I64>();
    //         for (int k = 0; k < fluidRegions[i][j].size(); ++k) {
    //             I64 cellIndex = fluidRegions[i][j][k]; // global
    //             I64 ix, iy, iz;
    //             indexToCoords(cellIndex, ix, iy, iz);
    //             minX = std::min(ix, minX);
    //             minY = std::min(iy, minY);
    //             minZ = std::min(iz, minZ);
    //             maxX = std::max(ix, maxX);
    //             maxY = std::max(iy, maxY);
    //             maxZ = std::max(iz, maxZ);
    //             preSize[j]++;
    //         }
    //         procMetadata[j].insert(procMetadata[j].end(), {minX, maxX, minY, maxY, minZ, maxZ});
    //         procSize[j] = ((maxZ - minZ)+1)*((maxY - minY)+1)*((maxX - minX)+1);
    //         sortedVector[j] = j;
    //     }
    //     // auto sorted by sorted map

    //     for (auto iter = procSize.begin(); iter != procSize.end(); ++iter) { //计数
    //         std::cout << "i: " << i << " " << iter->first << " " << iter->second;
    //         std::cout << " " << procMetadata[iter->first][0] << " " << procMetadata[iter->first][1] << " " <<
    //                   procMetadata[iter->first][2] << " " << procMetadata[iter->first][3] << " " << 
    //                   procMetadata[iter->first][4] << " " << procMetadata[iter->first][5] << std::endl;
            
    //     }
    //     for (auto iter = preSize.begin(); iter != preSize.end(); ++iter) {
    //         std::cout << "i: " << iter->first << " " << iter->second << std::endl;
    //     }
    //     std::sort(sortedVector.begin(), sortedVector.end(),
    //          [&](const I64& i, const I64& j) {return procSize[i] < procSize[j];});
    //     for (int k = 0; k < sortedVector.size(); ++k) {
    //         std::cout << sortedVector[k] << " " ;
    //     }
    //     std::cout << std::endl;

    //     for (int k1 = 0; k1 < sortedVector.size(); ++k1) {
    //         for (int k2 = k1 + 1; k2 < sortedVector.size(); ++k2) {
    //             auto j1 = sortedVector[k1];
    //             auto j2 = sortedVector[k2];
    //             auto p = std::make_pair(j1, j2);
    //             interaction.insert({p, std::vector<I64>()});
    //             auto left1 = std::max(procMetadata[j1][0], procMetadata[j2][0]);
    //             auto right1 = std::min(procMetadata[j1][1], procMetadata[j2][1]);
    //             if (right1 >= left1) {
    //                 interaction[p].push_back(left1);
    //                 interaction[p].push_back(right1);
    //             }
    //             auto left2 = std::max(procMetadata[j1][2], procMetadata[j2][2]);
    //             auto right2 = std::min(procMetadata[j1][3], procMetadata[j2][3]);
    //             if (right2 >= left2) {
    //                 interaction[p].push_back(left2);                    
    //                 interaction[p].push_back(right2);
    //             }
    //             auto left3 = std::max(procMetadata[j1][4], procMetadata[j2][4]);
    //             auto right3 = std::min(procMetadata[j1][5], procMetadata[j2][5]);                
    //             if (right3 >= left3) {
    //                 interaction[p].push_back(left3);                    
    //                 interaction[p].push_back(right3);
    //             }                                             
    //         }
    //     }
    //     for (auto &e: interaction) {
    //         if (e.second.size() < 6){
    //             e.second.clear();
    //         }
    //         std::cout << e.first.first << " " << e.first.second << " : value ";
    //         for (auto &v: e.second) std::cout << v << " " ;
    //         std::cout << std::endl;
    //     }
    //     // set baseline 
    //     for (int k1 = 0; k1 < sortedVector.size(); ++k1) {
    //          for (int k2 = k1 + 1; k2 < sortedVector.size(); ++k2) {
    //             auto j1 = sortedVector[k1];
    //             auto j2 = sortedVector[k2];    
    //             auto p = std::make_pair(j1, j2);
    //             if (interaction[p].size() == 0) 
    //                 continue;
    //             // auto d1 = interaction[p][1] - interaction[p][0];
    //             // auto d2 = interaction[p][3] - interaction[p][2];
    //             // auto d3 = interaction[p][5] - interaction[p][4];
    //             auto ciend = interaction[p][1], cjend = interaction[p][3], ckend = interaction[p][5];
    //             for (I64 ci = interaction[p][0]; ci <= ciend; ci++) {
    //                 for (I64 cj = interaction[p][2]; cj <= cjend; cj++) {
    //                     for (I64 ck = interaction[p][4]; ck <= ckend; ck++) {
    //                         auto c = getCellIndex(ci, cj, ck);
    //                         StCell& cell = getCell(c);
    //                         if (cell.cellFlag == 1 || cell.cellFlag == 0 && cell.groupId == j2) {
    //                             cell.groupId = j1;  
    //                             cell.procId = j1; 
    //                             cell.cellFlag = 2;      
    //                         }                  
    //                     }
    //                 }
    //             }
    //         }
    //     }    
    // }   
    // for(size_t i=0; i<fluidRegions.size(); ++i) {
    //     for (int j = 0; j < fluidRegions[i].size(); ++j)
    //         fluidRegions[i][j].clear();
    // }     
    // for(I64 i=1; i<=NCells; ++i) {
    //     StCell& cell = getCell(i);
    //     if(cell.type == CellType::Fluid) {
    //         if(cell.regionId >= 0) {
    //             if (fluidRegions[cell.regionId].size() <= cell.groupId) {
    //                 fluidRegions[cell.regionId].resize(cell.groupId + 1);
    //                 fluidRegions[cell.regionId][cell.groupId] = std::vector<I64>();
    //             }
    //             fluidRegions[cell.regionId][cell.groupId].push_back(i);
    //         } else {
    //             std::cerr << "Fluid cell without region ID." << std::endl;
    //         }
    //     } else if(cell.type == CellType::Obstacle) {
    //         if(cell.regionId >= 0) {
    //             mobsRegions[cell.regionId][0].push_back(i);
    //         } else {
    //             std::cerr << "Obstacle cell without region ID." << std::endl;
    //         }
    //     } else {
    //         std::cerr << "TODO: Hole cell." << std::endl;
    //     }
    // }      
}
void StructuredMesh::decomposeMesh(std::vector<std::vector<I64>>& fRegions,
    std::vector<std::vector<I64>>& mRegions,
    I32 nParts)
{

    // fluidRegions
    // 计算fluidRegions里面每个region占用几个processor
    I64 nFuildCells = 0;
    for(size_t i=0; i<fRegions.size(); ++i) {
        std::cout << "Info: region " << i << " has " << fRegions[i].size() << " cells." << std::endl;
        nFuildCells += fRegions[i].size();
    }

    std::vector<int> nProcsPerRegion(fRegions.size(), 0);
    int totalProcs = 0;
    // 按照每个region的cell数目给每个region分配处理器
    // TODO： 这里可以考虑使用更好的划分算法
    for(size_t i=0; i<fRegions.size(); ++i) {
        
        nProcsPerRegion[i] = std::round((double)fRegions[i].size() / (double)nFuildCells * nParts);
        if(nProcsPerRegion[i] < 1) {
            nProcsPerRegion[i] = 1;
        }
        if(nProcsPerRegion[i] > nParts) {
            nProcsPerRegion[i] = nParts;
        }
        totalProcs += nProcsPerRegion[i];
    }
    if(totalProcs > nParts) {
        std::cerr << "Error: total number of processors is greater than the number of processors." << std::endl;
        return;
    }
   
    if(nProcsPerRegion.size()>0 && totalProcs < nParts) {
        nProcsPerRegion[0] += nParts - totalProcs; // 将剩余的处理器分配给第一个区域
    }

    // 初始化容器
    fluidRegions.resize(fRegions.size());
    for(size_t i=0; i<fRegions.size(); ++i) {
        fluidRegions[i].resize(nProcsPerRegion[i]);
        for(size_t j=0; j<nProcsPerRegion[i]; ++j) {
            fluidRegions[i][j] = std::vector<I64>();
        }
    }

    for(size_t i=0; i<fRegions.size(); ++i) {

        // 只有一个处理器的region，直接分配处理器ID为0
        if(nProcsPerRegion[i] <= 1) {
            for(size_t j=0; j<fRegions[i].size(); ++j) {
                I64 cellIndex = fRegions[i][j];
                StCell& cell = getCell(cellIndex);
                cell.groupId = 0; // 设置处理器ID
                cell.regionId = i; // 设置区域ID
            }
            continue;
        }

        // 处理器个数大于1的region，使用METIS划分处理器
        std::vector<I64> visitedFaces(NFaces, -1);
        std::vector<I64> owners; // used to build adjacency list
        std::vector<I64> neighbours; // used to build adjacency list

        std::vector<I64>& region = fRegions[i];
        for(size_t j=0; j<region.size(); ++j) {
            I64 cellIndex = region[j];
            I64 ix, iy, iz;
            indexToCoords(cellIndex, ix, iy, iz);
            std::vector<I64> cellFaces;
            getFaceIndexes(ix, iy, iz, cellFaces);
            for(size_t k=0; k<cellFaces.size(); ++k) {
                if(visitedFaces[cellFaces[k]] == -1) {
                    visitedFaces[cellFaces[k]] = 0;
                } else if(visitedFaces[cellFaces[k]] == 0) {
                    visitedFaces[cellFaces[k]] = j;
                }
            }
        }

        for(size_t j=0; j<region.size(); ++j) {
            I64 cellIndex = region[j];
            I64 ix, iy, iz;
            indexToCoords(cellIndex, ix, iy, iz);
            std::vector<I64> cellFaces;
            getFaceIndexes(ix, iy, iz, cellFaces);
            for(size_t k=0; k<cellFaces.size(); ++k) {
                if(visitedFaces[cellFaces[k]] > 0
                    && faceTypes[cellFaces[k]] == FaceType::Internal) // Wall面没有连接
                {
                    owners.push_back(j);
                    neighbours.push_back(visitedFaces[cellFaces[k]]);
                    visitedFaces[cellFaces[k]] = -1; // 置为-1，表示已经使用过，壁面重复使用
                }
            }
        }

        std::vector<std::vector<I64>> adj(region.size());
        for(I64 k=0; k<owners.size(); ++k) {
            adj[owners[k]].push_back(neighbours[k]);
            adj[neighbours[k]].push_back(owners[k]);
        }
#ifdef METIS
        I64 nAdj = 0;
        for(I64 k=0; k<adj.size(); ++k) {
            nAdj += adj[k].size();
        }
        idx_t nVertices = region.size();
        idx_t connWeight = 1;
        idx_t np = nProcsPerRegion[i]; // 划分子域个数
        idx_t objVal;
        std::vector<idx_t> xAdj(nVertices+1);
        std::vector<idx_t> adjncy(nAdj);
        std::vector<idx_t> decomposedParts(nVertices, 0);
        // 组装邻接数组
        xAdj[0] = 0;
        for(idx_t k=0; k<nVertices; ++k) {
            xAdj[k+1] = xAdj[k] + adj[k].size();
            for(size_t j=0; j<adj[k].size(); j++) {
                adjncy[xAdj[k]+j] = adj[k][j];
            }
        }
        // 调用METIS进行划分
        if(np<=1) {
            std::cerr << "Error: number of processors is less than 1." << std::endl;
            return;
        }

        METIS_PartGraphKway(&nVertices, &connWeight, xAdj.data(), adjncy.data(),
                            NULL, NULL, NULL, &np,
                            NULL, NULL, NULL, &objVal, decomposedParts.data());
#else

        I64 nAdj = 0;
        for(I64 k=0; k<adj.size(); ++k) {
            nAdj += adj[k].size();
        }
        SCOTCH_Num nVertices = region.size();
        SCOTCH_Num connWeight = 1;
        SCOTCH_Num np = nProcsPerRegion[i]; // 划分子域个数
        SCOTCH_Num objVal;
        std::vector<SCOTCH_Num> xAdj(nVertices+1);
        std::vector<SCOTCH_Num> adjncy(nAdj);
        std::vector<SCOTCH_Num> decomposedParts(nVertices, 0);
        // 组装邻接数组
        xAdj[0] = 0;
        for(SCOTCH_Num k=0; k<nVertices; ++k) {
            xAdj[k+1] = xAdj[k] + adj[k].size();
            for(size_t j=0; j<adj[k].size(); j++) {
                adjncy[xAdj[k]+j] = adj[k][j];
            }
        }
        SCOTCH_Graph graph;
        SCOTCH_graphInit(&graph);
        SCOTCH_graphBuild(
            &graph,                // 图对象
            0,                     // 基准编号（0-based）
            nVertices,           // 顶点数
            xAdj.data(),        // 顶点指针数组
            nullptr,               // 顶点结束指针（可为nullptr）
            nullptr, nullptr,      // 顶点权重（无）
            adjncy.size(),        // 边数
            adjncy.data(),        // 边数组
            nullptr                // 边权重（无）
        );
        SCOTCH_Arch arch;
        SCOTCH_archInit(&arch);
        SCOTCH_archCmplt(&arch, np); // 创建完全连接的架构
    //     if(SCOTCH_archMesh3(&arch, 1, 2, 2)!=0) {
    //         std::cerr << "Error: SCOTCH_archMesh3 failed." << std::endl;
    //         //return;
	// } 
        // 创建划分策略
        SCOTCH_Strat strat;
        SCOTCH_stratInit(&strat);
        SCOTCH_stratGraphMapBuild(&strat, SCOTCH_STRATQUALITY, np, 0.05);

        // 执行划分
        int result = SCOTCH_graphMap(
            &graph,        // 图对象
            &arch,         // 目标架构
            &strat,        // 划分策略
            decomposedParts.data()    // 输出划分结果
        );
    
        if (result != 0) {
            std::cerr << "Scotch划分失败!错误代码: " << result << std::endl;
            return;
        }
#endif

        // 分配处理器ID
        for(size_t j=0; j<region.size(); ++j) {
            I64 cellIndex = region[j];
            StCell& cell = getCell(cellIndex);
            cell.groupId = decomposedParts[j]; // 设置处理器ID
            //cell.regionId = i; // 设置区域ID
        }
    }

    // mobsRegions
    // 计算mobsRegions里面每个region占用几个processor
    I64 nMobsCells = 0;
    for(size_t i=0; i<mRegions.size(); ++i) {
        std::cout << "Info: mobs region " << i << ": " << mRegions[i].size() << std::endl;
        nMobsCells += mRegions[i].size();
    }
    std::vector<int> nMobsProcsPerRegion(mRegions.size(), 0);
    int totalMobsProcs = 0;
    for(size_t i=0; i<mRegions.size(); ++i) {
        nMobsProcsPerRegion[i] = std::round((double)mRegions[i].size() / (double)nMobsCells * nParts);
	//std::cout << mRegions[i].size() << " " << nMobsCells << " " << nParts << " ";
        if(nMobsProcsPerRegion[i] < 1) {
            nMobsProcsPerRegion[i] = 1;
        }
        if(nMobsProcsPerRegion[i] > nParts) {
            nMobsProcsPerRegion[i] = nParts;
        }
        totalMobsProcs += nMobsProcsPerRegion[i];
	//std::cout << nMobsProcsPerRegion[i] << " " << totalMobsProcs << std::endl;
    }
    if(totalMobsProcs > nParts && totalMobsProcs - nParts == 1) {
	nMobsProcsPerRegion[0] -= 1;
	totalMobsProcs -= 1;
        //std::cerr << "Error: total number of processors is greater than the number of processors." << std::endl;
        //return;
    }
    if(nMobsProcsPerRegion.size()>0 && totalMobsProcs < nParts) {
        nMobsProcsPerRegion[0] += nParts - totalMobsProcs; // 将剩余的处理器分配给第一个区域
    }

    mobsRegions.resize(mRegions.size());
    for(size_t i=0; i<mRegions.size(); ++i) {
        mobsRegions[i].resize(nMobsProcsPerRegion[i]);
        for(size_t j=0; j<nMobsProcsPerRegion[i]; ++j) {
            mobsRegions[i][j] = std::vector<I64>();
        }
    }

    for(size_t i=0; i<mRegions.size(); ++i) {

        // 只有一个处理器的区域，直接分配处理器ID为0
        if(nMobsProcsPerRegion[i]<=1) {
            // 把所有groupId设为0
            for(size_t j=0; j<mRegions[i].size(); ++j) {
                I64 cellIndex = mRegions[i][j];
                StCell& cell = getCell(cellIndex);
                cell.groupId = 0; // 设置处理器ID
                cell.regionId = i; // 设置区域ID
            }
            continue;
        }

        std::vector<I64> visitedFaces(NFaces, -1);
        std::vector<I64> owners; // used to build adjacency list
        std::vector<I64> neighbours; // used to build adjacency list
        std::vector<I64>& region = mRegions[i];
        for(size_t j=0; j<region.size(); ++j) {
            I64 cellIndex = region[j];
            I64 ix, iy, iz;
            indexToCoords(cellIndex, ix, iy, iz);
            std::vector<I64> cellFaces;
            getFaceIndexes(ix, iy, iz, cellFaces);
            for(size_t k=0; k<cellFaces.size(); ++k) {
                if(visitedFaces[cellFaces[k]] == -1) {
                    visitedFaces[cellFaces[k]] = 0;
                } else if(visitedFaces[cellFaces[k]] == 0) {
                    visitedFaces[cellFaces[k]] = j;
                }
            }
        }
        for(size_t j=0; j<region.size(); ++j) {
            I64 cellIndex = region[j];
            I64 ix, iy, iz;
            indexToCoords(cellIndex, ix, iy, iz);
            std::vector<I64> cellFaces;
            getFaceIndexes(ix, iy, iz, cellFaces);
            for(size_t k=0; k<cellFaces.size(); ++k) {
                if(visitedFaces[cellFaces[k]] > 0
                    && faceTypes[cellFaces[k]] == FaceType::Internal) 
                {
                    
                    owners.push_back(j);
                    neighbours.push_back(visitedFaces[cellFaces[k]]);
                    visitedFaces[cellFaces[k]] = -1; // 置为-1，表示已经使用过，避免重复使用
                }
            }
        }

        std::vector<std::vector<I64>> adj(region.size());
        for(I64 k=0; k<owners.size(); ++k) {
            adj[owners[k]].push_back(neighbours[k]);
            adj[neighbours[k]].push_back(owners[k]);
        }

#ifdef METIS
        I64 nAdj = 0;
        for(I64 k=0; k<adj.size(); ++k) {
            nAdj += adj[k].size();
        }
        idx_t nVertices = region.size();
        idx_t connWeight = 1;
        idx_t np = nMobsProcsPerRegion[i]; // 划分子域个数
        idx_t objVal;
        std::vector<idx_t> xAdj(nVertices+1);
        std::vector<idx_t> adjncy(nAdj);
        std::vector<idx_t> decomposedParts(nVertices, 0);
        // 组装邻接数组
        xAdj[0] = 0;
        for(idx_t k=0; k<nVertices; ++k) {
            xAdj[k+1] = xAdj[k] + adj[k].size();
            for(size_t j=0; j<adj[k].size(); j++) {
                adjncy[xAdj[k]+j] = adj[k][j];
            }
        }

        // 调用METIS进行划分
        if(np<=1) {
            std::cerr << "Error: number of processors is less than 1." << std::endl;
            return;
        }
        METIS_PartGraphKway(&nVertices, &connWeight, xAdj.data(), adjncy.data(),
                            NULL, NULL, NULL, &np,
                            NULL, NULL, NULL, &objVal, decomposedParts.data());
#else
        I64 nAdj = 0;
        for(I64 k=0; k<adj.size(); ++k) {
            nAdj += adj[k].size();
        }
        SCOTCH_Num nVertices = region.size();
        SCOTCH_Num connWeight = 1;
        SCOTCH_Num np = nMobsProcsPerRegion[i]; // 划分子域个数
        SCOTCH_Num objVal;
        std::vector<SCOTCH_Num> xAdj(nVertices+1);
        std::vector<SCOTCH_Num> adjncy(nAdj);
        std::vector<SCOTCH_Num> decomposedParts(nVertices, 0);
        // 组装邻接数组
        xAdj[0] = 0;
        for(SCOTCH_Num k=0; k<nVertices; ++k) {
            xAdj[k+1] = xAdj[k] + adj[k].size();
            for(size_t j=0; j<adj[k].size(); j++) {
                adjncy[xAdj[k]+j] = adj[k][j];
            }
        }
        SCOTCH_Graph graph;
        SCOTCH_graphInit(&graph);
        SCOTCH_graphBuild(
            &graph,                // 图对象
            0,                     // 基准编号（0-based）
            nVertices,           // 顶点数
            xAdj.data(),        // 顶点指针数组
            nullptr,               // 顶点结束指针（可为nullptr）
            nullptr, nullptr,      // 顶点权重（无）
            adjncy.size(),        // 边数
            adjncy.data(),        // 边数组
            nullptr                // 边权重（无）
        );
        SCOTCH_Arch arch;
        SCOTCH_archInit(&arch);
        SCOTCH_archCmplt(&arch, np); // 创建完全连接的架构
        //SCOTCH_archMesh3D(&arch, 1, 2, 2);

        // 创建划分策略
        SCOTCH_Strat strat;
        SCOTCH_stratInit(&strat);
        SCOTCH_stratGraphMapBuild(&strat, SCOTCH_STRATQUALITY, np, 0.05);

        // 执行划分
        int result = SCOTCH_graphMap(
            &graph,        // 图对象
            &arch,         // 目标架构
            &strat,        // 划分策略
            decomposedParts.data()    // 输出划分结果
        );

        if (result != 0) {
            std::cerr << "Scotch划分失败!错误代码: " << result << std::endl;
            return;
        }
#endif
        // 分配处理器ID
        for(size_t j=0; j<region.size(); ++j) {
            I64 cellIndex = region[j];
            StCell& cell = getCell(cellIndex);
            cell.groupId = decomposedParts[j]; // 设置处理器ID
            //cell.regionId = i; // 设置区域ID
        }
    }

    // 遍历所有cell，分配到对应的region
    for(I64 i=1; i<=NCells; ++i) {
        StCell& cell = getCell(i);
        if(cell.type == CellType::Fluid) {
            if(cell.regionId >= 0) {
                fluidRegions[cell.regionId][cell.groupId].push_back(i);
            } else {
                std::cerr << "Fluid cell without region ID." << std::endl;
            }
        } else if(cell.type == CellType::Obstacle) {
            if(cell.regionId >= 0) {
                mobsRegions[cell.regionId][cell.groupId].push_back(i);
            } else {
                std::cerr << "Obstacle cell without region ID." << std::endl;
            }
        } else {
            std::cerr << "TODO: Hole cell." << std::endl;
        }
    }

    // 分配全局通信组中的ID
    for(int i=0; i<nParts; i++)
    {

        int proc =0;
        for(int j=0; j<fluidRegions.size(); j++)
        {
            if(i>=proc+fluidRegions[j].size())
            {
                proc += fluidRegions[j].size();
                continue;
            }

            for(int k=0; k<fluidRegions[j][i-proc].size(); ++k)
            {
                StCell& cell = getCell(fluidRegions[j][i-proc][k]);
                cell.procId = i; // 设置处理器ID
            }
        }

        proc = 0;
        for(int j=0; j<mobsRegions.size(); j++)
        {
            if(i>=proc+mobsRegions[j].size())
            {
                proc += mobsRegions[j].size();
                continue;
            }

            for(int k=0; k<mobsRegions[j][i-proc].size(); ++k)
            {
                StCell& cell = getCell(mobsRegions[j][i-proc][k]);
                cell.procId = i; // 设置处理器ID
            }
            break;
        }            
            
    }
    // flatten();
    // // 遍历所有cell，分配到对应的region
    // for(I64 i=1; i<=NCells; ++i) {
    //     StCell& cell = getCell(i);
    //     if(cell.type == CellType::Fluid) {
    //         if(cell.regionId >= 0) {
    //             fluidRegions[cell.regionId][cell.groupId].push_back(i);
    //         } else {
    //             std::cerr << "Fluid cell without region ID." << std::endl;
    //         }
    //     } else if(cell.type == CellType::Obstacle) {
    //         if(cell.regionId >= 0) {
    //             mobsRegions[cell.regionId][cell.groupId].push_back(i);
    //         } else {
    //             std::cerr << "Obstacle cell without region ID." << std::endl;
    //         }
    //     } else {
    //         std::cerr << "TODO: Hole cell." << std::endl;
    //     }
    // }

    // // 分配全局通信组中的ID
    // for(int i=0; i<nParts; i++)
    // {

    //     int proc =0;
    //     for(int j=0; j<fluidRegions.size(); j++)
    //     {
    //         if(i>=proc+fluidRegions[j].size())
    //         {
    //             proc += fluidRegions[j].size();
    //             continue;
    //         }

    //         for(int k=0; k<fluidRegions[j][i-proc].size(); ++k)
    //         {
    //             StCell& cell = getCell(fluidRegions[j][i-proc][k]);
    //             cell.procId = i; // 设置处理器ID
    //         }
    //     }

    //     proc = 0;
    //     for(int j=0; j<mobsRegions.size(); j++)
    //     {
    //         if(i>=proc+mobsRegions[j].size())
    //         {
    //             proc += mobsRegions[j].size();
    //             continue;
    //         }

    //         for(int k=0; k<mobsRegions[j][i-proc].size(); ++k)
    //         {
    //             StCell& cell = getCell(mobsRegions[j][i-proc][k]);
    //             cell.procId = i; // 设置处理器ID
    //         }
    //         break;
    //     }            
            
    // }
    // 初始化fluidFaces
    fluidFaces.resize(fluidRegions.size());
    fluidBoundaries.resize(fluidRegions.size());
    for(int i=0; i<fluidRegions.size(); ++i) {
        fluidFaces[i].resize(fluidRegions[i].size());
        fluidBoundaries[i].resize(fluidRegions[i].size());
    }

    // 遍历流体域中的每一个region
    for(int i=0; i<fluidRegions.size(); ++i) {
        // 遍历单个region中的每一个processor
        for(int j=0; j<fluidRegions[i].size(); ++j) {

            // 判断每一个processor的每一个face的类型
            std::vector<I64> visitedFaces(NFaces, -1);  // TODO：用一个map代替
            std::vector<StFace> faces;
            std::map<std::string, std::vector<StFace>> boundaryFaces; // 临时存储边界面，后续分类排序后存入faces

            // 遍历每一个cell，touch一下每个面，被touch 2次的是processor内部面
            for(int k=0; k<fluidRegions[i][j].size(); ++k) {
                I64 cellIndex = fluidRegions[i][j][k];
                StCell& cell = getCell(cellIndex);
                I64 ix, iy, iz;
                indexToCoords(cellIndex, ix, iy, iz);
                std::vector<I64> cellFaces;
                getFaceIndexes(ix, iy, iz, cellFaces);
                for(size_t l=0; l<cellFaces.size(); ++l) {
                    if(visitedFaces[cellFaces[l]] == -1) {
                        visitedFaces[cellFaces[l]] = 0;
                    } else if(visitedFaces[cellFaces[l]] == 0) {
                        visitedFaces[cellFaces[l]] = k; // 第二次访问，说明k是neighbourcell
                    }
                }
            }

            //std::vector<I64> owners;
            //std::vector<I64> neighbours;
            for(int k=0; k<fluidRegions[i][j].size(); ++k) {
                I64 cellIndex = fluidRegions[i][j][k];
                StCell& cell = getCell(cellIndex);
                I64 ix, iy, iz;
                indexToCoords(cellIndex, ix, iy, iz);
                std::vector<I64> cellFaces;
                getFaceIndexes(ix, iy, iz, cellFaces);

                std::vector<I64> neighbourCells(6, -1);
                neighbourCells[0] = cellIndex - NxNeigh;
                neighbourCells[1] = cellIndex + NxNeigh;
                neighbourCells[2] = cellIndex - NyNeigh;
                neighbourCells[3] = cellIndex + NyNeigh;
                neighbourCells[4] = cellIndex - NzNeigh;
                neighbourCells[5] = cellIndex + NzNeigh;

                for(size_t l=0; l<cellFaces.size(); ++l) {
                    if(visitedFaces[cellFaces[l]] > 0
                        && faceTypes[cellFaces[l]] == FaceType::Internal) // 内部面
                    {
                        //owners.push_back(k);
                        //neighbours.push_back(visitedFaces[cellFaces[l]]);
                        
                        StFace f;
                        f.owner = k;
                        f.globalId = cellFaces[l];
                        f.neighbour = visitedFaces[cellFaces[l]]; // 这里的neighbour是cellIndex在region中的索引
                        faces.push_back(f);
                        visitedFaces[cellFaces[l]] = -1; // 置为-1，表示已经访问过

                    } else { // 在processor意义上的边界面，对这些面进行分类
                        if(faceTypes[cellFaces[l]] == FaceType::Wall) {
                            std::string bdName = "wall";
                            if(boundaryFaces.find(bdName) == boundaryFaces.end()) {
                                boundaryFaces[bdName] = std::vector<StFace>();
                            }
                            StFace f;
                            f.owner = k;
                            f.globalId = cellFaces[l];
                            f.groupName = bdName;
                            boundaryFaces[bdName].push_back(f); // 墙面
                        } else if(faceTypes[cellFaces[l]] == FaceType::Boundary) {
                            std::string bdName = "boundary";
                            if(boundaryFaces.find(bdName) == boundaryFaces.end()) {
                                boundaryFaces[bdName] = std::vector<StFace>();
                            }
                            StFace f;
                            f.owner = k;
                            f.globalId = cellFaces[l];
                            f.groupName = bdName;
                            boundaryFaces[bdName].push_back(f); // 边界面
                        } else if(faceTypes[cellFaces[l]] == FaceType::Internal 
                                    && neighbourCells[l] > 0 
                                    && neighbourCells[l] <= NCells) {
                            // 判断隔壁的cell是否在同一个region和type，如果是，则设为Processor，否则设为Interface
                            StCell& cell2 = getCell(neighbourCells[l]);
                            if(cell2.type == cell.type) { // 同类型cell的processor边界
                                if(cell.regionId == cell2.regionId) {
                                    if(cell.groupId != cell2.groupId) {
                                        std::string bdName = "pr_"+std::to_string(cell.procId) + "_" + std::to_string(cell2.procId);
                                        if(boundaryFaces.find(bdName) == boundaryFaces.end()) {
                                            boundaryFaces[bdName] = std::vector<StFace>();
                                        }
                                        StFace f;
                                        f.owner = k;
                                        f.globalId = cellFaces[l];
                                        f.groupName = bdName;
                                        boundaryFaces[bdName].push_back(f); // 处理器面
                                    }
                                } else {
                                    // TODO: unknown situation
                                    std::cerr << "Error: cell region ID not match." << std::endl;
                                }
                            } else { // 不同类型cell的内部边界
                                std::string bdName = "ib_"+std::to_string(cell.procId) + "_" + std::to_string(cell2.procId);
                                if(boundaryFaces.find(bdName) == boundaryFaces.end()) {
                                    boundaryFaces[bdName] = std::vector<StFace>();
                                }
                                StFace f;
                                f.owner = k;
                                f.globalId = cellFaces[l];
                                f.groupName = bdName;
                                boundaryFaces[bdName].push_back(f); // 内部边界面
                            }
                        } else {
                            // TODO: unknown situation
                            std::cerr << "Error: unknown face type." << std::endl;
                        }

                    }
                }
            }
        
            // 对每一类boundary face根据globalId进行排序，保证通信的时候跟邻接进程一一对应
            // it是同一类face的集合，it2集合中的每一个面
            
            for(auto it=boundaryFaces.begin(); it != boundaryFaces.end(); ++it)
            {
                std::string bdName = it->first;
                std::vector<StFace>& bfaces = it->second;
                std::map<I64, StFace> faceMap;
                BoundaryFaces bdFaces;
                bdFaces.name = bdName;
                bdFaces.startFace = faces.size(); // 记录开始位置
                for(size_t k=0; k<bfaces.size(); ++k) {
                    faceMap[bfaces[k].globalId] = bfaces[k];
                }

                // 遍历faceMap，获取owners、neighbours和boundary
                for(auto it2=faceMap.begin(); it2 != faceMap.end(); ++it2) {
                    I64 faceIndex = it2->first;
                    StFace& f = it2->second;
                    faces.push_back(f); // 墙面或边界面
                }

                bdFaces.nFaces = faces.size() - bdFaces.startFace; // 记录结束位置
                if(bdName.substr(0, 3) == "pr_") {
                    bdFaces.type = "processor"; // 处理器面
                    bdFaces.myProcNo = std::stoi(bdName.substr(3, bdName.find("_", 3)-3)); // 获取处理器ID
                    bdFaces.neighbProcNo = std::stoi(bdName.substr(bdName.find("_", 3)+1)); // 获取邻接处理器ID
                } else if(bdName.substr(0, 3) == "ib_") {
                    bdFaces.type = "internal"; // 内部边界面
                    bdFaces.myProcNo = std::stoi(bdName.substr(3, bdName.find("_", 3)-3)); // 获取处理器ID
                    bdFaces.neighbProcNo = std::stoi(bdName.substr(bdName.find("_", 3)+1)); // 获取邻接处理器ID
                } else {
                    bdFaces.type = bdName; // 墙面或边界面
                }

                fluidBoundaries[i][j][bdName] = bdFaces; // 把bdFaces加入到类成员中
                
            }


            fluidFaces[i][j] = faces; // 把faces加入到类成员中
            
            // 把owners、neighbours和boundary加入到类成员中
        }
    }

    // 初始化mobsFaces
    mobsFaces.resize(mobsRegions.size());
    mobsBoundaries.resize(mobsRegions.size());
    for(int i=0; i<mobsRegions.size(); ++i) {
        mobsFaces[i].resize(mobsRegions[i].size());
        mobsBoundaries[i].resize(mobsRegions[i].size());
    }

    // 遍历障碍物中的每一个region
    for(int i=0; i<mobsRegions.size(); ++i) {
        // 遍历单个region中的每一个processor
        for(int j=0; j<mobsRegions[i].size(); ++j) {
            std::vector<I64> visitedFaces(NFaces, -1);  // TODO：用一个map代替
            std::vector<StFace> faces;
            std::map<std::string, std::vector<StFace>> boundaryFaces; // 临时存储边界面，后续分类排序后存入faces

            // 遍历每一个cell，touch一下每个面，被touch 2次的是processor内部面
            for(int k=0; k<mobsRegions[i][j].size(); ++k) {
                I64 cellIndex = mobsRegions[i][j][k];
                StCell& cell = getCell(cellIndex);
                I64 ix, iy, iz;
                indexToCoords(cellIndex, ix, iy, iz);
                std::vector<I64> cellFaces;
                getFaceIndexes(ix, iy, iz, cellFaces);
                for(size_t l=0; l<cellFaces.size(); ++l) {
                    if(visitedFaces[cellFaces[l]] == -1) {
                        visitedFaces[cellFaces[l]] = 0;
                    } else if(visitedFaces[cellFaces[l]] == 0) {
                        visitedFaces[cellFaces[l]] = k; // 第二次访问，说明k是neighbourcell
                    }
                }
            }

            for(int k=0; k<mobsRegions[i][j].size(); ++k) {
                I64 cellIndex = mobsRegions[i][j][k];
                StCell& cell = getCell(cellIndex);
                I64 ix, iy, iz;
                indexToCoords(cellIndex, ix, iy, iz);
                std::vector<I64> cellFaces;
                getFaceIndexes(ix, iy, iz, cellFaces);

                std::vector<I64> neighbourCells(6, -1);
                neighbourCells[0] = cellIndex - NxNeigh;
                neighbourCells[1] = cellIndex + NxNeigh;
                neighbourCells[2] = cellIndex - NyNeigh;
                neighbourCells[3] = cellIndex + NyNeigh;
                neighbourCells[4] = cellIndex - NzNeigh;
                neighbourCells[5] = cellIndex + NzNeigh;

                for(size_t l=0; l<cellFaces.size(); ++l) {
                    if(visitedFaces[cellFaces[l]] >
                        0 && faceTypes[cellFaces[l]] == FaceType::Internal) // 内部面
                    {
                        StFace f;
                        f.owner = k;
                        f.globalId = cellFaces[l];
                        f.neighbour = visitedFaces[cellFaces[l]]; // 这里的neighbour是cellIndex在region中的索引
                        faces.push_back(f);
                        visitedFaces[cellFaces[l]] = -1; // 置为-1，表示已经访问过

                    } else { // 在processor意义上的边界面，对这些面进行分类
                        if(faceTypes[cellFaces[l]] == FaceType::Wall) {
                            std::string bdName = "wall";
                            if(boundaryFaces.find(bdName) == boundaryFaces.end()) {
                                boundaryFaces[bdName] = std::vector<StFace>();
                            }
                            StFace f;
                            f.owner = k;
                            f.globalId = cellFaces[l];
                            f.groupName = bdName;
                            boundaryFaces[bdName].push_back(f); // 墙面
                        } else if(faceTypes[cellFaces[l]] == FaceType::Boundary) {
                            std::string bdName = "boundary";
                            if(boundaryFaces.find(bdName) == boundaryFaces.end()) {
                                boundaryFaces[bdName] = std::vector<StFace>();
                            }
                            StFace f;
                            f.owner = k;
                            f.globalId = cellFaces[l];
                            f.groupName = bdName;
                            boundaryFaces[bdName].push_back(f); // 边界面
                        } else if(faceTypes[cellFaces[l]] == FaceType::Internal 
                                    && neighbourCells[l] > 0 
                                    && neighbourCells[l] <= NCells) {
                            // 判断隔壁的cell是否在同一个region和type，如果是，则设为Processor，否则设为Interface
                            StCell& cell2 = getCell(neighbourCells[l]);
                            if(cell2.type == cell.type) { // 同类型cell的processor边界
                                if(cell.regionId == cell2.regionId) {
                                    if(cell.groupId != cell2.groupId) {
                                        std::string bdName = "pr_"+std::to_string(cell.procId) + "_" + std::to_string(cell2.procId);
                                        if(boundaryFaces.find(bdName) == boundaryFaces.end()) {
                                            boundaryFaces[bdName] = std::vector<StFace>();
                                        }
                                        StFace f;
                                        f.owner = k;
                                        f.globalId = cellFaces[l];
                                        f.groupName = bdName;
                                        boundaryFaces[bdName].push_back(f); // 处理器面
                                    }
                                } else {
                                    // TODO: unknown situation
                                    std::cerr << "Error: cell region ID not match." << std::endl;
                                }
                            } else { // 不同类型cell的内部边界
                                std::string bdName = "ib_"+std::to_string(cell.procId) + "_" + std::to_string(cell2.procId);
                                if(boundaryFaces.find(bdName) == boundaryFaces.end()) {
                                    boundaryFaces[bdName] = std::vector<StFace>();
                                }
                                StFace f;
                                f.owner = k;
                                f.globalId = cellFaces[l];
                                f.groupName = bdName;
                                boundaryFaces[bdName].push_back(f); // 内部边界面
                            }
                        } else {
                            // TODO: unknown situation
                            std::cerr << "Error: unknown face type." << std::endl;
                        }
                    }
                }
            }

            // 对每一类boundary face根据globalId进行排序，保证通信的时候跟邻接进程一一对应
            // it是同一类face的集合，it2集合中的每一个面
            for(auto it=boundaryFaces.begin(); it != boundaryFaces.end(); ++it)
            {
                std::string bdName = it->first;
                std::vector<StFace>& bfaces = it->second;
                std::map<I64, StFace> faceMap;
                BoundaryFaces bdFaces;
                bdFaces.name = bdName;
                bdFaces.startFace = faces.size(); // 记录开始位置
                for(size_t k=0; k<bfaces.size(); ++k) {
                    faceMap[bfaces[k].globalId] = bfaces[k];
                }

                // 遍历faceMap，获取owners、neighbours和boundary
                for(auto it2=faceMap.begin(); it2 != faceMap.end(); ++it2) {
                    I64 faceIndex = it2->first;
                    StFace& f = it2->second;
                    faces.push_back(f); // 墙面或边界面
                }

                bdFaces.nFaces = faces.size() - bdFaces.startFace; // 记录结束位置
                if(bdName.substr(0, 3) == "pr_") {
                    bdFaces.type = "processor"; // 处理器面
                    bdFaces.myProcNo = std::stoi(bdName.substr(3, bdName.find("_", 3)-3)); // 获取处理器ID
                    bdFaces.neighbProcNo = std::stoi(bdName.substr(bdName.find("_", 3)+1)); // 获取邻接处理器ID
                } else if(bdName.substr(0, 3) == "ib_") {
                    bdFaces.type = "internal"; // 内部边界面
                    bdFaces.myProcNo = std::stoi(bdName.substr(3, bdName.find("_", 3)-3)); // 获取处理器ID
                    bdFaces.neighbProcNo = std::stoi(bdName.substr(bdName.find("_", 3)+1)); // 获取邻接处理器ID
                } else {
                    bdFaces.type = bdName; // 墙面或边界面
                }
                // 把bdFaces加入到类成员中
                mobsBoundaries[i][j][bdName] = bdFaces; // 把bdFaces加入到类成员中
                
            }

            mobsFaces[i][j] = faces; // 把faces加入到类成员中
        }
    }
      
}

void StructuredMesh::BFS(I64 start, std::vector<std::vector<I64>>& cellAdjs,
    std::vector<bool>& visited, std::vector<I64>& region)
{
    std::queue<I64> q;
    q.push(start);
    visited[start] = true;
    region.push_back(start);

    while(!q.empty()) {
        I64 cellIndex = q.front();
        q.pop();

        for(U64 i=0; i<cellAdjs[cellIndex].size(); ++i) {
            I64 adjCellIndex = cellAdjs[cellIndex][i];
            if(!visited[adjCellIndex]) {
                visited[adjCellIndex] = true;
                q.push(adjCellIndex);
                region.push_back(adjCellIndex);
            }
        }
    }
}

void StructuredMesh::createMesh()
{

}

StructuredMesh::~StructuredMesh() 
{

}

void StructuredMesh::readMesh(std::string filename) {
    // 读取网格文件的实现
    // 这里可以使用类似于读取网格文件的方式来读取配置文件
    // 例如，使用std::ifstream打开文件并解析内容

    // 读入 xPoints, yPoints, zPoints
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return;
    }

    // dummy example
    /*
    xpoints = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]
    ypoints = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]
    zpoints = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]
    mobs = [0.4 0.6 0.2 0.5 0.4 0.6]
    wall = [0.5 0.5 0.2 0.3 0.4 0.5]
    */
    std::ostringstream ss;
    ss << file.rdbuf();  // 读取整个文件内容到 ss
    std::string content = ss.str();

    auto xgridPos = content.find("xgrid");
    auto ygridPos = content.find("ygrid"), ytemp = ygridPos;
    auto zgridPos = content.find("zgrid");
    auto mobsPos = content.find("mobs");
    auto wallsPos = content.find("walls");
    auto holesPos = content.find("holes");
    if (ygridPos == std::string::npos) {
        ygridPos = content.find("nyr(1)=");
    }
    std::vector<int> sortedVector;
    sortedVector.push_back(xgridPos);
    sortedVector.push_back(ygridPos);
    sortedVector.push_back(zgridPos);
    sortedVector.push_back(mobsPos);
    sortedVector.push_back(wallsPos);
    sortedVector.push_back(holesPos);
    std::sort(sortedVector.begin(), sortedVector.end());

    std::string sliceString, line, name;
    char eq, spliter;
    double value1;
    int value2;
    auto xstart = content.substr(xgridPos).find_first_of('=');
    int endPosition = content.size() - 1;
    for (int i = 0; i < sortedVector.size(); ++i) {
        if (sortedVector[i] == xgridPos && i != sortedVector.size() - 1) {
            endPosition = i + 1;
            break;
        }
    }
    auto length = (endPosition == content.size() - 1) ? content.size() - xstart - 2: sortedVector[endPosition] - xgridPos;
    auto s1 =  content.substr(xgridPos).substr(xstart + 1, length);
    std::istringstream iss1(s1);
    while(std::getline(iss1, line)) {
        sliceString.clear();
        sliceInputString(line, sliceString);
        std::istringstream lineiss(sliceString);
        while (lineiss >> value1 >> spliter)
            xPoints.push_back(value1);
    }
    if (ytemp == std::string::npos) {
        auto cylindricalPos = content.find("nyr(1)=");
        auto commaPos = content.find(",", cylindricalPos);
        std::string num_str = content.substr(cylindricalPos + 7, commaPos - (cylindricalPos + 7));
        auto num = std::stoi(num_str);

        cylindricalPos = content.find("yl(1)=");
        commaPos = content.find(",", cylindricalPos);
        num_str = content.substr(cylindricalPos + 6, commaPos - (cylindricalPos + 6));
        auto start_angle = std::stoi(num_str);

        cylindricalPos = content.find("yl(2)=");
        commaPos = content.find(",", cylindricalPos);
        num_str = content.substr(cylindricalPos + 6, commaPos - (cylindricalPos + 6));
        auto end_angle = std::stoi(num_str);
        for (int i = 0; i < num; ++i) {
            yPoints.push_back((end_angle - start_angle) / num * i);
        }
        std::cout << "size: " << yPoints.size() << std::endl;
        perioFlag = end_angle - start_angle;
    } else {
        auto ystart = content.substr(ygridPos).find_first_of('=');
        endPosition = content.size() - 1;
        for (int i = 0; i < sortedVector.size(); ++i) {
            if (sortedVector[i] == ygridPos && i != sortedVector.size() - 1) {
                endPosition = i + 1;
                break;
            }
        }
        auto length = (endPosition == content.size() - 1) ? content.size() - ystart - 2: sortedVector[endPosition] - ygridPos;                 
        std::istringstream iss2(content.substr(ygridPos).substr(ystart + 1, length));
        while(std::getline(iss2, line)) {
            sliceString.clear();
            sliceInputString(line, sliceString);
            std::istringstream lineiss(sliceString);
            while (lineiss >> value1 >> spliter)
                yPoints.push_back(value1);
        }          
    }
    endPosition = content.size() - 1;
    for (int i = 0; i < sortedVector.size(); ++i) {
        if (sortedVector[i] == zgridPos && i != sortedVector.size() - 1) {
            endPosition = i + 1;
            break;
        }
    } 
    auto zstart = content.substr(zgridPos).find_first_of('=');
    length = (endPosition == content.size() - 1) ? content.size() - zstart - 2: sortedVector[endPosition] - zgridPos;
    std::istringstream iss3(content.substr(zgridPos).substr(zstart + 1, length));
    while(std::getline(iss3, line)) {
        sliceString.clear();
        sliceInputString(line, sliceString);
        std::istringstream lineiss(sliceString);
        while (lineiss >> value1 >> spliter)
            zPoints.push_back(value1);
    }    
    endPosition = content.size() - 1;
    for (int i = 0; i < sortedVector.size(); ++i) {
        if (sortedVector[i] == mobsPos && i != sortedVector.size() - 1) {
            endPosition = i + 1;
            break;
        }
    }
    auto mobsstart = content.substr(mobsPos).find_first_of('=');
    length = (endPosition == (content.size() - 1)) ? content.size() - mobsstart - 2: sortedVector[endPosition] - mobsPos; 
    std::string s = content.substr(mobsPos).substr(mobsstart + 1, length);
    std::istringstream iss4(s);
    while(std::getline(iss4, line)) {
        std::vector<int> mobs;
        sliceString.clear();
        sliceInputString(line, sliceString);
        std::istringstream lineiss(sliceString);
        while (lineiss >> value2 >> spliter)
            mobs.push_back(value2);
        if(mobs.size() != 0)
            mobsPoints.push_back(mobs);
    } 
    endPosition = content.size() - 1;
    for (int i = 0; i < sortedVector.size(); ++i) {
        if (sortedVector[i] == wallsPos && i != sortedVector.size() - 1) {
            endPosition = i + 1;
            break;
        }
    }
    auto wallstart = content.substr(wallsPos).find_first_of('=');
    length = (endPosition == content.size() - 1) ? content.size() - wallstart - 2: sortedVector[endPosition] - wallsPos;
    std::istringstream iss5(content.substr(wallsPos).substr(wallstart + 1, length));
    while(std::getline(iss5, line)) {
        std::vector<int> walls;
        sliceString.clear();
        sliceInputString(line, sliceString);
        std::istringstream lineiss(sliceString);
        while (lineiss >> value2 >> spliter)
            walls.push_back(value2);
        if (walls.size() != 0 && walls.size() == 8)
            wallPoints.push_back(walls);
    }   
    endPosition = content.size() - 1;
    for (int i = 0; i < sortedVector.size(); ++i) {
        if (sortedVector[i] == holesPos && i != sortedVector.size() - 1) {
            endPosition = i + 1;
            break;
        }
    }    
    auto holestart = content.substr(holesPos).find_first_of('=');
    length = (endPosition == content.size() - 1) ? content.size() - holestart - 2: sortedVector[endPosition] - holesPos;
    std::istringstream iss6(content.substr(holesPos).substr(holestart + 1, length));
    while(std::getline(iss6, line)) {
        std::vector<int> holes;
        sliceString.clear();
        sliceInputString(line, sliceString);
        std::istringstream lineiss(sliceString);
        while (lineiss >> value2 >> spliter)
            holes.push_back(value2);
        if (holes.size() != 0 && holes.size() == 13)
            holePoints.push_back(holes);
    }    
    
    // std::string line, subline;
    // while (std::getline(file, line)) {
    //     line = ltrim(line);
    //     std::istringstream iss(line);

    //     if (line[0] == ';') 
    //         continue;
    //     if (iss >> name >> eq) {
    //         if (name == "xpoints") {
    //             while (iss >> value1 >> spliter) {
    //                 xPoints.push_back(value1);
    //             }
    //             while (std::getline(file, subline)) {
    //                 subline = ltrim(line);
    //                 std::istringstream subiss(subline);
    //                 if (subline[0] == ';') 
    //                     continue;
    //                 if (std::isalpha(subline[0]) || std::ispunct(subline[0]))
    //                     break;
    //                 while (subiss >> value1 >> spliter) {
    //                     xPoints.push_back(value1);
    //                 }
    //             }
    //         } else if (name == "ypoints") {
    //             while (iss >> value1 >> spliter) {
    //                 yPoints.push_back(value1);
    //             }
    //             while (std::getline(file, subline)) {
    //                 subline = ltrim(line);
    //                 std::istringstream subiss(subline);
    //                 if (subline[0] == ';') 
    //                     continue;
    //                 if (std::isalpha(subline[0]) || std::ispunct(subline[0]))
    //                     break;
    //                 while (subiss >> value1 >> spliter) {
    //                     yPoints.push_back(value1);
    //                 }
    //             }                
    //         } else if (name == "zpoints") {
    //             while (iss >> value1 >> spliter) {
    //                 zPoints.push_back(value1);
    //             }
    //             while (std::getline(file, subline)) {
    //                 subline = ltrim(line);
    //                 std::istringstream subiss(subline);
    //                 if (subline[0] == ';') 
    //                     continue;
    //                 if (std::isalpha(subline[0]) || std::ispunct(subline[0])) //字母开头或者"开头
    //                     break;
    //                 while (subiss >> value1 >> spliter) {
    //                     zPoints.push_back(value1);
    //                 }
    //             }                
    //         } else if (name == "mobs") {
    //             std::vector<int> mobs;
    //             while (iss >> value2 >> spliter) {
    //                 mobs.push_back(value2);
    //             }
    //             mobsPoints.push_back(mobs);
    //             while (std::getline(file, subline)) {
    //                 mobs.clear();
    //                 subline = ltrim(line);
    //                 std::istringstream subiss(subline);
    //                 if (subline[0] == ';') 
    //                     continue;
    //                 if (std::isalpha(subline[0]) || std::ispunct(subline[0])) //字母开头或者"开头
    //                     break;
    //                 while (subiss >> value2 >> spliter) {
    //                     mobs.push_back(value2);
    //                 }
    //                 mobsPoints.push_back(mobs);
    //             }                 
    //         } else if (name == "walls") {
    //             std::vector<int> wall;
    //             while (iss >> value2 >> spliter) {
    //                 wall.push_back(value2);
    //             }
    //             wallPoints.push_back(wall);
    //             while (std::getline(file, subline)) {
    //                 wall.clear();
    //                 subline = ltrim(line);
    //                 std::istringstream subiss(subline);
    //                 if (subline[0] == ';') 
    //                     continue;
    //                 if (std::isalpha(subline[0]) || std::ispunct(subline[0])) //字母开头或者"开头
    //                     break;
    //                 while (subiss >> value2 >> spliter) {
    //                     wall.push_back(value2);
    //                 }
    //                 wallPoints.push_back(wall);
    //             }                
    //         }
    //     }
    //     if (subline != "" && std::isalpha(subline[0]))
    // }
}

void StructuredMesh::writeMesh(std::string prefix) {
    // 写入网格文件的实现
    // 这里可以使用类似于读取网格文件的方式来读取配置文件
    // 例如，使用std::ofstream打开文件并解析内容

    std::vector<std::vector<Point>> fluidRegionPoints; // 流体区域点坐标
    std::vector<std::vector<Point>> mobsRegionPoints; // obstacle区域点坐标
    std::vector<std::vector<Cell>> fluidRegionCells; // 流体区域单元
    std::vector<std::vector<Cell>> mobsRegionCells; // obstacle区域单元

    std::string filename = prefix + ".vtk";
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return;
    }
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "vtk output" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    file << "POINTS " << NxP*NyP*NzP << " float" << std::endl;
    for(int k=0; k<NzP; ++k) {
        for(int j=0; j<NyP; ++j) {
            for(int i=0; i<NxP; ++i) {
                file << xPoints[i] << " "
                     << yPoints[j] << " "
                     << zPoints[k] << std::endl;
            }
        }
    }
    file << "CELLS " << NCells << " " << NCells*9 << std::endl;
    for(int k=0; k<Nz; ++k) {
        for(int j=0; j<Ny; ++j) {
            for(int i=0; i<Nx; ++i) {
                file << "8 ";
                file << i + j*NxP + k*NxP*NyP << " "
                     << i + 1 + j*NxP + k*NxP*NyP << " "
                     << i + 1 + (j+1)*NxP + k*NxP*NyP << " "
                     << i + (j+1)*NxP + k*NxP*NyP << " "
                     << i + j*NxP + (k+1)*NxP*NyP << " "
                     << i + 1 + j*NxP + (k+1)*NxP*NyP << " "
                     << i + 1 + (j+1)*NxP + (k+1)*NxP*NyP << " "
                     << i + (j+1)*NxP + (k+1)*NxP*NyP << std::endl;
            }
        }
    }
    file << "CELL_TYPES " << NCells << std::endl;
    for(int k=0; k<Nz; ++k) {
        for(int j=0; j<Ny; ++j) {
            for(int i=0; i<Nx; ++i) {
                file << "12" << std::endl; // VTK_HEXAHEDRON
            }
        }
    }
    file << "CELL_DATA " << NCells << std::endl;
    file << "SCALARS procId int" << std::endl;
        file << "LOOKUP_TABLE default\n";
    for(int i=0; i<NCells; ++i) {
        // file << getCell(i+1).procId << " " << static_cast<int>(getCell(i+1).type)*10+getCell(i+1).procId << " "<<getCell(i+1).groupId << std::endl;
        file << getCell(i+1).procId << std::endl;
    }

/*     // 遍历fluidFaces和mobsFaces，如果groupName不为空，则打印
    for(int i=0; i<fluidFaces.size(); ++i) {
        for(int j=0; j<fluidFaces[i].size(); ++j) {
            for(int k=0; k<fluidFaces[i][j].size(); ++k) {
                if(fluidFaces[i][j][k].groupName.substr(0,2) == "ib") {
                    std::cout << fluidFaces[i][j][k].groupName << ":" << fluidFaces[i][j][k].globalId << "  " << fluidFaces[i][j][k].owner<<std::endl;
                }
            }
        }
    }

    for(int i=0; i<mobsFaces.size(); ++i) {
        for(int j=0; j<mobsFaces[i].size(); ++j) {
            for(int k=0; k<mobsFaces[i][j].size(); ++k) {
                if(mobsFaces[i][j][k].groupName.substr(0,2) == "ib") {
                    std::cout << mobsFaces[i][j][k].groupName << ":" << mobsFaces[i][j][k].globalId << "  " << mobsFaces[i][j][k].owner<<std::endl;
                }
            }
        }
    } */
    
/*     for(int i=0; i<fluidBoundaries.size(); ++i) {
        for(int j=0; j<fluidBoundaries[i].size(); ++j) {
            for(auto it=fluidBoundaries[i][j].begin(); it != fluidBoundaries[i][j].end(); ++it) {
                std::cout << it->first << ":" << it->second.startFace << "  " << it->second.nFaces << " from " << it->second.myProcNo << " to " << it->second.neighbProcNo << std::endl;
            }
        }
    }

    for(int i=0; i<mobsBoundaries.size(); ++i) {
        for(int j=0; j<mobsBoundaries[i].size(); ++j) {
            for(auto it=mobsBoundaries[i][j].begin(); it != mobsBoundaries[i][j].end(); ++it) {
                std::cout << it->first << ":" << it->second.startFace << "  " << it->second.nFaces << " from " << it->second.myProcNo << " to " << it->second.neighbProcNo << std::endl;
            }
        }
    } */

}

// void StructuredMesh::flatten() {
//     for (int i = 0 ;i < fluidRegions.size(); ++i) {

//         std::vector<std::vector<I64>> procMetadata(fluidRegions[i].size());
//         std::map<I64, I64> procSize;
//         std::map<I64, I64> preSize; 
//         std::vector<I64> sortedVector;
//         std::map<std::pair<I64, I64>, std::vector<I64>> interaction; 
//         sortedVector.resize(fluidRegions[i].size());
//         for (int j = 0 ; j < fluidRegions[i].size(); ++j) {
//             I64 minX = NCells, minY = NCells, minZ = NCells;
//             I64 maxX = 0, maxY = 0, maxZ = 0;            
//             procMetadata[j] = std::vector<I64>();
//             for (int k = 0; k < fluidRegions[i][j].size(); ++k) {
//                 I64 cellIndex = fluidRegions[i][j][k]; // global
//                 I64 ix, iy, iz;
//                 indexToCoords(cellIndex, ix, iy, iz);
//                 minX = std::min(ix, minX);
//                 minY = std::min(iy, minY);
//                 minZ = std::min(iz, minZ);
//                 maxX = std::max(ix, maxX);
//                 maxY = std::max(iy, maxY);
//                 maxZ = std::max(iz, maxZ);
//                 preSize[j]++;
//             }
//             procMetadata[j].insert(procMetadata[j].end(), {minX, maxX, minY, maxY, minZ, maxZ});
//             procSize[j] = ((maxZ - minZ)+1)*((maxY - minY)+1)*((maxX - minX)+1);
//             sortedVector[j] = j;
//         }
//         // auto sorted by sorted map

//         for (auto iter = procSize.begin(); iter != procSize.end(); ++iter) { //计数
//             std::cout << "i: " << i << " " << iter->first << " " << iter->second;
//             std::cout << " " << procMetadata[iter->first][0] << " " << procMetadata[iter->first][1] << " " <<
//                       procMetadata[iter->first][2] << " " << procMetadata[iter->first][3] << " " << 
//                       procMetadata[iter->first][4] << " " << procMetadata[iter->first][5] << std::endl;
            
//         }
//         for (auto iter = preSize.begin(); iter != preSize.end(); ++iter) {
//             std::cout << "i: " << iter->first << " " << iter->second << std::endl;
//         }
//         std::sort(sortedVector.begin(), sortedVector.end(),
//              [&](const I64& i, const I64& j) {return procSize[i] < procSize[j];});
//         for (int k = 0; k < sortedVector.size(); ++k) {
//             std::cout << sortedVector[k] << " " ;
//         }
//         std::cout << std::endl;

//         for (int k1 = 0; k1 < sortedVector.size(); ++k1) {
//             for (int k2 = k1 + 1; k2 < sortedVector.size(); ++k2) {
//                 auto j1 = sortedVector[k1];
//                 auto j2 = sortedVector[k2];
//                 auto p = std::make_pair(j1, j2);
//                 interaction.insert({p, std::vector<I64>()});
//                 auto left1 = std::max(procMetadata[j1][0], procMetadata[j2][0]);
//                 auto right1 = std::min(procMetadata[j1][1], procMetadata[j2][1]);
//                 if (right1 >= left1) {
//                     interaction[p].push_back(left1);
//                     interaction[p].push_back(right1);
//                 }
//                 auto left2 = std::max(procMetadata[j1][2], procMetadata[j2][2]);
//                 auto right2 = std::min(procMetadata[j1][3], procMetadata[j2][3]);
//                 if (right2 >= left2) {
//                     interaction[p].push_back(left2);                    
//                     interaction[p].push_back(right2);
//                 }
//                 auto left3 = std::max(procMetadata[j1][4], procMetadata[j2][4]);
//                 auto right3 = std::min(procMetadata[j1][5], procMetadata[j2][5]);                
//                 if (right3 >= left3) {
//                     interaction[p].push_back(left3);                    
//                     interaction[p].push_back(right3);
//                 }                                             
//             }
//         }
//         for (auto &e: interaction) {
//             if (e.second.size() < 6){
//                 e.second.clear();
//             }
//             std::cout << e.first.first << " " << e.first.second << " : value ";
//             for (auto &v: e.second) std::cout << v << " " ;
//             std::cout << std::endl;
//         }
//         // set baseline 
//         for (int k1 = 0; k1 < sortedVector.size(); ++k1) {
//             for (int k2 = k1 + 1; k2 < sortedVector.size(); ++k2) {
//                 auto j1 = sortedVector[k1];
//                 auto j2 = sortedVector[k2];    
//                 auto p = std::make_pair(j1, j2);
//                 if (interaction[p].size() == 0) 
//                     continue;
//                 auto d1 = interaction[p][1] - interaction[p][0];
//                 auto d2 = interaction[p][3] - interaction[p][2];
//                 auto d3 = interaction[p][5] - interaction[p][4];
//                 auto ciend = interaction[p][1], cjend = interaction[p][3], ckend = interaction[p][5];
//                 // if (d1 < d2 && d1 < d3) {
//                 //     ciend = interaction[p][1]/2;
//                 // } else if (d2 < d1 && d2 < d3) {
//                 //     cjend = interaction[p][3]/2;
//                 // } else if (d3 < d1 && d3 < d2) {
//                 //     ckend = interaction[p][5]/2;
//                 // }
//                 for (I64 ci = interaction[p][0]; ci <= ciend; ci++) {
//                     for (I64 cj = interaction[p][2]; cj <= cjend; cj++) {
//                         for (I64 ck = interaction[p][4]; ck <= ckend; ck++) {
//                             auto c = getCellIndex(ci, cj, ck);
//                             StCell& cell = getCell(c);
//                             if (cell.cellFlag == 0 && cell.groupId == j2) {
//                                 cell.groupId = j1;  
//                                 cell.procId = j1; 
//                                 cell.cellFlag = 1;      
//                             }                  
//                         }
//                     }
//                 }
//             }
//         }    
//     }
//     // for(size_t i=0; i<fluidRegions.size(); ++i) {
//     //     for (int j = 0; j < fluidRegions[i].size(); ++j)
//     //         fluidRegions[i][j].clear();
//     // }     
//     // for(I64 i=1; i<=NCells; ++i) {
//     //     StCell& cell = getCell(i);
//     //     if(cell.type == CellType::Fluid) {
//     //         if(cell.regionId >= 0) {
//     //             if (fluidRegions[cell.regionId].size() <= cell.groupId) {
//     //                 fluidRegions[cell.regionId].resize(cell.groupId + 1);
//     //                 fluidRegions[cell.regionId][cell.groupId] = std::vector<I64>();
//     //             }
//     //             fluidRegions[cell.regionId][cell.groupId].push_back(i);
//     //         } else {
//     //             std::cerr << "Fluid cell without region ID." << std::endl;
//     //         }
//     //     } else if(cell.type == CellType::Obstacle) {
//     //         if(cell.regionId >= 0) {
//     //             mobsRegions[cell.regionId][0].push_back(i);
//     //         } else {
//     //             std::cerr << "Obstacle cell without region ID." << std::endl;
//     //         }
//     //     } else {
//     //         std::cerr << "TODO: Hole cell." << std::endl;
//     //     }
//     // }
//     // for (int i = 0 ;i < fluidRegions.size(); ++i) {

//     //     std::vector<std::vector<I64>> procMetadata(fluidRegions[i].size());
//     //     std::map<I64, I64> procSize;
//     //     std::map<I64, I64> preSize; 
//     //     std::vector<I64> sortedVector;
//     //     std::map<std::pair<I64, I64>, std::vector<I64>> interaction; 
//     //     sortedVector.resize(fluidRegions[i].size());
//     //     for (int j = 0 ; j < fluidRegions[i].size(); ++j) {
//     //         I64 minX = NCells, minY = NCells, minZ = NCells;
//     //         I64 maxX = 0, maxY = 0, maxZ = 0;            
//     //         procMetadata[j] = std::vector<I64>();
//     //         for (int k = 0; k < fluidRegions[i][j].size(); ++k) {
//     //             I64 cellIndex = fluidRegions[i][j][k]; // global
//     //             I64 ix, iy, iz;
//     //             indexToCoords(cellIndex, ix, iy, iz);
//     //             minX = std::min(ix, minX);
//     //             minY = std::min(iy, minY);
//     //             minZ = std::min(iz, minZ);
//     //             maxX = std::max(ix, maxX);
//     //             maxY = std::max(iy, maxY);
//     //             maxZ = std::max(iz, maxZ);
//     //             preSize[j]++;
//     //         }
//     //         procMetadata[j].insert(procMetadata[j].end(), {minX, maxX, minY, maxY, minZ, maxZ});
//     //         procSize[j] = ((maxZ - minZ)+1)*((maxY - minY)+1)*((maxX - minX)+1);
//     //         sortedVector[j] = j;
//     //     }
//     //     // auto sorted by sorted map

//     //     for (auto iter = procSize.begin(); iter != procSize.end(); ++iter) { //计数
//     //         std::cout << "i: " << i << " " << iter->first << " " << iter->second;
//     //         std::cout << " " << procMetadata[iter->first][0] << " " << procMetadata[iter->first][1] << " " <<
//     //                   procMetadata[iter->first][2] << " " << procMetadata[iter->first][3] << " " << 
//     //                   procMetadata[iter->first][4] << " " << procMetadata[iter->first][5] << std::endl;
            
//     //     }
//     //     for (auto iter = preSize.begin(); iter != preSize.end(); ++iter) {
//     //         std::cout << "i: " << iter->first << " " << iter->second << std::endl;
//     //     }
//     //     std::sort(sortedVector.begin(), sortedVector.end(),
//     //          [&](const I64& i, const I64& j) {return procSize[i] < procSize[j];});
//     //     for (int k = 0; k < sortedVector.size(); ++k) {
//     //         std::cout << sortedVector[k] << " " ;
//     //     }
//     //     std::cout << std::endl;

//     //     for (int k1 = 0; k1 < sortedVector.size(); ++k1) {
//     //         for (int k2 = k1 + 1; k2 < sortedVector.size(); ++k2) {
//     //             auto j1 = sortedVector[k1];
//     //             auto j2 = sortedVector[k2];
//     //             auto p = std::make_pair(j1, j2);
//     //             interaction.insert({p, std::vector<I64>()});
//     //             auto left1 = std::max(procMetadata[j1][0], procMetadata[j2][0]);
//     //             auto right1 = std::min(procMetadata[j1][1], procMetadata[j2][1]);
//     //             if (right1 >= left1) {
//     //                 interaction[p].push_back(left1);
//     //                 interaction[p].push_back(right1);
//     //             }
//     //             auto left2 = std::max(procMetadata[j1][2], procMetadata[j2][2]);
//     //             auto right2 = std::min(procMetadata[j1][3], procMetadata[j2][3]);
//     //             if (right2 >= left2) {
//     //                 interaction[p].push_back(left2);                    
//     //                 interaction[p].push_back(right2);
//     //             }
//     //             auto left3 = std::max(procMetadata[j1][4], procMetadata[j2][4]);
//     //             auto right3 = std::min(procMetadata[j1][5], procMetadata[j2][5]);                
//     //             if (right3 >= left3) {
//     //                 interaction[p].push_back(left3);                    
//     //                 interaction[p].push_back(right3);
//     //             }                                             
//     //         }
//     //     }
//     //     for (auto &e: interaction) {
//     //         if (e.second.size() < 6){
//     //             e.second.clear();
//     //         }
//     //         std::cout << e.first.first << " " << e.first.second << " : value ";
//     //         for (auto &v: e.second) std::cout << v << " " ;
//     //         std::cout << std::endl;
//     //     }
//     //     // set baseline 
//     //     for (int k1 = 0; k1 < sortedVector.size(); ++k1) {
//     //          for (int k2 = k1 + 1; k2 < sortedVector.size(); ++k2) {
//     //             auto j1 = sortedVector[k1];
//     //             auto j2 = sortedVector[k2];    
//     //             auto p = std::make_pair(j1, j2);
//     //             if (interaction[p].size() == 0) 
//     //                 continue;
//     //             // auto d1 = interaction[p][1] - interaction[p][0];
//     //             // auto d2 = interaction[p][3] - interaction[p][2];
//     //             // auto d3 = interaction[p][5] - interaction[p][4];
//     //             auto ciend = interaction[p][1], cjend = interaction[p][3], ckend = interaction[p][5];
//     //             for (I64 ci = interaction[p][0]; ci <= ciend; ci++) {
//     //                 for (I64 cj = interaction[p][2]; cj <= cjend; cj++) {
//     //                     for (I64 ck = interaction[p][4]; ck <= ckend; ck++) {
//     //                         auto c = getCellIndex(ci, cj, ck);
//     //                         StCell& cell = getCell(c);
//     //                         if (cell.cellFlag == 1 || cell.cellFlag == 0 && cell.groupId == j2) {
//     //                             cell.groupId = j1;  
//     //                             cell.procId = j1; 
//     //                             cell.cellFlag = 2;      
//     //                         }                  
//     //                     }
//     //                 }
//     //             }
//     //         }
//     //     }    
//     // }   
//     // for(size_t i=0; i<fluidRegions.size(); ++i) {
//     //     for (int j = 0; j < fluidRegions[i].size(); ++j)
//     //         fluidRegions[i][j].clear();
//     // }     
//     // for(I64 i=1; i<=NCells; ++i) {
//     //     StCell& cell = getCell(i);
//     //     if(cell.type == CellType::Fluid) {
//     //         if(cell.regionId >= 0) {
//     //             if (fluidRegions[cell.regionId].size() <= cell.groupId) {
//     //                 fluidRegions[cell.regionId].resize(cell.groupId + 1);
//     //                 fluidRegions[cell.regionId][cell.groupId] = std::vector<I64>();
//     //             }
//     //             fluidRegions[cell.regionId][cell.groupId].push_back(i);
//     //         } else {
//     //             std::cerr << "Fluid cell without region ID." << std::endl;
//     //         }
//     //     } else if(cell.type == CellType::Obstacle) {
//     //         if(cell.regionId >= 0) {
//     //             mobsRegions[cell.regionId][0].push_back(i);
//     //         } else {
//     //             std::cerr << "Obstacle cell without region ID." << std::endl;
//     //         }
//     //     } else {
//     //         std::cerr << "TODO: Hole cell." << std::endl;
//     //     }
//     // }      
// }

void StructuredMesh::writeProcFile(std::string prefix, int nParts) {
    // 写入处理器文件的实现
    // 这里可以使用类似于读取网格文件的方式来读取配置文件
    // 例如，使用std::ofstream打开文件并解析内容

    for(int i=0; i<nParts; i++)
    {
        // 创建文件夹
        std::string dirName = "processor" + std::to_string(i);
        std::string command = "mkdir -p " + dirName;
        system(command.c_str());
        // 创建文件
        std::string cellfilename = dirName + "/cells.txt";
        std::string facefilename = dirName + "/faces.txt";
        std::string boundaryfilename = dirName + "/boundary.txt";
        std::ofstream cellfile(cellfilename);
        std::ofstream facefile(facefilename);
        std::ofstream boundaryfile(boundaryfilename);
        if (!cellfile.is_open() || !facefile.is_open() || !boundaryfile.is_open()) {
            std::cerr << "无法打开文件！" << std::endl;
            return;
        }

        int proc =0;
        for(int j=0; j<fluidRegions.size(); j++)
        {
            if(i>=proc+fluidRegions[j].size())
            {
                proc += fluidRegions[j].size();
                continue;
            }

            cellfile << "meta fluid " << fluidRegions[j][i-proc].size() << " " << j << " " << i-proc << std::endl;
            for(int k=0; k<fluidRegions[j][i-proc].size(); ++k)
            {
                cellfile << fluidRegions[j][i-proc][k] << " ";
            }
            cellfile << std::endl;

            facefile << "meta fluid " << fluidFaces[j][i-proc].size() << " " << j << " " << i-proc << std::endl;
            for(int k=0; k<fluidFaces[j][i-proc].size(); ++k)
            {
                facefile << fluidFaces[j][i-proc][k].globalId << " "<< 
                    fluidFaces[j][i-proc][k].owner << " " << fluidFaces[j][i-proc][k].neighbour << std::endl;
            }
            facefile << std::endl;

            boundaryfile << "meta fluid " << fluidBoundaries[j][i-proc].size() << " " << j << " " << i-proc << std::endl;
            for(auto it=fluidBoundaries[j][i-proc].begin(); it != fluidBoundaries[j][i-proc].end(); ++it)
            {
                boundaryfile << it->first << " " << it->second.startFace << " " << it->second.nFaces << std::endl;
            }
            boundaryfile << std::endl;

            break;
        }

        proc = 0;
        for(int j=0; j<mobsRegions.size(); j++)
        {
            if(i>=proc+mobsRegions[j].size())
            {
                proc += mobsRegions[j].size();
                continue;
            }

            cellfile << "meta mobs " << mobsRegions[j][i-proc].size() << " " << j << " " << i-proc << std::endl;
            for(int k=0; k<mobsRegions[j][i-proc].size(); ++k)
            {
                cellfile << mobsRegions[j][i-proc][k] << " ";
            }
            cellfile << std::endl;

            facefile << "meta mobs " << mobsFaces[j][i-proc].size() << " " << j << " " << i-proc << std::endl;
            for(int k=0; k<mobsFaces[j][i-proc].size(); ++k)
            {
                facefile << mobsFaces[j][i-proc][k].globalId << " "
                    << mobsFaces[j][i-proc][k].owner << " " << mobsFaces[j][i-proc][k].neighbour << std::endl;
            }
            facefile << std::endl;

            boundaryfile << "meta mobs " << mobsBoundaries[j][i-proc].size() << " " << j << " " << i-proc << std::endl;
            for(auto it=mobsBoundaries[j][i-proc].begin(); it != mobsBoundaries[j][i-proc].end(); ++it)
            {
                boundaryfile << it->first << " " << it->second.startFace << " " << it->second.nFaces << std::endl;
            }
            boundaryfile << std::endl;

            break;
        }
        
        cellfile.close();
        facefile.close();
        boundaryfile.close();
            
            
    }
}

CYCASMesh::CYCASMesh(std::string filename, int iProc)
{
    // 初始化网格参数
    readMesh(filename);
    setDims();

}

void CYCASMesh::readMesh(std::string filename) {
    // 读取网格文件的实现
    // 这里可以使用类似于读取网格文件的方式来读取配置文件
    // 例如，使用std::ifstream打开文件并解析内容

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return;
    }

    // dummy example
    /*
    xpoints = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]
    ypoints = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]
    zpoints = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]
    mobs = [0.4 0.6 0.2 0.5 0.4 0.6]
    wall = [0.5 0.5 0.2 0.3 0.4 0.5]
    */
    std::ostringstream ss;
    ss << file.rdbuf();  // 读取整个文件内容到 ss
    std::string content = ss.str();
    auto xgridPos = content.find("xgrid");
    auto ygridPos = content.find("ygrid"), ytemp = ygridPos;
    auto zgridPos = content.find("zgrid");
    auto mobsPos = content.find("mobs");
    auto wallsPos = content.find("walls");
    auto holesPos = content.find("holes");
    if (ygridPos == std::string::npos) {
        ygridPos = content.find("nyr(1)=");
    }
    std::string sliceString, line, name;
    char eq, spliter;
    double value1;
    int value2;
    auto xstart = content.substr(xgridPos).find_first_of('=');
    std::istringstream iss1(content.substr(xstart + 1, ygridPos - xstart));
    while(std::getline(iss1, line)) {
        sliceString.clear();
        sliceInputString(line, sliceString);
        std::istringstream lineiss(sliceString);
        while (lineiss >> value1 >> spliter)
            xPoints.push_back(value1);
    }
    if (ygridPos != ytemp)
        ygridPos = ytemp;
    std::cout << "size: " << xPoints.size() << std::endl;
    if (ygridPos == std::string::npos) {
        auto cylindricalPos = content.find("nyr(1)=");
        auto commaPos = content.find(",", cylindricalPos);
        std::string num_str = content.substr(cylindricalPos + 7, commaPos - (cylindricalPos + 7));
        auto num = std::stoi(num_str);

        cylindricalPos = content.find("yl(1)=");
        commaPos = content.find(",", cylindricalPos);
        num_str = content.substr(cylindricalPos + 6, commaPos - (cylindricalPos + 6));
        auto start_angle = std::stoi(num_str);

        cylindricalPos = content.find("yl(2)=");
        commaPos = content.find(",", cylindricalPos);
        num_str = content.substr(cylindricalPos + 6, commaPos - (cylindricalPos + 6));
        auto end_angle = std::stoi(num_str);
        for (int i = 0; i < num; ++i) {
            yPoints.push_back((end_angle - start_angle) / num * i);
        }
        std::cout << "size: " << yPoints.size() << std::endl;
        // perioFlag = end_angle - start_angle;
    } else {
        auto ystart = content.substr(ygridPos).find_first_of('=');
        std::istringstream iss2(content.substr(ygridPos).substr(ystart + 1, zgridPos - ygridPos));
        while(std::getline(iss2, line)) {
            sliceString.clear();
            sliceInputString(line, sliceString);
            std::istringstream lineiss(sliceString);
            while (lineiss >> value1 >> spliter)
                yPoints.push_back(value1);
        }          
    }
 
    auto zstart = content.substr(zgridPos).find_first_of('=');
    std::istringstream iss3(content.substr(zgridPos).substr(zstart + 1, mobsPos - zgridPos));
    while(std::getline(iss3, line)) {
        sliceString.clear();
        sliceInputString(line, sliceString);
        std::istringstream lineiss(sliceString);
        while (lineiss >> value1 >> spliter)
            zPoints.push_back(value1);
    }    

    // auto mobsstart = content.substr(mobsPos).find_first_of('=');
    // std::istringstream iss4(content.substr(mobsPos).substr(mobsstart + 1, wallsPos - mobsPos));
    // while(std::getline(iss4, line)) {
    //     std::vector<int> mobs;
    //     sliceString.clear();
    //     sliceInputString(line, sliceString);
    //     std::istringstream lineiss(sliceString);
    //     while (lineiss >> value2 >> spliter)
    //         mobs.push_back(value2);
    //     if(mobs.size() != 0)
    //         mobsPoints.push_back(mobs);
    // } 

    file.close();

}

CYCASMesh::~CYCASMesh() 
{

}

void CYCASMesh::setDims() {
    NxP = xPoints.size();
    NyP = yPoints.size();
    NzP = zPoints.size();
    Nx = NxP - 1;
    Ny = NyP - 1;
    Nz = NzP - 1;
    NCells = Nx * Ny * Nz;
    NFaces = 2 * (Nx * Ny + Ny * Nz + Nz * Nx);
    NPoints = NxP * NyP * NzP;
}

void CYCASMesh::readCells(int iProc)
{
    // 读取单元格的实现
    // 这里可以使用类似于读取网格文件的方式来读取配置文件
    // 例如，使用std::ifstream打开文件并解析内容

    std::string filename = "processor" + std::to_string(iProc) + "/cells.txt";
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // 检查line是否以“meta”开头
        if (line.substr(0, 4) == "meta") {
            std::istringstream iss(line);
            std::string meta, type;
            int nCells, regionId, groupId;
            iss >> meta >> type >> nCells >> regionId >> groupId;
            std::getline(file, line); // 读取下一行，包含单元格索引
            std::istringstream iss2(line);
            std::vector<I64> regionCells;
            for (int i = 0; i < nCells; ++i) {
                I64 cellIndex;
                iss2 >> cellIndex;
                regionCells.push_back(cellIndex);
            }
            if (type == "fluid") {
                fluidRegions[regionId] = regionCells;
            } else if (type == "mobs") {
                mobsRegions[regionId] = regionCells;
            }
        }
        
    }

    file.close();
}

void CYCASMesh::readFaces(int iProc)
{
    // 读取面信息的实现
    // 这里可以使用类似于读取网格文件的方式来读取配置文件
    // 例如，使用std::ifstream打开文件并解析内容

    std::string filename = "processor" + std::to_string(iProc) + "/faces.txt";
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // 检查line是否以“meta”开头
        if (line.substr(0, 4) == "meta") {
            std::istringstream iss(line);
            std::string meta, type;
            int nFaces, regionId, groupId;
            iss >> meta >> type >> nFaces >> regionId >> groupId;
            
            std::vector<StFace> regionFaces;
            for (int i = 0; i < nFaces; ++i) {
                std::getline(file, line); // 读取下一行，包含单元格索引
                std::istringstream iss2(line);
                StFace face;
                iss2 >> face.globalId >> face.owner >> face.neighbour;
                regionFaces.push_back(face);
            }
            if (type == "fluid") {
                fluidFaces[regionId] = regionFaces;
            } else if (type == "mobs") {
                mobsFaces[regionId] = regionFaces;
            }
        }
        
    }

    file.close();
}

void CYCASMesh::readBoundaries(int iProc)
{
    // 读取边界信息的实现
    // 这里可以使用类似于读取网格文件的方式来读取配置文件
    // 例如，使用std::ifstream打开文件并解析内容

    std::string filename = "processor" + std::to_string(iProc) + "/boundary.txt";
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // 检查line是否以“meta”开头
        if (line.substr(0, 4) == "meta") {
            std::istringstream iss(line);
            std::string meta, type;
            int nBds, regionId, groupId;
            iss >> meta >> type >> nBds >> regionId >> groupId;
            
            for(int i = 0; i < nBds; ++i) {
                std::getline(file, line); // 读取下一行，包含单元格索引
                std::istringstream iss2(line);
                std::string bdName;
                I64 startFace, nFaces;
                iss2 >> bdName >> startFace >> nFaces;
                BoundaryFaces bdFaces;
                bdFaces.name = bdName;
                bdFaces.startFace = startFace;
                bdFaces.nFaces = nFaces;

                if (type == "fluid") {
                    if(fluidBoundaries.find(regionId) == fluidBoundaries.end()) {
                        fluidBoundaries[regionId] = std::map<std::string, BoundaryFaces>();
                    }
                    fluidBoundaries[regionId][bdName] = bdFaces; // 把bdFaces加入到类成员中
                } else if (type == "mobs") {
                    if(mobsBoundaries.find(regionId) == mobsBoundaries.end()) {
                        mobsBoundaries[regionId] = std::map<std::string, BoundaryFaces>();
                    }
                    mobsBoundaries[regionId][bdName] = bdFaces; // 把bdFaces加入到类成员中
                }
            }


        }
        
    }

    file.close();
}

void CYCASMesh::printInfo()
{
    // 打印fluidRegions和mobsRegions, fluidFaces和mobsFaces, fluidBoundaries和mobsBoundaries的信息
    std::cout << "Fluid Regions:" << std::endl;
    for (const auto& region : fluidRegions) {
        std::cout << "Region ID: " << region.first << ", Cells: ";
        for (const auto& cell : region.second) {
            std::cout << cell << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Mobs Regions:" << std::endl;
    for (const auto& region : mobsRegions) {
        std::cout << "Region ID: " << region.first << ", Cells: ";
        for (const auto& cell : region.second) {
            std::cout << cell << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Fluid Faces:" << std::endl;
    for (const auto& region : fluidFaces) {
        std::cout << "Region ID: " << region.first << ", Faces: ";
        for (const auto& face : region.second) {
            std::cout << face.globalId << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Mobs Faces:" << std::endl;
    for (const auto& region : mobsFaces) {
        std::cout << "Region ID: " << region.first << ", Faces: ";
        for (const auto& face : region.second) {
            std::cout << face.globalId << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Fluid Boundaries:" << std::endl;
    for (const auto& region : fluidBoundaries) {
        std::cout << "Region ID: " << region.first << ", Boundaries: ";
        for (const auto& boundary : region.second) {
            std::cout << boundary.first << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Mobs Boundaries:" << std::endl;
    for (const auto& region : mobsBoundaries) {
        std::cout << "Region ID: " << region.first << ", Boundaries: ";
        for (const auto& boundary : region.second) {
            std::cout << boundary.first << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Number of Cells: " << NCells << std::endl;
}

