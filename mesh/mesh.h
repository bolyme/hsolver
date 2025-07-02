#ifndef MESH_H
#define MESH_H

#include <vector>
#include <map>
#include <string>

#include "../defs.h"

enum class CellType { 
    Fluid, 
    Obstacle, 
    Hole
};

enum class FaceType { 
    Internal, 
    Boundary, 
    Processor,
    Wall,
    Hole,
    Interface
};

struct Point {
    double x, y, z;
};

struct Face {
    FaceType type = FaceType::Internal; // 面的类型
    I64 owner = -1; // 面的邻居cell ID, -1表示没有邻居
    I64 neighbour = -1; // 面的邻居cell ID, -1表示没有邻居
    I64 globalId = -1; // 面的全局ID, 0开始
    int fromProcId = -1; // 面的处理器ID, 0开始
    int toProcId = -1; // 面的邻居处理器ID, 0开始
    std::vector<int> vertices;
};

struct Cell {
    CellType type; // cell的类型
    int regionId; // cell的区域ID, 0开始
    int groupId; // 通信组ID, 0开始
    int procId; // cell的处理器ID， 0开始
    std::vector<int> vertices; // cell的顶点索引
};

struct StFace
{
    I64 owner = 0; // 面的邻居cell ID, -1表示没有邻居
    I64 neighbour = 0; // 面的邻居cell ID, -1表示没有邻居
    I64 globalId = 0; // 面的全局ID, 0开始
    FaceType type = FaceType::Internal; // 面的类型
    std::string groupName=""; // 面的组名
};

struct StCell
{
    I64 globalId = 0; // cell的全局ID, 0开始
    int regionId = 0; // cell的区域ID, 0开始
    int groupId = 0; // region通信组中的processor id, 0开始
    int procId = 0; // 全局通信组中的processor id, 0开始
    CellType type = CellType::Fluid; // cell的类型
    int cellFlag = 0;
};



struct Cells {
    std::vector<Cell> cells; // cell的列表

    Cell& operator[](int index) {
        return cells[index-1]; // index从1开始
    }

    Cell& operator[](I64 index) {
        return cells[index-1]; // index从1开始
    }

    int size() const {
        return cells.size(); // 返回cell的数量
    }
};

struct ProcFaces {
    int startFace;
    int nFaces;
};

struct BoundaryFaces
{
    std::string name;
    std::string type;
    int nFaces;
    int startFace;
    int myProcNo;
    int neighbProcNo;

    BoundaryFaces(std::string nm, std::string t, int n, int s, int mp, int np) : name(nm), type(t), nFaces(n), startFace(s), myProcNo(mp), neighbProcNo(np) {}
    BoundaryFaces() : name(""), type(""), nFaces(0), startFace(-1), myProcNo(-1), neighbProcNo(-1) {}
};

// OpenFOAM网格数据结构
struct OFMesh {
    std::vector<Point> points;
    std::vector<Face> faces;
    std::vector<int> owner;
    std::vector<int> neighbour;
    std::map<std::string, BoundaryFaces> boundaries;
};



typedef std::vector<I64> ProcCells;
typedef std::vector<ProcCells> RegionProcs;
typedef std::vector<RegionProcs> Regions;


struct StructuredMesh
{
    
    I64 NxP, NyP, NzP; // 网格尺寸
    I64 Nx, Ny, Nz; // 体网格尺寸
    I64 NxF, NyF, NzF; // x, y, z方向面网格数量
    I64 NxNeigh, NyNeigh, NzNeigh; // 邻居网格尺寸
    I64 NCells; // 网格单元数
    I64 NFaces; // 网格面数
    I64 NPoints; // 网格点数
    I64 perioFlag = 0;

    // 从文件中读取网格数据
    std::vector<double> xPoints; // 网格点坐标
    std::vector<double> yPoints; // 网格点坐标
    std::vector<double> zPoints; // 网格点坐标
    std::vector<std::vector<int>> mobsPoints; // 障碍物点坐标
    std::vector<std::vector<int>> wallPoints; // 墙面点坐标
    std::vector<std::vector<int>> holePoints; // 墙面点坐标


    std::vector<StCell> domain; // 整个模拟区域
    //std::vector<Face> domainFaces; // 整个模拟区域的面

    StCell& getCell(I64 i, I64 j, I64 k) {
        return domain[getCellIndex(i, j, k)-1];
    }

    StCell& getCell(I64 index) {
        return domain[index-1]; // index从1开始
    }

    Regions fluidRegions; // 流体区域
    Regions mobsRegions; // 障碍物区域
    std::vector<std::vector<std::vector<StFace>>> fluidFaces; // 流体区域面
    std::vector<std::vector<std::vector<StFace>>> mobsFaces; // 障碍物区域面
    std::vector<std::vector<std::map<std::string, BoundaryFaces>>> fluidBoundaries; // 边界面
    std::vector<std::vector<std::map<std::string, BoundaryFaces>>> mobsBoundaries; // 边界面
    
    //std::vector<CellType> cellTypes; // 网格单元 idx(i,j,k) = (i-ONE) + (j-ONE)*Nx + (k-ONE)*Nx*Ny
    std::vector<FaceType> faceTypes; // 网格面 xidx(i,j,k) = (i-ONE) + (j-ONE)*NxP + (k-ONE)*NxP*Ny
                                 // yidx(i,j,k) = NxF + (i-ONE) + (j-ONE)*Nx + (k-ONE)*NyP*Nx
                                 // zidx(i,j,k) = NxF + NyF + (i-ONE) + (j-ONE)*Nx + (k-ONE)*Ny*Nx
    int nFluidProcs, nMobsProcs;

/*     std::vector<std::vector<I64>> fRegions; // 流体区域
    std::vector<std::vector<I64>> mRegions; // obstacle区域 */

/*     std::vector<std::vector<Point>> fluidRegionPoints; // 流体区域点坐标
    std::vector<std::vector<Point>> mobsRegionPoints; // obstacle区域点坐标 */

/*     std::vector<std::vector<Face>> fluidRegionFaces; // 流体区域面
    std::vector<std::vector<Face>> mobsRegionFaces; // obstacle区域面 */

/*     std::vector<std::vector<Cell>> fluidRegionCells; // 流体区域单元
    std::vector<std::vector<Cell>> mobsRegionCells; // obstacle区域单元 */

/*     std::vector<std::vector<std::vector<I64>>> fluidRegionOwners; // 流体区域单元owner
    std::vector<std::vector<std::vector<I64>>> mobsRegionOwners; // obstacle区域单元owner
    std::vector<std::vector<std::vector<I64>>> fluidRegionNeighbours; // 流体区域单元neighbour
    std::vector<std::vector<std::vector<I64>>> mobsRegionNeighbours; // obstacle区域单元neighbour */

    
    // constructor
    StructuredMesh(std::string filename);

    // destructor
    ~StructuredMesh();

    void getConnectedRegions(std::vector<std::vector<I64>>& fRegions,
        std::vector<std::vector<I64>>& mRegions);

    void decomposeMesh(std::vector<std::vector<I64>>& fRegions,
        std::vector<std::vector<I64>>& mRegions, I32 nParts);

    void flatten();

    void writeProcFile(std::string filename, int nParts);  

    void createMesh();

    void writeMesh(std::string prefix);

    // cell index start from 1
    inline I64 getCellIndex(I64 i, I64 j, I64 k) const {
        return (i-ONE) + (j-ONE)*Nx + (k-ONE)*Nx*Ny + ONE; // index from 1
    }

    inline void getCellIndex(I64 i, I64 j, I64 k, I64& index) const {
        index = (i-ONE) + (j-ONE)*Nx + (k-ONE)*Nx*Ny + ONE; // index from 1
    }

    inline I64 getFaceIndexX(I64 i, I64 j, I64 k) const {
        return (i-ONE) + (j-ONE)*NxP + (k-ONE)*NxP*Ny;
    }

    inline I64 getFaceIndexY(I64 i, I64 j, I64 k) const {
        return NxF + (i-ONE) + (j-ONE)*Nx + (k-ONE)*NyP*Nx;
    }

    inline I64 getFaceIndexZ(I64 i, I64 j, I64 k) const {
        return NxF + NyF + (i-ONE) + (j-ONE)*Nx + (k-ONE)*Ny*Nx;
    }

    inline std::vector<I64> getNeighbours(I64 idx) const {
        std::vector<I64> neighbours(6, -1);
        if(idx -NxNeigh > 0 && idx - NxNeigh <= NCells) {
            neighbours[0] = idx - NxNeigh; // x-1
        }
        if(idx + NxNeigh > 0 && idx + NxNeigh <= NCells) {
            neighbours[1] = idx + NxNeigh; // x+1
        }
        if(idx - NyNeigh > 0 && idx - NyNeigh <= NCells) {
            neighbours[2] = idx - NyNeigh; // y-1
        }
        if(idx + NyNeigh > 0 && idx + NyNeigh <= NCells) {
            neighbours[3] = idx + NyNeigh; // y+1
        }
        if(idx - NzNeigh > 0 && idx - NzNeigh <= NCells) {
            neighbours[4] = idx - NzNeigh; // z-1
        }
        if(idx + NzNeigh > 0 && idx + NzNeigh <= NCells) {
            neighbours[5] = idx + NzNeigh; // z+1
        }
        return neighbours;
    }

    inline void getFaceIndexes(I64 i, I64 j, I64 k, std::vector<I64>& faceIndexes) const {
        faceIndexes.resize(6);
        faceIndexes[0] = getFaceIndexX(i, j, k); // x-1
        faceIndexes[1] = getFaceIndexX(i+1, j, k); // x+1
        faceIndexes[2] = getFaceIndexY(i, j, k); // y-1
        faceIndexes[3] = getFaceIndexY(i, j+1, k); // y+1
        faceIndexes[4] = getFaceIndexZ(i, j, k); // z-1
        faceIndexes[5] = getFaceIndexZ(i, j, k+1); // z+1
    }

    inline void getPointIndexes(I64 i, I64 j, I64 k, std::vector<I64>& pointIndexes) const {
        pointIndexes.resize(8);
        pointIndexes[0] = (i-ONE) + (j-ONE)*NxP + (k-ONE)*NxP*NyP; // x-1 y-1 z-1
        pointIndexes[1] = (i) + (j-ONE)*NxP + (k-ONE)*NxP*NyP; // x+1 y-1 z-1
        pointIndexes[2] = (i) + (j)*NxP + (k-ONE)*NxP*NyP; // x+1 y+1 z-1
        pointIndexes[3] = (i-ONE) + (j)*NxP + (k-ONE)*NxP*NyP; // x-1 y+1 z-1
        pointIndexes[4] = (i-ONE) + (j-ONE)*NxP + (k)*NxP*NyP; // x-1 y-1 z+1
        pointIndexes[5] = (i) + (j-ONE)*NxP + (k)*NxP*NyP; // x+1 y-1 z+1
        pointIndexes[6] = (i) + (j)*NxP + (k)*NxP*NyP; // x+1 y+1 z+1
        pointIndexes[7] = (i-ONE) + (j)*NxP + (k)*NxP*NyP; // x-1 y+1 z+1
    }

    inline void getPointCoords(I64 i, I64 j, I64 k, std::vector<Point>& coords) const {
        coords.resize(8);
        coords[0] = {xPoints[i-1], yPoints[j-1], zPoints[k-1]}; // x-1 y-1 z-1
        coords[1] = {xPoints[i], yPoints[j-1], zPoints[k-1]}; // x+1 y-1 z-1
        coords[2] = {xPoints[i], yPoints[j], zPoints[k-1]}; // x+1 y+1 z-1
        coords[3] = {xPoints[i-1], yPoints[j], zPoints[k-1]}; // x-1 y+1 z-1
        coords[4] = {xPoints[i-1], yPoints[j-1], zPoints[k]}; // x-1 y-1 z+1
        coords[5] = {xPoints[i], yPoints[j-1], zPoints[k]}; // x+1 y-1 z+1
        coords[6] = {xPoints[i], yPoints[j], zPoints[k]}; // x+1 y+1 z+1
        coords[7] = {xPoints[i-1], yPoints[j], zPoints[k]}; // x-1 y+1 z+1
    }


    void indexToCoords(I64 index, I64& i, I64& j, I64& k) const {
        index -= ONE; // index from 1
        k = index / (Nx * Ny);
        j = (index - k * Nx * Ny) / Nx;
        i = index - k * Nx * Ny - j * Nx;
        i += ONE;
        j += ONE;
        k += ONE;
    }


    
    
private:
    void readMesh(std::string filename);

    void setDims(); // 设置网格维度

    void initCellsType();

    void initFacesType();

    void BFS(I64 start, std::vector<std::vector<I64>>& cellAdjs,
             std::vector<bool>& visited, std::vector<I64>& region);

};

struct CYCASMesh
{
    I64 NxP, NyP, NzP; // 网格尺寸
    I64 Nx, Ny, Nz; // 体网格尺寸
    I64 NxF, NyF, NzF; // x, y, z方向面网格数量
    I64 NxNeigh, NyNeigh, NzNeigh; // 邻居网格尺寸
    I64 NCells; // 网格单元数
    I64 NFaces; // 网格面数
    I64 NPoints; // 网格点数

    I64 LocalCells; // 本地cell数量
    I64 LocalFaces; // 本地face数量
    I64 LocalPoints; // 本地point数量

    // 从文件中读取网格数据
    std::vector<double> xPoints; // 网格点坐标
    std::vector<double> yPoints; // 网格点坐标
    std::vector<double> zPoints; // 网格点坐标

    std::map<int, std::vector<I64>> fluidRegions; // 流体区域
    std::map<int, std::vector<I64>> mobsRegions; // 障碍物区域

    std::map<int, std::vector<StFace>> fluidFaces; // 流体区域面
    std::map<int, std::vector<StFace>> mobsFaces; // 障碍物区域面

    std::map<int, std::map<std::string, BoundaryFaces>> fluidBoundaries; // 流体区域边界面
    std::map<int, std::map<std::string, BoundaryFaces>> mobsBoundaries; // 障碍物区域边界面

    CYCASMesh(std::string filename, int iProc);
    ~CYCASMesh();
    void readMesh(std::string filename); // 读取网格文件
    void setDims(); // 设置网格维度
    void readCells(int iProc); // 读取cell数据
    void readFaces(int iProc); // 读取face数据
    void readBoundaries(int iProc); // 读取边界数据

    void printInfo(); // 打印网格信息
};


#endif // MESH_H
