#include <iostream>

#include "mesh/mesh.h"

int main(int argc, char** argv) {
    
    std::string filename = "example13.in";
    int nParts = 4; // 处理器数

    StructuredMesh mesh(filename);

    std::vector<std::vector<I64>> fRegions; // 流体区域
    std::vector<std::vector<I64>> mRegions; // obstacle区域
    std::cout << "Info: getting connected regions." << std::endl;
    mesh.getConnectedRegions(fRegions, mRegions); // 获取单元之间的连接关系
    
    std::cout << "Info: decomposing mesh." << std::endl;
    mesh.decomposeMesh(fRegions, mRegions, nParts); // 划分子域
    
    mesh.writeProcFile("test", nParts); // 写入proc文件

    mesh.writeMesh("output"); // 写入网格文件

    std::cout << "Info: fuildRegions: " << mesh.fluidRegions.size() << std::endl;
    std::cout << "Info: mobsRegions: " << mesh.mobsRegions.size() << std::endl;


    //mesh.decomposeMesh
    //mesh.writeMesh("output");

    return 0;
}
