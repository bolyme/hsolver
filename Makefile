# 编译器配置
CXX      = g++
MPICXX  = mpicxx
MPICXXFLAGS = -g -O3 -I. -I../ -I/usr/include/x86_64-linux-gnu/mpi -I/usr/include/scotch /home/dennis/hypre-2.27.0/src/hypre/include
CXXFLAGS = -g -O3 -I. -I../ -I/usr/include/x86_64-linux-gnu/mpi -I/usr/include/scotch -I/home/dennis/hypre-2.27.0/src/hypre/include
MPILDFLAGS  = -g -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/scotch -L/home/dennis/hypre-2.27.0/src/hypre/lib
MPILIBS   = -lscotch -lscotcherr -lHYPRE

# 目标可执行文件
EXEC = main
DECOMP = decompose

# 源文件和对象文件
SRCS = main.cpp mesh/mesh.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = mesh/mesh.h distributed_mesh.h fluid_solver.h mobs_solver.h solver.h

SRCS2 = decompose.cpp mesh/mesh.cpp
OBJS2 = $(SRCS2:.cpp=.o)
DEPS2 = mesh/mesh.h

# 默认构建目标
all: $(EXEC) $(DECOMP)

# 链接可执行文件
$(EXEC): $(OBJS)
	$(MPICXX) $(MPILDFLAGS) $^ -o $@ $(MPILIBS)

# 生成分解器可执行文件
$(DECOMP): $(OBJS2)
	$(MPICXX) $(MPILDFLAGS) $^ -o $@ $(MPILIBS)

# 通用编译规则
%.o: %.cpp $(DEPS) $(DEPS2)
	$(CXX)  -c $< -o $@ $(CXXFLAGS)

# 清理生成文件
clean:
	rm -f $(OBJS) $(EXEC) $(OBJS2) $(DECOMP)

.PHONY: all clean
