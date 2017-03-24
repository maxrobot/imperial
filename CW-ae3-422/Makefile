# Declaration of variables
CC = mpicxx
CC_FLAGS = -O2 -std=c++11
LD_FLAGS = -lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -llapack -lblas
DEBUG_FLAGS = -g

COTEXT  = "\033[1;31m Compiling\033[0m\033[1m $(<F)\033[0m\n"
LITEXT  = "\033[1;31m Building \033[0m\033[1m $@\033[0m\n"

# File names
EXEC    = $(BINDIR)hpc_solver
SRCDIR = ./src/
OBJDIR = ./.Obj/
BINDIR = ./bin/
MKDIRS = $(OBJDIR) $(BINDIR)
CPP_FILES = $(wildcard $(SRCDIR)*.cpp)
OBJ_FILES = $(addprefix $(OBJDIR), $(notdir $(CPP_FILES:.cpp=.o)))

# Main target
compile: author compile_info $(MKDIRS) $(OBJ_FILES) 
	@echo
	@printf $(LITEXT) 
	@$(CC) $(OBJ_FILES) -o $(EXEC) $(LD_FLAGS)

# To obtain object files
$(OBJDIR)%.o: $(SRCDIR)%.cpp 
	@printf $(COTEXT) 
	@$(CC) -c $(CC_FLAGS) $< -o  $@

all: task1 task2 task3 task4 task5

# To remove generated files
clean: 
	@printf "\033[1;31m Cleaning object, executable and archive files\033[0m\n"
	@rm -f $(EXEC) $(OBJ_FILES) solver.tgz
	@./clear

# To create a tar archive of the project
tar: clean
	@printf "\033[1;31m Creating tar archive of the code \033[0m\n"
	@mkdir -p solver_01220997
	@cp -r Makefile src solver_01220997
	@tar czf solver_01220997.tgz --exclude .git --exclude .Obj --exclude output solver_01220997
	@rm -rf solver_01220997

$(MKDIRS): 
	@mkdir -p $@

compile_info:
	@printf "\033[1;31m Compiler used \033[0m \033[1m    [$(CC)]\033[0m\n"
	@printf "\033[1;31m Compiling options\033[0m \033[1m [$(CC_FLAGS)]\033[0m\n" 
	@printf "\033[1;31m Linking options \033[0m \033[1m  [$(LD_FLAGS)]\033[0m\n"
	@echo

author:
	@echo
	@echo    "       *************************************************************"
	@echo    "       *                                                           *"
	@printf    "       *  \033[0m \033[1mAE3-422\033[0m High-Performance Computing                      *\n"
	@printf    "       *  \033[0m \033[1mCoursework Assignment\033[0m                                   *\n"
	@echo    "       *                                                           *"
	@printf    "       *  \033[0m \033[1mAuthor\033[0m: FM Grabner                                      *\n"
	@printf    "       *  \033[0m \033[1mCID\033[0m:    01220997                                        *\n"
	@echo    "       *                                                           *"
	@echo    "       *************************************************************"
	@printf    "       *                  \033[0m \033[1mImperial College London\033[0m                 *\n"
	@printf "       *************************************************************\033[0m\n"
	@echo

# To print the help message
help: author
	@echo
	@printf "\033[1;31m Make options of solver code\033[0m\n"
	@echo
	@printf "\033[1;31m Compiler used \033[0m \033[1m    [$(CC)]\033[0m\n"
	@printf "\033[1;31m Compiling options\033[0m \033[1m [$(CC_FLAGS)]\033[0m\n" 
	@printf "\033[1;31m Linking options \033[0m \033[1m  [$(LD_FLAGS)]\033[0m\n"
	@echo	
	@printf "\033[1;31m Provided Rules: \033[0m\n"
	@printf "\033[1;31m  all          =>\033[0m\033[1m run all tasks for this assignment\033[0m\n"
	@printf "\033[1;31m  clean        =>\033[0m\033[1m cleaning object, executable and archive files\033[0m\n"
	@printf "\033[1;31m  compile      =>\033[0m\033[1m compiles solver\033[0m\n"
	@printf "\033[1;31m  help         =>\033[0m\033[1m printing this help message\033[0m\n"
	@printf "\033[1;31m  tar          =>\033[0m\033[1m creating a tar archive of the project\033[0m\n"
	@printf "\033[1;31m  task1        =>\033[0m\033[1m runs task 1\033[0m\n"
	@printf "\033[1;31m  task2        =>\033[0m\033[1m runs task 2\033[0m\n"
	@printf "\033[1;31m  task3        =>\033[0m\033[1m runs task 3\033[0m\n"
	@printf "\033[1;31m  task4        =>\033[0m\033[1m runs task 4\033[0m\n"
	@printf "\033[1;31m  task5        =>\033[0m\033[1m runs task 5\033[0m\n"
	@echo

task1:
	@mpiexec -np 1 bin/hpc_solver param_file/task1.txt && python output/scripts/task1_sol.py

task2:
	@mpiexec -np 1 bin/hpc_solver param_file/task2.txt && python output/scripts/task2_sol.py

task3:
	@mpiexec -np 1 bin/hpc_solver param_file/task3.txt && python output/scripts/task3_sol.py

task4:
	@mpiexec -np 2 bin/hpc_solver param_file/task4.txt && python output/scripts/task4_sol.py

task5:
	@mpiexec -np 4 bin/hpc_solver param_file/task5.txt && python output/scripts/task5_sol.py