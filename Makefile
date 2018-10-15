N=100
MAX_ITER=100

MPI=1
N_PROCS=2

COMPILER=g++
EXECUTER=

test: bin/life_game bin/life_game_mpi prepare test/compare_all.py
	mkdir -p test/expected
	mkdir -p test/got
	bin/life_game $(N) $(N) $(MAX_ITER) test/test_100_100.input test/expected/
	mpirun -n $(N_PROCS) bin/life_game_mpi $(N) $(N) $(MAX_ITER) test/test_100_100.input test/got/
	python3 test/compare_all.py $(N) $(N) test/expected/ test/got/

bench: bin/life_game bin/life_game_mpi prepare
	# bin/life_game 1000 1000 1000
	mpirun -n 1 bin/life_game_mpi 1000 1000 100
	mpirun -n 2 bin/life_game_mpi 1000 1000 100
	mpirun -n 4 bin/life_game_mpi 1000 1000 100

run: bin/life_game prepare 
	bin/life_game $(N) $(N) $(MAX_ITER)
	python3 visualize.py $(N) $(N)

run_mpi: bin/life_game_mpi prepare
	mpirun -n $(N_PROCS) bin/life_game_mpi $(N) $(N) $(MAX_ITER)
	python3 visualize.py $(N) $(N)

prepare:
	mkdir -p rounds
	rm -f rounds/*

bin/life_game: src/simple_life.cpp build/life_scene.o
	mkdir -p bin
	g++ -o $@ src/simple_life.cpp build/life_scene.o

bin/life_game_mpi: src/life_mpi.cpp build/life_scene_mpi.o
	mkdir -p bin
	mpic++ -std=c++11 -o $@ src/life_mpi.cpp build/life_scene_mpi.o

build/life_scene.o: src/life_scene.h src/life_scene.cpp Makefile
	mkdir -p build
	$(COMPILER) -o $@ -c src/life_scene.cpp

build/life_scene_mpi.o: src/life_scene.h src/life_scene.cpp Makefile
	mkdir -p build
	mpic++ -o $@ -c src/life_scene.cpp

visualize:
	python3 visualize.py $(N) $(N) 

clean:
	rm -rf rounds
	rm -rf bin
	rm -rf build
	rm -rf test/expected
	rm -rf test/got

mpi_hello_world:
	mpicc -o mpi_hello_world mpi_hw.c
	mpirun -n 2 mpi_hello_world

