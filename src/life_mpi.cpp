#include <mpi.h>
#include <stdio.h>
#include <sstream>

#include "life_scene.h"

#define DEBUG_PRINT 0

std::string default_output_folder = "rounds/";
std::string output_folder = default_output_folder;

void set_output_folder(std::string new_output_folder) {
    output_folder = new_output_folder;
}

std::string row_to_string(int len, int*row) {
    std::string res;
    for (int i = 0; i < len; i++) {
        res += std::to_string(row[i]);
    }

    return res;
}

int counter = 0;

struct TimeStat {
    double common = 0.0;
    double io = 0.0;
    double exchange = 0.0;
    double round_avg = 0.0;
    double compute = 0.0;

    std::string Format() {
        std::ostringstream stringStream;
        stringStream << "common: " << common << std::endl;
        stringStream << "io: " << io << std::endl;
        stringStream << "exchange: " << exchange << std::endl;
        stringStream << "compute: " << compute << std::endl;
        stringStream << "round average: " << round_avg << std::endl;

        return stringStream.str();
    }
};

class LifeProc {
    int empty = 1;
    int comm_size, my_rank;

    int M, N;
    int my_N, my_M;
    int prev, next;
    int my_first_row_n, my_last_row_n;

    LifeScene myLifeScene;
    int myStatus = 0;
    int mySteps = 0;

    // for master gather
    int *all_data = nullptr;

public:
    TimeStat ts;

    LifeProc(){}
    LifeProc(int cs, int r, int n, int m) : comm_size(cs), my_rank(r), M(m), N(n) {
        my_N = n / comm_size;
        my_M = m;
        
        if ( (comm_size > 1) && (my_rank > 0) ) {
            prev = my_rank - 1;
            my_first_row_n = 1;
        } else {
            prev = -1;
            my_first_row_n = 0;
        }

        if ( (comm_size > 1) && (my_rank < comm_size - 1) ) {
            next = my_rank + 1;
            my_last_row_n = my_N + my_first_row_n;
        } else {
            next = -1;
            my_last_row_n = my_N + my_first_row_n;
        }

        myLifeScene = LifeScene(my_M, my_N + (prev == -1 ? 0 : 1) + (next == -1 ? 0 : 1));


        // std::cout << "Proc " << my_rank << " inited LifeProc: my_N = " << my_N <<  
        //                                     " next = " << next << " prev = " << prev <<
        //                                     " my_first_row_n = " << my_first_row_n << 
        //                                     " my_last_row_n = " << my_last_row_n << std::endl;
    }

    void ExchangeNeighsRows() {
        double exhange_time = 0.0;
        double tmp_time;

        // pair exchange 1
        if (my_rank % 2 == 0) {

            // send & recv last rows
            if (next != -1) {
                int *my_last_row = myLifeScene.GetRowPtr(my_last_row_n-1);
                int *next_row = myLifeScene.GetRowPtr(my_last_row_n);

                if (DEBUG_PRINT) {
                    std::cout << my_rank << " sending row (last) [ " + row_to_string(M, my_last_row) + "] to " << next << std::endl;
                }

                tmp_time = MPI_Wtime();
                MPI_Status status;
                MPI_Sendrecv(my_last_row, my_M, MPI_INT, next, 0, next_row, my_M, MPI_INT, next, 0, MPI_COMM_WORLD, &status);
                exhange_time += MPI_Wtime() - tmp_time;


                if (DEBUG_PRINT) {
                    std::cout << my_rank << " received row (next) [ " + row_to_string(M, next_row) + "] from " << next << std::endl;
                }
            }
        } else {
            // send & recv first rows
            if (prev != -1) {
                int *my_first_row = myLifeScene.GetRowPtr(my_first_row_n);
                int *prev_row = myLifeScene.GetRowPtr(my_first_row_n-1);

                if (DEBUG_PRINT) {
                    std::cout << my_rank << " sending row (first) [ " + row_to_string(M, my_first_row) + "] to " << prev << std::endl;
                }

                tmp_time = MPI_Wtime();
                MPI_Status status;
                MPI_Sendrecv(my_first_row, my_M, MPI_INT, prev, 0, prev_row, my_M, MPI_INT, prev, 0, MPI_COMM_WORLD, &status);
                exhange_time += MPI_Wtime() - tmp_time;

                if (DEBUG_PRINT) {
                    std::cout << my_rank << " received row (prev) [ " + row_to_string(M, prev_row) + "] from " << prev << std::endl;
                }
            }
        }

        // pair exchange 2
        if (my_rank % 2 == 0) {
            // send & recv first rows
            if (prev != -1) {
                int *my_first_row = myLifeScene.GetRowPtr(my_first_row_n);
                int *prev_row = myLifeScene.GetRowPtr(my_first_row_n-1);
                
                if (DEBUG_PRINT) {
                    std::cout << my_rank << " sending row (first) [ " + row_to_string(M, my_first_row) + "] to " << prev << std::endl;
                }

                tmp_time = MPI_Wtime();
                MPI_Status status;
                MPI_Sendrecv(my_first_row, my_M, MPI_INT, prev, 0, prev_row, my_M, MPI_INT, prev, 0, MPI_COMM_WORLD, &status);
                exhange_time += MPI_Wtime() - tmp_time;

                if (DEBUG_PRINT) {
                    std::cout << my_rank << " received row (prev) [ " + row_to_string(M, prev_row) + "] from " << prev << std::endl;
                }
            }
        } else {
            // send & recv last rows
            if (next != -1) {
                int *my_last_row = myLifeScene.GetRowPtr(my_last_row_n-1);
                int *next_row = myLifeScene.GetRowPtr(my_last_row_n);

                if (DEBUG_PRINT) {
                    std::cout << my_rank << " sending row (last) [ " + row_to_string(M, my_last_row) + "] to " << next << std::endl;
                }

                tmp_time = MPI_Wtime();
                MPI_Status status;
                MPI_Sendrecv(my_last_row, my_M, MPI_INT, next, 0, next_row, my_M, MPI_INT, next, 0, MPI_COMM_WORLD, &status);
                exhange_time += MPI_Wtime() - tmp_time;

                if (DEBUG_PRINT) {
                    std::cout << my_rank << " received row (next) [ " + row_to_string(M, next_row) + "] from " << next << std::endl;
                }
            }
        }

        ts.exchange += exhange_time;
    }

    void Round() {
        double round_time = MPI_Wtime();

        // dump
        // myLifeScene.WriteTo(std::string("test/procs/") + std::to_string(my_rank) + std::to_string(counter) + ".txt");
        // counter++;

        myStatus = myLifeScene.Round(my_first_row_n, my_last_row_n);

        round_time = MPI_Wtime() - round_time;
        ts.compute += round_time;
        mySteps++;
    }

    void ReadFromAndScatter(std::string filename) {
        double io_time = 0.0;
        double exchange_time = 0.0;
        double tmp_time_1, tmp_time_2;

        if (my_rank == 0) {
            tmp_time_1 = MPI_Wtime();
            if (all_data == nullptr) {
                all_data = new int[M*N];
            }

            std::ifstream fin(filename);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M; j++) {
                    fin >> all_data[i * M + j];
                }
            }
            fin.close();
            io_time += MPI_Wtime() - tmp_time_1;
        }
        
        int *my_data = myLifeScene.GetRowPtr(my_first_row_n);

        tmp_time_2 = MPI_Wtime();
        MPI_Scatter(all_data, my_N*my_M, MPI_INT, my_data, my_N*my_M, MPI_INT, 0, MPI_COMM_WORLD);
        exchange_time += MPI_Wtime() - tmp_time_2;

        ts.io += io_time;        
        ts.exchange += exchange_time;
    }

    void GatherAndWriteTo(std::string filename) {
        double io_time = 0.0;
        double exchange_time = 0.0;
        double tmp_time_1, tmp_time_2;

        int *my_data = myLifeScene.GetRowPtr(my_first_row_n);
        
        if (my_rank == 0) {
            if (all_data == nullptr) {
                all_data = new int[M*N];
            }
        }
        
        tmp_time_2 = MPI_Wtime();
        MPI_Gather(my_data, my_N*my_M, MPI_INT, all_data, my_N*my_M, MPI_INT, 0, MPI_COMM_WORLD);
        exchange_time += MPI_Wtime() - tmp_time_2;

        if (my_rank == 0) {
            tmp_time_1 = MPI_Wtime();
            LifeScene ls = LifeScene(M, N, all_data);
            ls.WriteTo(filename);
            io_time += MPI_Wtime() - tmp_time_1;
        }

        ts.io += io_time;
        ts.exchange += exchange_time;
    }

    int ReduceStatus() {
        int all_status;

        double exchange_time = MPI_Wtime();
        MPI_Allreduce(&myStatus, &all_status, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

        exchange_time = MPI_Wtime() - exchange_time;
        ts.exchange += exchange_time;

        return all_status;
    }

    std::string FormatTimeStat() {
        ts.round_avg = ts.compute / mySteps;
        return "Proc: " + std::to_string(my_rank) + "\n" + 
                "Scene size: " + std::to_string(my_N) + " x " + std::to_string(my_M) + "\n" +
                ts.Format();
    }
};

int live_and_die(LifeProc& lp, int max_iters);

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int comm_size;
    int my_rank;
    int name_len;

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (argc < 4) {
        MPI_Finalize();
        std::cout << "usage: live <N> <M> <MAX_ITERS>\n";
        return 0;
    } else if ((argc > 4) && (argc < 6)) {
        MPI_Finalize();
        std::cout << "usage: live <N> <M> <MAX_ITERS> [<test.input> <test output folder>]\n";
        return 0;
    }

    int N = std::atoi(argv[1]);
    int M = std::atoi(argv[2]);
    int max_iters = std::atoi(argv[3]);

    LifeProc lp(comm_size, my_rank, N, M);

    if (argc > 4) {
        std::string input_file = std::string(argv[4]);
        lp.ReadFromAndScatter(input_file);
        std::string out_folder = std::string(argv[5]);
        set_output_folder(out_folder);
    }

    double time_start = MPI_Wtime();
    int steps = live_and_die(lp, max_iters);
    lp.ts.common = MPI_Wtime() - time_start;

    std::cout << lp.FormatTimeStat();

    MPI_Finalize();
    return 0;
}



int live_and_die(LifeProc& lp, int max_iters) {
    int status = 0;

    std::string folder = output_folder;
    std::string round_file = std::string("00000.txt");

    for (int rounds = 0; rounds <= max_iters; rounds++) {
        std::string rounds_str = std::to_string(rounds);
        round_file = std::string(5-rounds_str.length(), '0') + rounds_str + ".txt";
        lp.GatherAndWriteTo(folder + round_file);

        lp.ExchangeNeighsRows();
        lp.Round();

        if (lp.ReduceStatus() != 0) {
            return rounds;
        }
    }

    return max_iters;
}
