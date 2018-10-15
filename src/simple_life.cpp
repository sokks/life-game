#include "life_scene.h"

const std::string default_output_folder = "rounds/";
std::string output_folder = default_output_folder;

void set_output_folder(std::string new_folder) {
    output_folder = new_folder;
}

int live_and_die(LifeScene& ls, int max_iters) {
    std::string folder = output_folder;
    int n_of_rounds = 0;
    std::string round_file = std::string("00000.txt");
    // ls.WriteTo(folder + round_file);

    for (int i = 0; i <= max_iters; i++) {
        std::string i_str = std::to_string(i);
        round_file = std::string(5-i_str.length(), '0') + i_str + ".txt";
        ls.WriteTo(folder + round_file);
        int status = ls.Round(-1, -1);

        if (status == 1) {
            n_of_rounds = i;
            break;
        }
        usleep(1);
    }
    return -1;
}

int *read_data(int N, int M, std::string filename) {
    int *data = new int[M*N];
    std::ifstream fin(filename);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            fin >> data[i * M + j];
        }
    }
    fin.close();
    return data;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "usage: live <N> <M> <MAX_ITERS> [<test (0/1)> <test.input> <test.output folder>]\n";
        return 0;
    }
    int N = std::atoi(argv[1]);
    int M = std::atoi(argv[2]);
    int max_iters = std::atoi(argv[3]);

    LifeScene ls;

    if ((argc > 4) && (argc < 6)) {
        std::cout << "usage: live <N> <M> <MAX_ITERS> [<test.input> <test.output folder>]\n";
        return 0;
    } else if (argc > 4) {
        // test
        std::string input_file = std::string(argv[4]);
        std::string out_folder = std::string(argv[5]);
        int *cells = read_data(N, M, input_file);
        ls = LifeScene(M, N, cells);
        set_output_folder(out_folder);
    } else {
        ls = LifeScene(M, N);
    }

    int steps = live_and_die(ls, max_iters);
    std::cout << "all dead in " << steps << " steps\n";
    return 0;
}

