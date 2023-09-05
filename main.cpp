#include <iostream>
#include <fstream>

#include "graph.hpp"

using node_t = uint32_t;

int main(int argc, char* argv[]) {


    if(argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <eta> <b> <edgeList file> "
                  << std::endl;
        exit(-1);
    }

    std::ifstream file(argv[3]);

    double eta = std::stod(argv[1]);
    unsigned int b = (unsigned int) std::stoul(argv[2]);
    size_t N;

    file >> N;

    std::cout << "Params: eta = " << eta << ", b = " << b << ", N = " << N << std::endl;

    Graph g(N, eta, b);

    while(!file.eof()) {
        node_t u, v;

        file >> u >> v;
        std::cout << "Leggo " << u << ", " << v << " ";

        g.add_undirected_edge(u, v);
        std::cout << "Nuova max dens: " << (double) g.delta().second / b << " (nodo: " << g.delta().first << ")" << std::endl;
        // Forse si può riutilizzare il codice di costruzione per lo heap......
        // E poi però cosa ci sta nelle foglie? Pensiamoci oggi....
    }
    
    g.update_densest_subgraph();

    std::cout << "Grafo risultante: " << std::endl;

    //g.print_all();

    std::cout << "Nodo max: " << g.delta().first << " density: " << (double) g.delta().second/b << std::endl;

    // auto densest_subgraph = g.density_grouping(.15);

    // std::ofstream outfile("densest.txt");

    // for(auto& v : densest_subgraph) {
    //     outfile << v << " ";
    // }

    // outfile.close();


    // auto rasterized = g.rasterize();
    // for(auto i = 0;i<rasterized.size();i++) {
    //     std::cout << "\t[" << i << "]: ";
    //     for(auto& neigh : rasterized[i]) {
    //     //std::cout << "(" << neigh.first << ", b: " << neigh.second << ") ";
    //         std::cout << neigh << " ";
    //     }
    //     std::cout << std::endl;
    // }

    return 0;
}