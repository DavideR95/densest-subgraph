#include <vector>
#include <utility>
#include <list>
#include <numeric>

using node_t = uint32_t;

struct bucket_t {
    std::list<node_t> neighs;
    size_t degree;

    int id() { return degree; }

    bool operator < (const bucket_t& cmp) const {
        return this->degree < cmp.degree;
    }

    bucket_t(std::list<node_t> n, size_t deg) : neighs(n), degree(deg) {}
};

class Graph {
    private:
        size_t N; // Number of vertices
        unsigned int b; // Edge duplication parameter
        double eta;
        // Will represent d+(u) sorted by the multiplicity of the (outward) edges
        std::vector<std::vector<std::pair<node_t, unsigned int>>> outward_neighs; 
        // d+(v) = adjlist[v].size()
        // Vector of doubly linked lists for N-(u)
        std::vector<std::list<bucket_t>> inward_neighs;

        size_t edge_count;

        
        // Direction: u -> v
        void add_edge(node_t u, node_t v) {
            unsigned int d_plus_x = INT_MAX;
            node_t x = 0;
            for(auto neigh : outward_neighs[u]) {
                node_t w = neigh.first;
                if(d_plus(w) < d_plus_x) { // argmin (line 1 of Alg. 1)
                    x = w;
                    d_plus_x = d_plus(x);
                }
            }
            

            auto d_plus_u = d_plus(u);
            if((double) d_plus_u + 1 > (double) (1 + eta/b) * d_plus_x + 0) { // Theta = 0
                remove_edge(u, x); // Flip the edge ux to xu
                add_edge(x, u);
            }
            // else {
            append_neigh(u, v); // Effectively edit the adjlist of u with uv
            for(auto w : outward_neighs[u]) {
                // Update d+(u) in the linked lists of N-(w);
                update_n_minus(w.first, u, d_plus_u); // alla sua linked list (di w)
            }
            // }

            edge_count++;

        }

    public: 
        Graph(size_t N_, double eta_, unsigned int b_) : 
            N(N_), b(b_), eta(eta_), outward_neighs(N_), inward_neighs(N_),
            edge_count(0)
        { }

        size_t size() { return N; }
        
        inline size_t d_plus(node_t v) { // Compute d+(v) 
            return std::accumulate(outward_neighs[v].begin(), 
                outward_neighs[v].end(), 0, 
                [](unsigned int a,std::pair<uint32_t, unsigned int> b) {
                    return a + b.second; 
                });
                // d+(u) = sum of all the (duplicated) outward edges of u
        }

        // Adds v to the adjList of u
        void append_neigh(node_t u, node_t v) {
            for(auto& neigh : outward_neighs[u]) {
                if(neigh.first == v) {
                    (neigh.second)++;
                    return;
                }
            }
            // If we got here, v was not in the adjacency list of v
            // Push 1 copy of the edge (u, v)
            outward_neighs[u].push_back(std::make_pair(v, 1));
        }

        void add_undirected_edge(node_t u, node_t v) {
            std::cout << "Hello devo aggiungere " << u << ", " << v << std::endl;
            for(auto i=0;i<b;i++) {
                add_edge(u, v); // Duplicate the edge b times
            }
        }

        void remove_edge(node_t u, node_t v) {
            // First remove it from the outward neighbors of u
            auto d_plus_u = d_plus(u);
            for(auto it=outward_neighs[u].begin();it<outward_neighs[u].end();it++) {
                auto neigh = *it;
                if(neigh.first == v) {
                    if(neigh.second == 1) // There will bo no copies left
                        outward_neighs[u].erase(it); // Remove the pair entirely
                    else 
                        (neigh.second)--; // Decrease the count by one

                    break;
                }
            } 
            // Then, remove it from the inward neighbors of v
            bool inserted = false;
            for(auto& bkt : inward_neighs[v]) {
                if(bkt.id() == d_plus_u) { // Old d+(u)
                    bkt.neighs.remove(u);
                }
                else if(bkt.id() == d_plus_u - 1) {
                    bkt.neighs.push_back(u);
                    inserted = true;
                }
            }

            if(!inserted && d_plus_u - 1 > 0) {
                std::list<node_t> l;
                l.push_back(u);
                bucket_t new_bucket(l, d_plus_u-1);
                inward_neighs[v].emplace_back(new_bucket);
                assert(!inward_neighs[v].back().neighs.empty());
                //inward_neighs[who].front().neighs.push_back(neigh);
                inward_neighs[v].sort(); // VERY inefficient
            }

            edge_count--;
        }        

        void update_n_minus(node_t who, node_t neigh, size_t d_plus_u) {
            bool inserted = false;
            for(bucket_t& bkt : inward_neighs[who]) {
                if(bkt.id() == d_plus_u) {
                    bkt.neighs.remove(neigh);
                }
                else if(bkt.degree == d_plus_u+1) { // Move u in the d(u)+1 bucket
                    bkt.neighs.push_back(neigh);
                    inserted = true;
                    break;
                }
            }

            if(!inserted) { // IF there was no bucket with j = d(u)+1, create it
                // Append a new bucket at the end of the list of w (neigh)
                // Then sort again the buckets
                // Workaround for the debugger I guess
                std::list<node_t> l;
                l.push_back(neigh);
                bucket_t new_bucket(l, d_plus_u+1);
                inward_neighs[who].push_back(new_bucket);
                assert(!inward_neighs[who].back().neighs.empty());
                //inward_neighs[who].front().neighs.push_back(neigh);
                inward_neighs[who].sort(); // VERY inefficient
            }

        }

        void print_all() {
            std::cout << "------- N+ -------" << std::endl;
            for(auto i=0;i<outward_neighs.size();i++) {
                std::cout << "\t[" << i << "]: ";
                for(auto& neigh : outward_neighs[i]) {
                    std::cout << "(" << neigh.first << ", b: " << neigh.second << ") ";
                }
                std::cout << std::endl;
            }

            std::cout << "------- N- -------" << std::endl;

            for(auto i=0;i<inward_neighs.size();i++) {
                std::cout << "\t[" << i << "]: ";
                for(auto& bkt : inward_neighs[i]) {
                    std::cout << "(deg: " << bkt.id() << ": ";
                    for(auto& n : bkt.neighs) std::cout << n << ",";
                    std::cout << " ), ";
                }
                std::cout << std::endl;
            }

            std::cout << "Edge count: " << edge_count << std::endl;

        }


};