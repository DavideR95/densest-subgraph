#include <vector>
#include <utility>
#include <list>
#include <numeric>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <limits.h>
#include <cassert>
#include <unistd.h>
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

struct parent_node_t { 
    size_t node_count; // How many vertices have degree at least T_i
    size_t first_leaf; // Pointer to first leaf

    parent_node_t() : node_count(0), first_leaf(0) {

    }

};

struct leaf_node_t {
    std::unordered_set<node_t> node_list; // For leaves only
    size_t prev, next;
    
    leaf_node_t() : node_list({}), prev(INT_MAX), next(INT_MAX) {}

    bool isEmpty() { return node_list.empty(); }
};

class tree_t {
    private:
        std::vector<parent_node_t> internal_nodes;
        std::vector<leaf_node_t> leaves;
        std::vector<size_t> reverse_index;
        size_t depth;
        size_t maximum_leaf;
    public:

        tree_t(size_t sz) : internal_nodes(sz-1), reverse_index(sz, 0), maximum_leaf(0) {
            auto number_of_nodes = sz;
            sz--;
            sz |= sz >> 1;
            sz |= sz >> 2;
            sz |= sz >> 4;
            sz |= sz >> 8;
            sz |= sz >> 16;
            sz++; // Next power of 2
            internal_nodes.resize(sz-1);
            update_first_leaf(0);
            leaves.resize(sz);
            depth = std::log2(sz);   
            for(node_t v=0;v<number_of_nodes;v++) leaves[0].node_list.insert(v); // All nodes start with outdegree = 0
        }

        void update_first_leaf(size_t start=0) {
            for(size_t i=start;i<internal_nodes.size();i++) {
                internal_nodes[i].first_leaf = (left(i) < i) ? left(i) : i * std::pow(2, depth-std::log2(i)) - 1 - internal_nodes.size();
            }
        }

        void update_node(node_t who, int count=1) {
            //std::cout << "\tAggiungo " << count << " al nodo " << who;
            size_t bucket = reverse_index[who]; // Move who to the +count bucket
            //std::cout << " che prima si trovava in posizione " << bucket << std::endl;
            leaves[bucket].node_list.erase(who);
            if(bucket == 1047) {
                std::cout << "!!!";
            }
            if(bucket + count >= leaves.size()) {
                internal_nodes.resize((internal_nodes.size() + 1) * 2 - 1);
                leaves.resize(leaves.size() * 2);
                depth++;
                update_first_leaf();
            }
            // STO ASSUMENDO CHE BUCKET+COUNT SIA SEMPRE MAGGIORE DI BUCKET, MA NON È VERO
            leaves[bucket+count].node_list.insert(who);
            if(leaves[bucket].isEmpty()) { // If now this leaf becomes empty, skip it in the list
                if(leaves[bucket].prev < INT_MAX) {
                    leaves[leaves[bucket].prev].next = bucket+count;
                }
                if(leaves[bucket].next < INT_MAX) {
                    leaves[leaves[bucket].next].prev = leaves[bucket].prev;
                }
            }
            if(leaves[bucket+count].prev >= leaves.size()) {
                leaves[bucket+count].prev = (leaves[bucket].isEmpty()) ? leaves[bucket].prev : bucket;
            }
            if(leaves[bucket].next >= leaves.size()) {
                leaves[bucket].next = (leaves[bucket].isEmpty()) ? INT_MAX : bucket+count;
            }
            if(bucket+count > maximum_leaf) maximum_leaf = bucket+count;

            

            // Traverse the tree and fix the count
            size_t curr = parent(bucket, true);
            while(curr != 0) { // Ignores the root
                internal_nodes[curr].node_count--;
                curr = parent(curr);
            }

            curr = parent(bucket+count, true);
            while(curr != 0) { // Ignores the root
                internal_nodes[curr].node_count++;
                curr = parent(curr);
            }

            // Counter of the root does not change (ever, because we do +1 and -1)

            reverse_index[who] = bucket+count;
        }

        size_t parent(size_t index, bool isLeaf=false) {
            return (isLeaf) ? (index + internal_nodes.size() - 1) / 2 : (index-1) / 2;
        }

        size_t left(size_t index, bool isLeaf=false) {
            if(2 * index + 1 >= internal_nodes.size()) { // It means we are the parent of a leaf
                return (isLeaf) ? 0 : index * 2 + 1 - internal_nodes.size();
            }
            return (isLeaf) ? 0 : index * 2 + 1;
        }

        size_t right(size_t index, bool isLeaf=false) {
            if(2 * index + 2 >= internal_nodes.size()) { // It means we are the parent of a leaf
                return (isLeaf) ? 0 : index * 2 + 2 - internal_nodes.size();
            }
            return (isLeaf) ? 0 : index * 2 + 2;
        }

        std::pair<size_t, size_t> tree_leaves(size_t tree_node) {
            int layer = std::log2(tree_node+1);

            size_t leaves_start = (tree_node + 1) * (1 << (depth-layer)) - 1; // n * (2^(depth-layer(n))) - 1
            size_t leaves_end = (tree_node + 2) * (1 << (depth-layer)) - 2; // (n+1) * (2^(depth-layer(n))) - 2

            return std::make_pair(leaves_start - internal_nodes.size(), leaves_end - internal_nodes.size());
        }

        size_t get_ti_size(size_t degree) {
            size_t count = 0;
            // for(auto i=degree;i<leaves.size();i++) {
            //     count += leaves[i].node_list.size();
            // }
            // Ora che so chi è la first leaf di ogni nodo
            // Devo cambiare questo algoritmo in modo da risalire l'albero
            // E poi scansionare le foglie a partire dalla sua prima foglia
            // E stampare i nodi che contengono le foglie
            // Parto dalla foglia degree, prendo il padre, e risalgo
            // Stavolta è più facile, l'unico check da fare è se la foglia
            // è una foglia "destra" o no e ricordarselo
            // Poi la foglia estremo dx è sempre inclusa, quindi facile anche lì
            // Perché la foglia estremo dx è leaves.size()

            //auto savage = count;
            //count = 0;

            std::vector<size_t> result;


            // Align start and end to the actualy nodes in the heap
            auto start = degree;
            auto end = (maximum_leaf % 2 == 0) ? std::min(maximum_leaf + 1, leaves.size() - 1) : maximum_leaf;
            size_t a = start, b = end; // start, end will be modified in the loops
            if(a == b) return leaves[a].node_list.size();

            if(a % 2 != 0) {
                count += leaves[a].node_list.size();
                a++;
                start++;
            }

            size_t curr;

            while(start <= end) {
                std::pair<size_t, size_t> my_leaves;
                curr = start;
                bool isLeaf = true;
                while(true) {
                    if(curr == 0) break;
                    my_leaves = tree_leaves(parent(curr, isLeaf));
                    if(my_leaves.first < a || my_leaves.second > b)
                        break;

                    curr = parent(curr, isLeaf);
                    isLeaf = false; 
                    
                }
                my_leaves = tree_leaves(curr);
                result.push_back(curr);
                if(start < my_leaves.first) {
                    end = std::min(end, my_leaves.first);
                }
                else if(start == my_leaves.first && end >= my_leaves.second) { start = my_leaves.second + 1; }
                if(end > my_leaves.second) { start = std::max(start, my_leaves.second); }
                else if(end == my_leaves.second && start <= my_leaves.first) { end = my_leaves.first + 1; }
            }


            for(auto& index : result) {
                count += internal_nodes[index].node_count;
            }

            //assert(savage == count);
            size_t i = degree;
            while(i<leaves.size()) {
                for(auto& elem : leaves[i].node_list) {
                    std::cout << elem << " ";
                }
                i = leaves[i].next;
            }
            std::cout << std::endl;

            return count;
        }
};

class Graph {
    private:
        size_t N; // Number of vertices
        unsigned int b; // Edge duplication parameter
        double eta;
        // Will represent d+(u) sorted by the multiplicity of the (outward) edges
        std::vector<std::vector<std::pair<node_t, size_t>>> outward_neighs;
        std::vector<size_t> d_pluses;
        // d+(v) = adjlist[v].size()
        // Vector of doubly linked lists for N-(u)
        std::vector<std::list<bucket_t>> inward_neighs;
        size_t maximal_outdegree = 0;
        node_t maximal_node = 0;
        tree_t vertex_tree;
        size_t edge_count;

        bool is_isolated(node_t v) {
            if(outward_neighs[v].empty()) { // If v has at least one outward neighbor
                if(inward_neighs[v].empty()) {
                    return true;
                }
                else {
                    for(auto& bkt : inward_neighs[v]) { // or it has at least one incoming neighbor
                        if(!bkt.neighs.empty())
                            return false;
                    }
                }
            }

            return false;
        }

        // Checks if the node needle is in N-(haystack)
        bool find_v_in_incoming_neighs(node_t haystack, node_t needle) {
            auto d_plus_needle = d_plus(needle);
            for(auto& bkt : inward_neighs[haystack]) { // Traverse the list of neighs
                if(bkt.id() == d_plus_needle) {
                    for(auto& neigh : bkt.neighs) {
                        if(neigh == needle)
                            return true;
                    }
                    return false;
                }
            }
            return false;
        }

        // adds count to the number of incoming edges to who from v
        // N-(who) = N-(who) \cup {v}
        void update_incoming_neighs(node_t who, node_t v, unsigned int count) {
            bool inserted = false;
            // Assumption: d_plus(v) has already been updated
            auto d_plus_v = d_plus(v);

            unsigned short checks = 0;

            for(auto& bkt : inward_neighs[who]) {
                if(bkt.id() == d_plus_v - count) { // Old d+(u)
                    bkt.neighs.remove(v);
                    checks++;
                }
                else if(bkt.id() == d_plus_v) {
                    bkt.neighs.push_back(v);
                    inserted = true;
                    checks++;
                }
                if(checks >= 2) return;
            }

            if(!inserted) {
                std::list<node_t> l;
                l.push_back(v);
                bucket_t new_bucket(l, d_plus_v);
                auto it = inward_neighs[who].begin();
                while(it != inward_neighs[who].end() && (*it).id() < d_plus_v) {
                    it++;
                }
                inward_neighs[who].insert(it, new_bucket);
                //inward_neighs[who].insert(std::lower_bound(inward_neighs[who].begin(), inward_neighs[who].end(), { {}, d_plus_v }), d_plus_v);

                //inward_neighs[who].push_back(new_bucket);
                //assert(!inward_neighs[who].back().neighs.empty());
                //inward_neighs[who].front().neighs.push_back(neigh);
                //inward_neighs[who].sort(); // VERY inefficient
            }
        }

        // Adds the edge (u, v)
        void append_edge(node_t u, node_t v) {
            vertex_tree.update_node(u); // Adds +1 to the outdegree of 1
            for(auto& neigh : outward_neighs[u]) {
                if(neigh.first == v) {
                    neigh.second++;
                    return;
                }
            }
            outward_neighs[u].push_back(std::make_pair(v, 1));
            
            //outward_neighs[u].push_back(v);
            //std::sort(outward_neighs[u].begin(), outward_neighs[u].end());
        }

        // Removes the edge u->v and updates data structures accordingly
        void remove_edge(node_t u, node_t v) {
            // First, remove v from N+(u)
           // auto old_d_plus_u = d_plus(u);
            for(auto it=outward_neighs[u].begin();it<outward_neighs[u].end();it++) {
                auto& neigh = *it;
                /*if(neigh == v) {
                    outward_neighs[u].erase(it);
                    //if(d_pluses[u] > 0) d_pluses[u]--;
                    break;
                } */
                if(neigh.first == v) {
                    // QUI SE È GIA ZERO CI SONO PROBLEMI
                    if(neigh.second == 1) { // If it is zero now remove the pair entirely
                        outward_neighs[u].erase(it);
                    }
                    else {
                        neigh.second--; // Decrease the degree by one 
                    }
                    break;
                }
            }
            /*bool inserted = false;
            for(auto& bkt : inward_neighs[v]) {
                if(bkt.id() == old_d_plus_u) {
                    bkt.neighs.remove(u);
                }
                if(old_d_plus_u > 1) {
                    if(bkt.id() == old_d_plus_u - 1) {
                        bkt.neighs.push_back(u);
                        inserted = true;
                    }
                }
            }

            if(!inserted && old_d_plus_u > 1) {
                std::list<node_t> l;
                l.push_back(u);
                bucket_t new_bucket(l, old_d_plus_u-1);
                inward_neighs[v].push_back(new_bucket);
                assert(!inward_neighs[v].back().neighs.empty());
                //inward_neighs[who].front().neighs.push_back(neigh);

                inward_neighs[v].sort(); // VERY inefficient
            }*/

            vertex_tree.update_node(u, -1); // Decrease outdegree by one
            edge_count--;

        }

        // insErt(u, v)
        void insert_edge(node_t u, node_t v) {
            auto old_d_plus_u = d_plus(u);
            bool skip_checks = false;
            if(!find_v_in_incoming_neighs(v, u)) { // if u \notin N-(v)
                // append (u,v) but don't update the counts yet
                /*bool ok = false;
                for(auto& neigh : outward_neighs[u]) {
                    if(neigh.first == v) {
                        ok = true;
                        break;
                    }
                }
                if(!ok)
                    outward_neighs[u].push_back(std::make_pair(v, 0)); // 0 not to count towards d+(u)*/
                
                append_edge(u, v);
                skip_checks = true;
            }
            else {
                append_edge(u, v);
            }

            // Pseudocode of Alg. 1
            if(!skip_checks) {
                // x = argmin{d+(w) | w \in N+(u)}
                unsigned int d_plus_x = INT_MAX;
                node_t x = INT32_MAX;
                for(auto neigh : outward_neighs[u]) {
                    node_t w = neigh.first;
                    if(d_plus(w) < d_plus_x) { // argmin (line 1 of Alg. 1)
                        x = w;
                        d_plus_x = d_plus(x);
                    }
                }
                
                if(old_d_plus_u + 1 > (((double) 1 + (double) (eta/b)) * d_plus_x + 2*0)) { // Theta = 0
                    remove_edge(u, x); // Flip the edge ux to xu
                    insert_edge(x, u);
                }
                else {
                    // Update d+(u)
                    d_pluses[u]++;
                    if(d_pluses[u] > maximal_outdegree) {
                        maximal_outdegree = d_pluses[u];
                        maximal_node = u;
                    }

                    for(auto& neigh : outward_neighs[u]) {
                        update_incoming_neighs(neigh.first, u, 1); // Add 1 to d+(u) in every N-(neigh)
                    }
                }
            }
            else {

                d_pluses[u]++;
                if(d_pluses[u] > maximal_outdegree) {
                    maximal_outdegree = d_pluses[u];
                    maximal_node = u;
                }

                for(auto& neigh : outward_neighs[u]) {
                    update_incoming_neighs(neigh.first, u, 1); // Add 1 to d+(u) in every N-(neigh)
                }

            }

            edge_count++;
        }

    public: 
        std::vector<std::vector<node_t>> buckets;
        Graph(size_t N_, double eta_, unsigned int b_) : 
                N(N_), b(b_), eta(eta_), outward_neighs(N_),
                d_pluses(N_), inward_neighs(N_), vertex_tree(N_), edge_count(0), 
                buckets(std::log2(N_)  / .15 + 1)
        { }

        ~Graph() {
            size_t length = 0;
            for(auto& list : outward_neighs) {
                length += list.size();
            }
            std::cout << "Lungh. media: " << length / outward_neighs.size() << std::endl;
        }

        const size_t size() const { return N; }

        void update_densest_subgraph(double epsilon=.5) { // eps' = eps/10, gamma = eps', eta > 1280, b = O(eps'^(-2)*eta*log n)
            size_t k = 0;
            std::cout << "Iterazioni: " << (int) (std::log2(size())/epsilon+1) << std::endl;
            size_t count;
            for(auto i=0;i<(int) (std::log2(size())/epsilon+1);i++) {
                double v_i = (double) maximal_outdegree / std::pow((1. + (double)eta/b), i);
                // for(auto v=0;v<size();v++) {
                //     if(d_pluses[v] > v_i) {
                //         buckets[i].push_back(v);
                //     }
                // }
                auto old_count = count;
                count = vertex_tree.get_ti_size(v_i);
                // gamma = 0.05
                if(i > 0 && (double) count < (1.05)*(double)old_count) { k = i-1; break; }
            }

            std::cout << "k è: " << k << std::endl;
            double dens = 0;
            for(auto& v : buckets[k]) {
                std::cout << v << " (deg: " << (double)d_pluses[v]/b << ") ";
                if((double)d_pluses[v]/b > dens) dens = (double)d_pluses[v]/b;
            }
            std::cout << std::endl << "Subgraph density: " << dens;
            // std::cout << "(soglia: " << (double) maximal_outdegree / (double) std::pow((1. + (double)eta/b), k) << " )" << std::endl;
            std::cout << std::endl;

            std::cout << "VERA Densità intrabucket: ";
            double edges = 0.;
            for(auto& v : buckets[k]) {
                for(auto& neigh : outward_neighs[v]) {
                    if(std::find(buckets[k].begin(), buckets[k].end(), neigh.first) != buckets[k].end()) {
                        edges += neigh.second;
                    }
                }
            }
            std::cout << ((double) edges / (double) b) / (double) buckets[k].size()  << std::endl;
            std::cout << "Dimensione: " << buckets[k].size() << std::endl;
            for(auto& b : buckets) b.clear(); 

        }

        std::vector<node_t> density_grouping(double epsilon=.25, double gamma = .01) {
            std::vector<std::vector<node_t>> buckets(std::log2(size())  / epsilon + 1);
            auto v_i = delta();
            double i = 0.;
            while(i < (double) std::log2(size())  / epsilon) { // i \in [0, log n / epsilon]
                std::cout << "e alla meno uno log n: " << std::log2(size())  / epsilon << std::endl;
                double density = (double) v_i.second / std::pow((1. + (double) eta/b), i); // Delta(Gb) * (1+eta/b)^-i
                std::cout << "La densità è: " << (double) v_i.second / b << std::endl;
                std::cout << "La densità diviso " << i << " è " << (double) density << " e i vertici sono: " << std::endl;
                std::cout << "epsilon/log n fa: " << epsilon << "/" << std::log2(size()) << " = " << (double) std::log2(size())  / epsilon << std::endl;

                for(node_t v=0;v<size();v++) {
                    if(d_pluses[v] >= density) {
                        //std::cout << v << ", ";
                        buckets[i].push_back(v);
                    }
                }
                // Calcolare densita intrabucket
                std::cout << std::endl;
                std::cout << "VERA Densità intrabucket: ";
                double edges = 0.;
                for(auto& v : buckets[i]) {
                    for(auto& neigh : outward_neighs[v]) {
                        if(std::find(buckets[i].begin(), buckets[i].end(), neigh.first) != buckets[i].end()) {
                            edges++;
                        }
                    }
                }
                std::cout << (double) edges / buckets[i].size() / (double) b << std::endl;
                i++;
                //if(i > 10 && buckets[i].size() < (1.+gamma)*(double)buckets[i-1].size()) break;
            }
            std::cout << "No. of iterations: " << i << std::endl;

            return buckets[i-1];
        }

        // Maximum OUT DEGREE 
        std::pair<node_t, size_t> delta() {
            // node_t who;
            // size_t max = 0;
            // for(auto i=0;i<N;i++) {
            //     if(d_pluses[i] > max) {
            //         who = i;
            //         max = d_pluses[i];
            //     }
            // }
            //assert(who == maximal_node && max == maximal_outdegree);
            return std::make_pair(maximal_node, maximal_outdegree);
        }

        std::vector<std::vector<node_t>> rasterize() {
            // Metodo che trasforma il multigrafo in grafo normale con un solo arco orientato in base al max tra (a, b) e (b, a)
            std::vector<std::vector<std::pair<node_t, size_t>>> edges(N);
            std::vector<std::vector<node_t>> rasterized(N);
            for(auto i=0;i<N;i++) { // Loop through the vertices
                for(auto& neigh : outward_neighs[i]) {
                    // if(edges[i].empty() || edges[i].back().first != neigh) {
                    //     edges[i].push_back(std::make_pair(neigh, 1));
                    // }
                    // else {
                    //     edges[i].back().second++;
                    // }
                }
                // outward_neighs[i].clear();
            }

            for(auto i=0;i<N;i++) {
                auto& list = edges[i];
                //std::cout << "[" << i << "] ";
                for(auto it=list.begin();it<list.end();it++) {
                    auto& pair = *it;
                    auto v = std::find_if(edges[pair.first].begin(), edges[pair.first].end(), [&](const std::pair<node_t, size_t>& element) { return element.first == i; } );
                    if(v != edges[pair.first].end()) {
                        //std::cout << "Trovato " << i << " nella lista di " << pair.first << " con grado = " << v->second << "\t";
                        if(pair.second >= v->second)
                            rasterized[i].push_back(pair.first); // Choose (a->b)
                        else 
                            rasterized[pair.first].push_back(i); // Choose (b->a)

                        // list.erase(it);
                        // edges[pair.first].erase(std::remove_if(edges[pair.first].begin(), edges[pair.first].end(), [&](const std::pair<node_t, size_t>& element) { return element.first == i; }), edges[pair.first].end());
                        
                    }
                    else {
                        //assert(false);
                        rasterized[i].push_back(pair.first); // Just a->b
                        // list.erase(it);
                        // edges[pair.first].erase(std::remove_if(edges[pair.first].begin(), edges[pair.first].end(), [&](const std::pair<node_t, size_t>& element) { return element.first == i; }), edges[pair.first].end());
                    }
                }
            }

            for(auto& list : rasterized) {
                std::sort(list.begin(), list.end());
                list.erase(std::unique(list.begin(), list.end()), list.end());
            }

            return rasterized;
        }

        // INSORT(u, v)
        void add_undirected_edge(node_t u, node_t v) {

            if(is_isolated(u) && is_isolated(v)) {
                // outward_neighs[u].push_back(std::make_pair(v, std::ceil((double)b/2)));
                // outward_neighs[v].push_back(std::make_pair(u, b/2));
                for(auto i=0;i<std::ceil((double)b/2);i++)
                    append_edge(u, v);
                for(auto i=0;i<b/2;i++)
                    append_edge(v, u);

                d_pluses[u] += std::ceil((double)b/2);
                if(d_pluses[u] > maximal_outdegree) {
                    maximal_outdegree = d_pluses[u];
                    maximal_node = u;
                }

                d_pluses[v] += b/2;
                if(d_pluses[v] > maximal_outdegree) {
                    maximal_outdegree = d_pluses[v];
                    maximal_node = v;
                }

                update_incoming_neighs(u, v, b/2); // add v to N-(u)
                update_incoming_neighs(v, u, std::ceil((double)b/2)); // add u to N-(v)
                edge_count += b;
            }
            else {
                for(auto i=0;i<b;i++) {
                    if(d_plus(u) <= d_plus(v)) {
                        insert_edge(u, v);
                    }
                    else {
                        insert_edge(v, u);
                    }
                }
            }
        }

        inline size_t d_plus(node_t v) { // Compute d+(v) 
            return d_pluses[v];
            // return std::accumulate(outward_neighs[v].begin(), 
            //     outward_neighs[v].end(), 0, 
            //     [](unsigned int a,std::pair<node_t, unsigned int> b) {
            //         return a + b.second; 
            //     });
                // d+(u) = sum of all the (duplicated) outward edges of u
        }

        void print_all() {
            std::cout << "------- N+ -------" << std::endl;
            for(auto i=0;i<outward_neighs.size();i++) {
                std::cout << "\t[" << i << " d+ = " << d_plus(i) << ", length=" << outward_neighs[i].size() << "]: ";
                for(auto& neigh : outward_neighs[i]) {
                    //std::cout << "(" << neigh.first << ", b: " << neigh.second << ") ";
                    std::cout << neigh.first << " ";
                }
                std::cout << std::endl;
            }

            /* std::cout << "------- N- -------" << std::endl;

            for(auto i=0;i<inward_neighs.size();i++) {
                std::cout << "\t[" << i << "]: ";
                for(auto& bkt : inward_neighs[i]) {
                    std::cout << "(deg: " << bkt.id() << " (size: " << bkt.neighs.size() << "): ";
                    for(auto& n : bkt.neighs) std::cout << n << ",";
                    std::cout << " ), ";
                }
                std::cout << std::endl;
            }
            */

            std::cout << "Edge count: " << edge_count << std::endl;

        }

};