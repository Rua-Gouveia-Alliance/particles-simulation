#include <vector>
#include <unordered_map>

class UnionFind {
private:
    std::vector<long long> parent;
    std::vector<long long> rank;

public:
    UnionFind(long long size) : parent(size), rank(size, 1) {
        for (long long i = 0; i < size; ++i) {
            parent[i] = i;
        }
    }

    long long find(long long x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]); 
        }
        return parent[x];
    }

    void unite(long long x, long long y) {
        long long rootX = find(x);
        long long rootY = find(y);
        if (rootX != rootY) {
            // Union by rank
            if (rank[rootX] > rank[rootY]) {
                parent[rootY] = rootX;
            } else if (rank[rootX] < rank[rootY]) {
                parent[rootX] = rootY;
            } else {
                parent[rootY] = rootX;
                rank[rootX]++;
            }
        }
    }
};
