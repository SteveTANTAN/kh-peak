#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <chrono>
#include<time.h>
#include<ctime>


using namespace std;
int vertexnum = 0, edgenum = 0;

typedef int nI;  //节点的编号类型int
typedef int eI;  //边的编号类型int
typedef map<nI, nI> MII;
typedef vector<vector<int>> VVI;

struct Vertex {
    vector<int> neighbor;   //1-hop邻点集
    vector<int> hneighbor;  //h-hop邻点集
    vector<int> shell;      //shell内邻点(opt)
    int hdegree = 0;
    int sdegree = 0;        //k-shell层内度数(opt)
    bool deleted = false;
    bool visited = false;
    int hcoreness = 0;
    int hpeakness = 0;
    int kbound = 0;
};

struct Edge {
    nI id0 = 0;  //边中小id端点
    nI id1 = 0;  //边中大id端点
    bool deleted = false;
    int support = 0;
    int trussness = 0;
};

struct Graph {
    int n = 0;
    int m = 0;
    vector<Vertex> node;
    vector<Edge> edge;
};

void deleteNode(Graph& g, int v) {
    g.node[v].deleted = true;
    g.n--;
    if (g.node[v].neighbor.empty())
        return;
    for (int i = 0; i < g.node[v].neighbor.size(); i++) {
        int u = g.node[v].neighbor[i];
        for (int j = 0; j < g.node[u].neighbor.size(); j++) {
            if (g.node[u].neighbor[j] == v) {
                g.node[u].neighbor[j] = g.node[u].neighbor.back();
                break;
            }
        }
        g.node[u].neighbor.pop_back();
    }
}

void get_hneighbor(Graph& g, int v, int h) {
    g.node[v].hneighbor.clear();
    g.node[v].hdegree = 0;

    vector<bool> hvisited(g.node.size(), false);
    queue<int> Q1;
    hvisited[v] = true;
    for (auto u : g.node[v].neighbor) {
        Q1.push(u); //先载入1-hop邻居
        hvisited[u] = true;
    }

    int d = 1;
    while (d <= h) {
        queue<int> Q2;  //Q2每轮循环初始为空
        while (!Q1.empty()) {
            int u = Q1.front();
            Q1.pop();
            g.node[v].hneighbor.push_back(u);
            g.node[v].hdegree++;
            for (auto w : g.node[u].neighbor) {
                if (!hvisited[w]) {
                    Q2.push(w);
                    hvisited[w] = true;
                }
            }
        }
        d++;
        Q1.swap(Q2);
    }
}

void GraphAToB(Graph g1, Graph& g2) {
    g2.n = g1.n;
    g2.node.resize(g1.node.size());
    for (int v = 1; v < g1.node.size(); ++v) {
        g2.node[v].deleted = g1.node[v].deleted;
        if (!g1.node[v].deleted) {
            g2.node[v].neighbor = move(g1.node[v].neighbor);
        } else {
            g2.node[v].neighbor = {0};
        }
        g2.node[v].hneighbor = {0};
        g2.node[v].hdegree = 0;
        g2.node[v].hcoreness = 0;
        g2.node[v].hpeakness = 0;
    }
}

void get_khcore(Graph& g, int h) {
    int n = g.n;
    VVI bin(g.node.size()); // Resize bin based on the size of the node vector

    // Each round of hpeak decomposition deletes vertices and reorganizes the remaining vertices
    for (int v = 1; v < g.node.size(); ++v) {
        if (g.node[v].deleted) continue;
        get_hneighbor(g, v, h);
        bin[g.node[v].hdegree].push_back(v);
    }

    for (int k = 0; k < bin.size(); ++k) {
        while (!bin[k].empty()) {
            int v = bin[k].back();
            bin[k].pop_back();
            g.node[v].hcoreness = k;
            deleteNode(g, v);

            for (auto u : g.node[v].hneighbor) {
                int pre_hdegree = max(g.node[u].hdegree, k);
                get_hneighbor(g, u, h);
                int temp = max(g.node[u].hdegree, k);
                bin[temp].push_back(u);

                // Remove u from the original bin layer
                for (int i = 0; i < bin[pre_hdegree].size(); i++) {
                    if (bin[pre_hdegree][i] == u) {
                        bin[pre_hdegree][i] = bin[pre_hdegree].back();
                        bin[pre_hdegree].pop_back();
                        break;
                    }
                }
            }
        }
    }
}

void BaselineAlgorithm(Graph& g, int h) {
    while (g.n > 0) {
        Graph g_copy;
        GraphAToB(g, g_copy);
        auto StartTime = chrono::steady_clock::now();

        get_khcore(g_copy, h);

        auto EndTime = chrono::steady_clock::now();
        if (g.n == g.node.size() - 1) {
            auto duration = chrono::duration_cast<chrono::milliseconds>(EndTime - StartTime).count();
            cout << " Core run time: " << duration << " ms" << endl;
        }

        int maxhcore = 0;
        for (int v = 1; v < g.node.size(); v++) {
            if (g.node[v].deleted) continue;
            if (maxhcore < g_copy.node[v].hcoreness) {
                maxhcore = g_copy.node[v].hcoreness;
            }
        }

        if (maxhcore == 0) {
            for (int u = 1; u < g.node.size(); u++) {
                if (g.node[u].deleted) continue;
                g.node[u].hpeakness = 0;
            }
            break;
        }

        queue<int> Q;
        for (int i = 1; i < g.node.size(); i++) {
            if (g.node[i].deleted) continue;
            if (g_copy.node[i].hcoreness == maxhcore) {
                Q.push(i);
            }
        }

        while (!Q.empty()) {
            int v = Q.front();
            Q.pop();
            g.node[v].hpeakness = maxhcore;
            deleteNode(g, v);
        }
    }
}

int main(int argc, char* argv[]) {
    //for (int h = 1; h <= 1; h++) {
    //cout << " -----   H = " << h << '\n' << endl;
        if (argc < 2) {
            cout << "Input file path missing!" << endl;
            return 1;
        }
        int h = atoi(argv[2]);
        string filePath = argv[1];
        ifstream inputFile(filePath);
        if (!inputFile) {
            cerr << "Input file could not be opened!" << endl;
            return 1;
        }


        int maxVertexIndex = 0;
        inputFile >> vertexnum >> edgenum;
        
        for (int i = 1; i <= edgenum; i++)  {
            int u, v;
            inputFile >> u >> v;
            maxVertexIndex = max(maxVertexIndex, max(u, v));
        }
        inputFile.close();
        Graph g;

        g.n = maxVertexIndex;
        g.m = edgenum;
        g.node.resize(maxVertexIndex + 1);

        inputFile.open(filePath);
        if (!inputFile) {
            cerr << "Input file could not be opened!" << endl;
            exit(EXIT_FAILURE);
        }

        inputFile >> vertexnum >> edgenum;
        for (int i = 1; i <= edgenum; i++) {
            int u, v;
            inputFile >> u >> v;
            if (u == v) continue;
            g.node[u].neighbor.push_back(v);
            g.node[v].neighbor.push_back(u);
        }
        inputFile.close();
        vertexnum = maxVertexIndex;
        cout << " Graph g loaded successfully!" << '\n' << endl;
        cout << " -----   H = " << h << " -----" << endl;

        struct timeval start, end;
        gettimeofday(&start, NULL);

        BaselineAlgorithm(g, h);

        gettimeofday(&end, NULL);
        double timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec) / 1000000.0;
        cout << "Baseline run time: " << timeuse << " s" << endl;

        /* Output the Result */
        ofstream outputFile("sample.txt");

        
        if (!outputFile) {
            cerr << "Output file could not be opened!" << endl;
            exit(EXIT_FAILURE);
        }
        for (int i = 1; i <= vertexnum; i++) {
            outputFile << i << ' ' << g.node[i].hpeakness << '\n';
        }
        outputFile.close();
    // }
    return 0;
}
