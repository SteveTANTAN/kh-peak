#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include<queue>
#include<map>
#include<time.h>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<cstdlib>
#include<ctime>
#include<sys/time.h>
#include<omp.h> 
#include <chrono>

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


void deleteNode(Graph& g, int v)  //删除顶点及其相关边（仅进行在1-hop原图中）
{
	g.node[v].deleted = true;
	g.n--;
	if (g.node[v].neighbor.size() == 0)  return;
	for (auto& u : g.node[v].neighbor)
	{
		
		for (int i = 0; i < g.node[u].neighbor.size(); i++)
		{
			if (g.node[u].neighbor[i] == v)
			{
				g.node[u].neighbor[i] = g.node[u].neighbor.back();
				break;
			}
		}
		g.node[u].neighbor.pop_back();
	}  
	
}


void get_hneighbor(Graph& g, int v, int h)  //在1-hop原图上做BFS计算顶点v的h-hop邻居
{

	g.node[v].hneighbor.clear();
	g.node[v].hdegree = 0;

	vector<bool> hvisited(g.node.size(), false);
	queue<int> Q1;
	hvisited[v] = true;
	for (auto& u : g.node[v].neighbor)
	{
		Q1.push(u); //先载入1-hop邻居
		hvisited[u] = true;
	}

	int d = 1;
	while (d <= h)
	{
		queue<int> Q2;  //Q2每轮循环初始为空
		while (!Q1.empty())
		{
			int u = Q1.front();
			Q1.pop();
			g.node[v].hneighbor.push_back(u);
			g.node[v].hdegree++;
			//cout << "neighbor: " << u << '\n' << endl;
			for (auto w : g.node[u].neighbor)
			{
				if (hvisited[w] == false)
				{
					Q2.push(w);
					hvisited[w] = true;
				}
			}
		}
		d++;
		Q1.swap(Q2);
	}

}


void GraphAToB(Graph g1, Graph& g2)   //因为peak和core涉及两层对图的删除操作，需要制作一个副本在core中执行
{
	g2.n = g1.n;
	g2.node.resize(vertexnum + 1);
	for (int v = 1; v <= vertexnum; ++v)
	{
		//cout << v << '\n' << endl;
		g2.node[v].deleted = g1.node[v].deleted;
		if (g1.node[v].deleted == false)
		{
			g2.node[v].neighbor = move(g1.node[v].neighbor);
		}  
		else
			g2.node[v].neighbor = { 0 };
		g2.node[v].hneighbor = { 0 };
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
                int pre_hdegree = std::max(g.node[u].hdegree, k);
                get_hneighbor(g, u, h);
                int temp = std::max(g.node[u].hdegree, k);
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

// /* 基础的higher-order core decompsition */
// void get_khcore(Graph& g, int h)
// {   
// 	int n = g.n;                     //剩余图顶点数量为n
// 	VVI bin(10000);

// 	/* 每轮hpeak decomposition删点后重新整合剩下的顶点*/
// 	for (int v = 1; v <= vertexnum; ++v)
// 	{
// 		if (g.node[v].deleted == true) continue;
// 		get_hneighbor(g, v, h);      //对于剩下的顶点重新计算其hneighbor
// 		//cout << v << " hdegree = " << g.node[v].hdegree << '\n' << endl;
// 		bin[g.node[v].hdegree].push_back(v);
// 	}
// 	// cout << " vertexnum = " << n << '\n' << endl;

// 	for (int k = 0; k < 10000; ++k)      //hcore decomposition主体部分
// 	{
// 		while (bin[k].empty() == false)
// 		{
// 			int v = bin[k].back();
// 			bin[k].pop_back();         
// 			g.node[v].hcoreness = k;
// 			deleteNode(g, v);
// 			// cout << k << " delete " << v << '\n' << endl;
// 			for (auto u : g.node[v].hneighbor)
// 			{                              
// 				int pre_hdegree = max(g.node[u].hdegree, k);  //注意目前最小bin层值只能为k
// 				get_hneighbor(g, u, h);
// 				int temp = max(g.node[u].hdegree, k); //若删点v使得u的hdegree比当前k还小，那么u的coreness就是k                
// 				bin[temp].push_back(u);
// 				for (int i = 0; i < bin[pre_hdegree].size(); i++)  //并将u从原bin层删除
// 				{
// 					if (bin[pre_hdegree][i] == u)
// 					{
// 						bin[pre_hdegree][i] = bin[pre_hdegree].back();
// 						break;
// 					}
// 				}
// 				//cout << u << ' ' << pre_hdegree << ' '<< g.node[u].hdegree << " binsize " << bin[pre_hdegree].size() << '\n' << endl;
// 				bin[pre_hdegree].pop_back();        
// 			}
// 		}
// 	}  

// 	 /*
// 	for (int v = 1; v <= vertexnum; v++)
// 	{
// 		cout << v << ' ' << g.node[v].hcoreness << '\n' << endl;
// 	}
// 	// */
// }





// // Function to perform higher-order core decomposition
// void get_khcore(Graph& g, int h) {
//     int n = g.n; // Remaining vertex count is n
//     std::vector<std::vector<int>> bin(10000); // Dynamic vector for bin

//     // Perform hpeak decomposition and reorganize remaining vertices
//     for (int v = 1; v <= vertexnum; ++v) {
//         if (g.node[v].deleted) continue;
//         get_hneighbor(g, v, h);
//         bin[g.node[v].hdegree].push_back(v);
//     }

//     for (int k = 0; k < 10000; ++k) {
//         while (!bin[k].empty()) {
//             int v = bin[k].back();
//             bin[k].pop_back();
//             g.node[v].hcoreness = k;
//             deleteNode(g, v);
            
//             for (auto u : g.node[v].hneighbor) {
//                 int pre_hdegree = std::max(g.node[u].hdegree, k);
//                 get_hneighbor(g, u, h);
//                 int temp = std::max(g.node[u].hdegree, k);
//                 bin[temp].push_back(u);

//                 for (int i = 0; i < bin[pre_hdegree].size(); i++) {
//                     if (bin[pre_hdegree][i] == u) {
//                         bin[pre_hdegree][i] = bin[pre_hdegree].back();
//                         break;
//                     }
//                 }
//                 bin[pre_hdegree].pop_back();
//             }
//         }
//     }
// }









/* H-index计算：使用计数排序，复杂度O(n) */
int hIndex(int v, Graph g, vector<int> Reach)  //注意传入为数组Reach
{
	int neighborSize = g.node[v].neighbor.size();
	if (neighborSize < 1) return 0;
	vector<int> record(neighborSize + 1, 0); //下标从0开始，所以要开“+1”数组

	for (auto& u : g.node[v].neighbor)
	{
		if (Reach[u] >= neighborSize)
			++record[neighborSize];
		else
			++record[Reach[u]];
	}

	int h = neighborSize, sum = record[neighborSize];
	while (sum < h)
	{
		sum += record[--h];
	}
	//cout << h << '\n' << endl;
	return h;
}

int ComputeH(Graph g, int v, int h, vector<int> pre_H_order)  //未验证
{
	vector<int> Reach1(vertexnum + 1, 0);
	vector<int> Reach2(vertexnum + 1, 0);
	vector<bool> pathvisited(vertexnum + 1, 0);
	queue<int> Q1;
	for (auto& u : g.node[v].neighbor)
	{
		Reach1[u] = pre_H_order[u];
		Reach2[u] = pre_H_order[u];
		Q1.push(u);
	}
	int d = 2;
	while (d <= h)
	{
		queue<int> Q2;  //Q2每轮循环初始为空
		while (Q1.empty() == false)
		{
			int u = Q1.front();
			Q1.pop();
			pathvisited[u] = true;
			for (auto w : g.node[u].neighbor)
			{
				int min = g.node[w].hcoreness;
				if (Reach2[u] < g.node[w].hcoreness) min = Reach2[u];
				if (pathvisited[w] == false || Reach2[w] < min)
				{
					Reach1[w] = min;
					Q2.push(w);
				}
			}
		}
		d++;
		Q1.swap(Q2);
		for (auto& u : g.node[v].neighbor)
		{
			
			Reach2[u] = Reach1[u];
		}
	}
	int h_result = hIndex(v, g, Reach1);
	return h_result;
}

int MyComputeH(Graph g, int v, int h)  //应该是因为原算法有点错误，我自己改了一点写的
{
	vector<int> Reach1(vertexnum + 1, 0);
	vector<int> Reach2(vertexnum + 1, 0);
	vector<bool> pathvisited(vertexnum + 1, 0);
	queue<int> Q1;
	for (auto& u : g.node[v].neighbor)
	{
		Reach1[u] = g.node[u].hcoreness;
		Reach2[u] = g.node[u].hcoreness;
		Q1.push(u); 
		pathvisited[u] = true;
	}
	int d = 2;
	while (d <= h)
	{
		queue<int> Q2;  //Q2每轮循环初始为空
		while (Q1.empty() == false)
		{
			int u = Q1.front();
			Q1.pop();
			//pathvisited[u] = true;          
			for (auto w : g.node[u].neighbor)
			{
				int min = g.node[w].hcoreness;
				if (Reach2[u] < g.node[w].hcoreness) min = Reach2[u];
				if (pathvisited[w] == false || Reach2[w] < min)
				{
					Reach1[w] = min;
					Q2.push(w);
					pathvisited[w] = true;
				}        
			}
			for (auto u : g.node[v].hneighbor)
			{
				Reach2[u] = Reach1[u];
			}
		}
		d++;
		Q1.swap(Q2);
	 
	}
	int h_result = hIndex(v, g, Reach1);
	return h_result;
}


/* Advanced Local Update：更新可能被影响的邻点 */ 
void UpdateNB_Plus(int u, int pre_cnu, Graph& g, int h)
{
	for (auto v : g.node[u].hneighbor)
	{
		//if (g.node[v].visited == true && visit_order[v] == 0) continue; //Q中的顶点集不需要更新，因为每轮Q.pop时会更新
		if (g.node[u].hcoreness >= g.node[v].hcoreness) continue;
		else
		{
			if (pre_cnu < g.node[v].hcoreness) continue;
			else
			{
				int new_cnv = MyComputeH(g, v, h);
				if (new_cnv != g.node[v].hcoreness)
				{
					int pre_cnv = g.node[v].hcoreness;
					g.node[v].hcoreness = new_cnv;
					//cout << v << " update from " << pre_cnv << " to " << g.node[v].coreness << '\n' << endl;
					UpdateNB_Plus(v, pre_cnv, g, h);
				}
			}
		}
	}
}
/* 并行的higher-order core decompsition */
void para_get_khcore(Graph& g, int h)
{
	vector<int> n0_order(vertexnum + 1, 0);
	vector<int> n_order(vertexnum + 1, 0);
	for (int v = 1; v <= vertexnum; ++v)
	{
		if (g.node[v].deleted == true) continue;
		get_hneighbor(g, v, h);  
		//cout << v << ' ' << g.node[v].hdegree << '\n' << endl;
		n0_order[v] = g.node[v].hdegree;  
	}

	bool tag = true;
	int n = 0;
	while (tag == true)
	{
		tag = false;
		n++;
		#pragma omp parallel for schedule(static) num_threads(40)
		for (int v = 1; v <= vertexnum; v++)
		{
			/* 原句：n_order[v] = ComputeH(g, v, h); */  //OpenMP内不使用复杂函数调用，要直接写
			vector<int> Reach1(vertexnum + 1, 0);  
			vector<int> Reach2(vertexnum + 1, 0);
			vector<bool> pathvisited(vertexnum + 1, 0);
			queue<int> Q1;
			for (auto& u : g.node[v].neighbor)
			{
				Reach1[u] = n0_order[u];
				Reach2[u] = n0_order[u];
				Q1.push(u);
				pathvisited[u] = true;
			}
			int d = 2;
			while (d <= h)
			{
				queue<int> Q2;  //Q2每轮循环初始为空
				while (Q1.empty() == false)
				{
					int u = Q1.front();
					Q1.pop();
					//pathvisited[u] = true;
					for (auto w : g.node[u].neighbor)
					{
						int min = n0_order[w];
						if (Reach2[u] < n0_order[w]) min = Reach2[u];
						if (pathvisited[w] == false || Reach2[w] < min)
						{
							Reach1[w] = min;
							Q2.push(w);
							pathvisited[w] = true;
						}
					}
					for (auto w : g.node[v].hneighbor)
					{
						Reach2[w] = Reach1[w];
					}
				}
				d++;
				Q1.swap(Q2);
				
			}
			/* 原句：int h_result = hIndex(v, g, Reach1); */
			int h_index = 0;
			int neighborSize = g.node[v].hneighbor.size();
			if (neighborSize < 1) h_index = 0;
			else
			{
				vector<int> record(neighborSize + 1, 0); //下标从0开始，所以要开“+1”数组
				for (auto u : g.node[v].hneighbor)
				{
					if (Reach1[u] >= neighborSize)
						++record[neighborSize];
					else
						++record[Reach1[u]];

				}
				int h = neighborSize, sum = record[neighborSize];
				while (sum < h)
				{
					sum += record[--h];
				}
				h_index = h;
			}
			n_order[v] = h_index;
			if (n_order[v] != n0_order[v]) tag = true;
			n0_order[v] = n_order[v];  //优化
		}
		for (int v = 1; v <= vertexnum; ++v)  n0_order[v] = n_order[v];
		
	}
	//cout << " order：" << n << '\n' << endl;

	for (int v = 1; v <= vertexnum; ++v)
	{
		g.node[v].hcoreness = n_order[v];
		//cout << v << ' ' << g.node[v].hcoreness << '\n' << endl;
	}
}
/* (1)迭代get_khcore算法 */
void BaselineAlgorithm(Graph& g, int h)
{
	while (g.n > 0)
	{
		
		Graph g_copy;
		GraphAToB(g, g_copy);
		clock_t StartTime, EndTime;
		StartTime = clock();
		//cout << "vertexnum = 1" << '\n' << endl;

		get_khcore(g_copy, h);  
		//para_get_khcore(g_copy, h);  
		EndTime = clock(); 
		if(g.n == vertexnum)  cout << " Core run time: " << (double)(EndTime - StartTime) / double(CLOCKS_PER_SEC) << '\n' << endl;

		/*
		for (int v = 1; v <= vertexnum; v++)  //test
		{
			if (g.node[v].deleted == true) continue;
			cout << v << " : " << g_copy.node[v].hcoreness << '\n' << endl;
		}*/

		int maxhcore = 0;
		for (int v = 1; v <= vertexnum; v++)  //找当前存在最大的hcoreness
		{
			if (g.node[v].deleted == true) continue;
			if (maxhcore < g_copy.node[v].hcoreness)  maxhcore = g_copy.node[v].hcoreness;
		}
		//cout << "maxhcore = " << maxhcore << '\n' << endl;
		if (maxhcore == 0)  //maxcore变为0后可以直接退出
		{
			for (int u = 1; u <= vertexnum; u++)
			{
				if (g.node[u].deleted == true) continue;
				g.node[u].hpeakness = 0;
			}
			break;
		}
		queue<int> Q;
		for (int i = 1; i <= vertexnum; i++)
		{
			if (g.node[i].deleted == true) continue;
			if (g_copy.node[i].hcoreness == maxhcore) Q.push(i);
		}
		while (Q.empty() == false)
		{
			int v = Q.front();  //v为每轮需要删除的顶点
			Q.pop();
			g.node[v].hpeakness = maxhcore;
			deleteNode(g, v);
		}    
		//cout << g.n << '\n' << endl;
		//free(&g_copy);  //这个用起来有问题
	}
}



/* (2)迭代para_get_khcore算法 */
void AllParalleAlgorithm(Graph& g, int h)
{
	while (g.n > 0)
	{

		Graph g_copy;
		GraphAToB(g, g_copy);
		struct timeval t1, t2;
		double cost;
		gettimeofday(&t1, NULL);
		para_get_khcore(g_copy, h);
		gettimeofday(&t2, NULL);
		cost = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
		if (g.n == vertexnum)  cout << " Core run time: " << cost << '\n' << endl;

		/*
		for (int v = 1; v <= vertexnum; v++)  //test
		{
			if (g.node[v].deleted == true) continue;
			cout << v << " : " << g_copy.node[v].hcoreness << '\n' << endl;
		}*/

		//cout << "vertexnum = " << g.n << '\n' << endl;
		int maxhcore = 0;
		for (int v = 1; v <= vertexnum; v++)  //找当前存在最大的hcoreness
		{
			if (g.node[v].deleted == true) continue;
			if (maxhcore < g_copy.node[v].hcoreness)  maxhcore = g_copy.node[v].hcoreness;
		}
		//cout << "maxhcore = " << maxhcore << '\n' << endl;
		if (maxhcore == 0)  //maxcore变为0后可以直接退出
		{
			for (int u = 1; u <= vertexnum; u++)
			{
				if (g.node[u].deleted == true) continue;
				g.node[u].hpeakness = 0;
			}
			break;
		}
		queue<int> Q;
		for (int i = 1; i <= vertexnum; i++)
		{
			if (g.node[i].deleted == true) continue;
			if (g_copy.node[i].hcoreness == maxhcore) Q.push(i);
		}
		while (Q.empty() == false)
		{
			int v = Q.front();  //v为每轮需要删除的顶点
			Q.pop();
			g.node[v].hpeakness = maxhcore;
			deleteNode(g, v);
		}
		//cout << g.n << '\n' << endl;

	}
}

/* (3/4)仅一次（para）get_khcore + Local Update算法 */
void LocalUpdateAlgorithm(Graph& g, int h)  //局部并行方案
{
	Graph g_copy;
	GraphAToB(g, g_copy);
	struct timeval t1, t2;
	double cost;
	gettimeofday(&t1, NULL);
	//get_khcore(g_copy, h);     //(3)基础get_khcore
	para_get_khcore(g_copy, h);  //(4)并行优化para_get_khcore
	gettimeofday(&t2, NULL);
	cost = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
	//if (g.n == vertexnum)  cout << " Core run time: " << cost << '\n' << endl;

	for (int v = 1; v <= vertexnum; ++v)
	{
		get_hneighbor(g, v, h);  //刚开始需要算一遍hneighbor，后续在此基础上更新即可
		g.node[v].hcoreness = g_copy.node[v].hcoreness;  //仅需要将算得的hcoreness记录下
	}

	while (g.n > 0)
	{     
		/*
		for (int v = 1; v <= vertexnum; v++)  //test
		{
			if (g.node[v].deleted == true) continue;
			cout << v << " : " << g_copy.node[v].hcoreness << '\n' << endl;
		}*/
		
		//cout << "vertexnum = " << g.n << '\n' << endl;
		int maxhcore = 0;
		for (int v = 1; v <= vertexnum; v++)  //找当前存在最大的coreness
		{
			if (g.node[v].deleted == true) continue;
			if (maxhcore < g.node[v].hcoreness)  maxhcore = g.node[v].hcoreness;
		}
		//cout << "maxhcore = " << maxhcore << '\n' << endl;
		if (maxhcore == 0)  //maxcore变为0后可以直接退出
		{
			for (int u = 1; u <= vertexnum; u++)
			{
				if (g.node[u].deleted == true) continue;
				g.node[u].hpeakness = 0;
			}
			break;
		}

		for (int i = 1; i <= vertexnum; i++)  //每轮初始化visited数组
		{
			g.node[i].visited = false;
			//visit_order[i] = 0;
		}

		vector<int> cand;
		for (int i = 1; i <= vertexnum; i++)  //找到当前hcoreness最大的顶点i，将其cn赋值给pn并删除这些顶点
		{
			if (g.node[i].deleted == true) continue;
			if (g.node[i].hcoreness == maxhcore)  //准备删除顶点i
			{
				g.node[i].hpeakness = maxhcore;
				//g.node[i].deleted = true;
				//cout << "delete " << i << '\n' << endl;
				
				for (auto j : g.node[i].hneighbor)
				{
					if (g.node[j].hcoreness < maxhcore)  //改成＜试试看，原为≠
					{
						if (g.node[j].visited == false)
						{
							cand.push_back(j); //将被删除的顶点i的邻点j加入cand[]“可能被影响的点集”
							//cout << "cand: " << j << '\n' << endl;
							g.node[j].visited = true;
						}
					}
				}
				deleteNode(g, i);
			}
		}

		for (int i = 0; i < cand.size(); i++)  //先更新一遍可能受影响顶点的hneighbor(只有这些顶点的hneighbor会改变)
		{
			get_hneighbor(g, cand[i], h);
		}

		/* 并行的对cand中顶点的coreness进行更新 */
		#pragma omp parallel for schedule(dynamic) num_threads(30)
		for (int i = 0; i < cand.size(); i++)
		{
			int u = cand[i];
			//printf("%d is from thread = %d\n", u, omp_get_thread_num());

			queue<int> Q;  //用队列实现递归转迭代
			Q.push(u);
			while (Q.empty() == false)
			{
				int v = Q.front();
				Q.pop();
				/* 原句：int new_cnv = MyComputeH(g, v, h); */
				vector<int> Reach1(vertexnum + 1, 0);
				vector<int> Reach2(vertexnum + 1, 0);
				vector<bool> pathvisited(vertexnum + 1, false);
				queue<int> Q1;
				for (auto w : g.node[v].neighbor)
				{
					Reach1[w] = g.node[w].hcoreness;
					Reach2[w] = g.node[w].hcoreness;
					Q1.push(w);
					pathvisited[w] = true;
				}
				int d = 2;
				while (d <= h)
				{
					queue<int> Q2;  //Q2每轮循环初始为空
					while (Q1.empty() == false)
					{
						int q = Q1.front();
						Q1.pop();
						//pathvisited[q] = true;
						for (auto w : g.node[q].neighbor)
						{
							int min = 0;
							if (Reach2[q] < g.node[w].hcoreness) min = Reach2[q];
							else  min = g.node[w].hcoreness;

							if (pathvisited[w] == false || Reach2[w] < min)
							{
								Reach1[w] = min;
								Q2.push(w);
								pathvisited[w] = true;
							}
						}
						
						for (auto w : g.node[v].hneighbor)
						{
							Reach2[w] = Reach1[w];
						}
					}
					d++;
					Q1.swap(Q2);
					/*
					for (auto w : g.node[v].neighbor)
					{
						Reach2[w] = Reach1[w];
					}*/
				}
				
				/* 原句：int h_result = hIndex(v, g, Reach1); */
				int h_index = 0;
				int neighborSize = g.node[v].hneighbor.size();
				if (neighborSize < 1) h_index = 0;
				else
				{
					vector<int> record(neighborSize + 1, 0); //下标从0开始，所以要开“+1”数组
					for (auto w : g.node[v].hneighbor)
					{
						if (Reach1[w] >= neighborSize)
							++record[neighborSize];
						else
							++record[Reach1[w]];
					}
					int h = neighborSize, sum = record[neighborSize];
					while (sum < h)
					{
						sum += record[--h];
					}
					h_index = h;
				}
				int new_cnv = h_index;
				int pre_cnv = g.node[v].hcoreness;

				if (new_cnv != pre_cnv)  //若v点的coreness变化了，将消息传递给其所有邻点w(此处操作为将w加入待处理点集Q)
				{
					#pragma omp atomic write
					g.node[v].hcoreness = new_cnv;  //赋值处需互斥进行写新数据

					/* 原句：UpdateNB_Plus(v, pre_cnv, g); */ 
					
					for (auto w : g.node[v].hneighbor)  //对于higher order还是可以使用的优化
					{
						//Q.push(w);
						// /* OPT: 可以提前判定跳过的点
						if (new_cnv > g.node[w].hcoreness) continue;
						else
						{
							if (pre_cnv < g.node[w].hcoreness)  continue;
							else  Q.push(w);
						}
						// */
					}
					

				}
			}
			//visit_order[u] = 1;  //标记Q中的顶点访问顺序，使得不被重复计算
			//cout << new_cnu << '\n' << endl;        
		}

	}
}

void get_shell(Graph& g, int v, int h, vector<int> ispeak)  //在1-hop原图上做BFS计算顶点v的shell层h-hop邻居
{

	g.node[v].shell.clear();
	g.node[v].sdegree = 0;

	vector<bool> hvisited(vertexnum + 1, false);
	queue<int> Q1;
	hvisited[v] = true;
	for (auto& u : g.node[v].neighbor)
	{
		if (g.node[u].hcoreness == g.node[v].hcoreness && ispeak[u] == 1)  //u为v的(k,h)-shell层1-hop邻居
		{
			Q1.push(u); //先载入1-hop邻居
			hvisited[u] = true;
		}
	}

	int d = 1;
	while (d <= h)
	{
		queue<int> Q2;  //Q2每轮循环初始为空
		while (Q1.empty() == false)
		{
			int u = Q1.front();
			Q1.pop();
			g.node[v].shell.push_back(u);
			g.node[v].sdegree++;
			//cout << "neighbor: " << u << '\n' << endl;
			for (auto w : g.node[u].neighbor)
			{
				if (hvisited[w] == false && g.node[w].hcoreness == g.node[v].hcoreness && ispeak[w] == 1)
				{
					Q2.push(w);
					hvisited[w] = true;
				}
			}
		}
		d++;
		Q1.swap(Q2);
	}

}

void EarlyCovergeOpt(Graph& g, int h, vector<int>& ispeak)
{
	//此前已获取全局顶点的hcoreness
	for (int v = 1; v <= vertexnum; v++)  //初始计算k-shell子图中每个顶点的shell度数
	{
		for (auto u : g.node[v].hneighbor)
		{
			if (g.node[u].hcoreness == g.node[v].hcoreness)
			{
				g.node[v].shell.push_back(u);
				g.node[v].sdegree++;
			}
		}
	}

	queue<int> Q;  // Q存储不满足shell内sdegree小于coreness的顶点
	for (int v = 1; v <= vertexnum; v++)
	{
		if (g.node[v].sdegree < g.node[v].hcoreness)
		{
			Q.push(v);
		}
	}

	while (Q.empty() == false)
	{
		int u = Q.front();
		Q.pop();
		ispeak[u] = 0;  //先把u删除，再重新计算它shell层内每个点的h-hop邻居
		for (auto v : g.node[u].shell)
		{
			get_shell(g, v, h, ispeak);  //重新计算它shell层内每个点的h-hop邻居
			if (g.node[v].sdegree < g.node[v].hcoreness)  Q.push(v);
		}       
	}
	 /*
	for (int v = 1; v <= vertexnum; v++)
	{
		if (ispeak[v] == 1)  cout << v << '\n' << endl;
	}
	// */

}

/* (5)提前收敛优化 + (3/4)仅一次（para）get_khcore + Local Update算法 */
void PlusLocalAlgorithm(Graph& g, int h, vector<int> ispeak)  //局部并行方案
{
	Graph g_copy;
	GraphAToB(g, g_copy);
	struct timeval t1, t2;
	double cost;
	gettimeofday(&t1, NULL);
	//get_khcore(g_copy, h);     //(3)基础get_khcore
	para_get_khcore(g_copy, h);  //(4)并行优化para_get_khcore
	gettimeofday(&t2, NULL);
	cost = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
	if (g.n == vertexnum)  cout << " Core run time: " << cost << '\n' << endl;
	int earlynum = 0;
	for (int v = 1; v <= vertexnum; ++v)
	{
		get_hneighbor(g, v, h);  //刚开始需要算一遍hneighbor，后续在此基础上更新即可
		g.node[v].hcoreness = g_copy.node[v].hcoreness;  //仅需要将算得的hcoreness记录下
		/* OPT：预处理找提前收敛点 */
		if (ispeak[v] == 1)  
		{
			g.node[v].hpeakness = g.node[v].hcoreness;
			earlynum++;
		}  
	}
	cout << " early pruning number: " << earlynum << '\n' << endl;

	while (g.n > 0)
	{
		//cout << "vertexnum = " << g.n << '\n' << endl;
		int maxhcore = 0;
		for (int v = 1; v <= vertexnum; v++)  //找当前存在最大的coreness
		{
			if (g.node[v].deleted == true) continue;
			if (maxhcore < g.node[v].hcoreness)  maxhcore = g.node[v].hcoreness;
		}
		//cout << "maxhcore = " << maxhcore << '\n' << endl;
		if (maxhcore == 0)  //maxcore变为0后可以直接退出
		{
			for (int u = 1; u <= vertexnum; u++)
			{
				if (g.node[u].deleted == true) continue;
				g.node[u].hpeakness = 0;
			}
			break;
		}

		for (int i = 1; i <= vertexnum; i++)  //每轮初始化visited数组
		{
			g.node[i].visited = false;
			//visit_order[i] = 0;
		}

		vector<int> cand;
		for (int i = 1; i <= vertexnum; i++)  //找到当前hcoreness最大的顶点i，将其cn赋值给pn并删除这些顶点
		{
			if (g.node[i].deleted == true) continue;
			if (g.node[i].hcoreness == maxhcore)  //准备删除顶点i
			{
				g.node[i].hpeakness = maxhcore;
				//g.node[i].deleted = true;
				//cout << "delete " << i << '\n' << endl;

				for (auto j : g.node[i].hneighbor)
				{
					if (ispeak[j] == 1)  continue;  //k-shell优化
					if (g.node[j].hcoreness < maxhcore)  //改成＜试试看，原为≠
					{
						if (g.node[j].visited == false)
						{
							cand.push_back(j); //将被删除的顶点i的邻点j加入cand[]“可能被影响的点集”
							//cout << "cand: " << j << '\n' << endl;
							g.node[j].visited = true;
						}
					}
				}
				deleteNode(g, i);
			}
		}

		for (int i = 0; i < cand.size(); i++)  //先更新一遍可能受影响顶点的hneighbor(只有这些顶点的hneighbor会改变)
		{
			get_hneighbor(g, cand[i], h);
		}

		/* 并行的对cand中顶点的coreness进行更新 */
		#pragma omp parallel for schedule(dynamic,40) 
		for (int i = 0; i < cand.size(); i++)
		{
			int u = cand[i];
			if (ispeak[u] == 1)  continue;  //k-shell优化
			//printf("%d is from thread = %d\n", u, omp_get_thread_num());
			queue<int> Q;  //用队列实现递归转迭代
			Q.push(u);
			while (Q.empty() == false)
			{
				int v = Q.front();
				Q.pop();
				/* 原句：int new_cnv = MyComputeH(g, v, h); */
				vector<int> Reach1(vertexnum + 1, 0);
				vector<int> Reach2(vertexnum + 1, 0);
				vector<bool> pathvisited(vertexnum + 1, false);
				queue<int> Q1;
				for (auto w : g.node[v].neighbor)
				{
					Reach1[w] = g.node[w].hcoreness;
					Reach2[w] = g.node[w].hcoreness;
					Q1.push(w);
					pathvisited[w] = true;
				}
				int d = 2;
				while (d <= h)
				{
					queue<int> Q2;  //Q2每轮循环初始为空
					while (Q1.empty() == false)
					{
						int q = Q1.front();
						Q1.pop();
						//pathvisited[q] = true;
						for (auto w : g.node[q].neighbor)
						{
							int min = 0;
							if (Reach2[q] < g.node[w].hcoreness) min = Reach2[q];
							else  min = g.node[w].hcoreness;

							if (pathvisited[w] == false || Reach2[w] < min)
							{
								Reach1[w] = min;
								Q2.push(w);
								pathvisited[w] = true;
							}
						}

						for (auto w : g.node[v].hneighbor)
						{
							Reach2[w] = Reach1[w];
						}
					}
					d++;
					Q1.swap(Q2);
					/*
					for (auto w : g.node[v].neighbor)
					{
						Reach2[w] = Reach1[w];
					}*/
				}

				/* 原句：int h_result = hIndex(v, g, Reach1); */
				int h_index = 0;
				int neighborSize = g.node[v].hneighbor.size();
				if (neighborSize < 1) h_index = 0;
				else
				{
					vector<int> record(neighborSize + 1, 0); //下标从0开始，所以要开“+1”数组
					for (auto w : g.node[v].hneighbor)
					{
						if (Reach1[w] >= neighborSize)
							++record[neighborSize];
						else
							++record[Reach1[w]];
					}
					int h = neighborSize, sum = record[neighborSize];
					while (sum < h)
					{
						sum += record[--h];
					}
					h_index = h;
				}
				int new_cnv = h_index;
				int pre_cnv = g.node[v].hcoreness;

				if (new_cnv != pre_cnv)  //若v点的coreness变化了，将消息传递给其所有邻点w(此处操作为将w加入待处理点集Q)
				{
					#pragma omp critical
					g.node[v].hcoreness = new_cnv;  //赋值处需互斥进行写新数据

					/* 原句：UpdateNB_Plus(v, pre_cnv, g); */
					for (auto w : g.node[v].hneighbor)  //对于higher order还是可以使用的优化
					{
						if (ispeak[w] == 1)  continue;  //k-shell优化
						if (new_cnv > g.node[w].hcoreness) continue;
						else
						{
							if (pre_cnv < g.node[w].hcoreness)  continue;
							else  Q.push(w);
						}
					}
				}
			}
			//visit_order[u] = 1;  //标记Q中的顶点访问顺序，使得不被重复计算
			//cout << new_cnu << '\n' << endl;        
		}

	}
}


/* 计算Core模块度 */
void CoreModularity(Graph& g, int h)
{
	para_get_khcore(g, h);
  
	double bound = 0.05;
	for(int k = 1; k <= 6; k++)
	{
		int count = vertexnum * bound * k;
		int maxhcore = 0;
		VVI bin(g.node.size()); 
		for (int v = 1; v <= vertexnum; v++)  //找当前存在最大的coreness
		{
			bin[g.node[v].hcoreness].push_back(v);
			if (maxhcore < g.node[v].hcoreness)  maxhcore = g.node[v].hcoreness;
		}
		for (int i = maxhcore; i >= 0; i--)
		{
			for (int j = 0; j < bin[i].size(); j++) 
			{
				g.node[bin[i][j]].kbound = 1;
				count--;
				if(count == 0) break;
			}
			if(count == 0) break;
		}
	
		int inedgenum = 0;
		int indegreesum = 0;
		for (int v = 1; v <= vertexnum; v++)
		{
			if (g.node[v].kbound == 1)
			{
				indegreesum += g.node[v].neighbor.size();
				for (auto& u : g.node[v].neighbor)
				{
					if (u > v && g.node[v].kbound == 1)  inedgenum++;  //边不重复计算
				}
			}
		}

		double coreresult = 0;
		double m = edgenum;
		coreresult = (inedgenum / m) - (indegreesum / (2 * m)) * (indegreesum / (2 * m));
		if(k == 1) cout << " kh-core 模块度：" << '\n' << endl;
		cout << coreresult << '\n' << endl;
	}

}

/* 计算Peak模块度 */
void PeakModularity(Graph& g, int h)
{
	LocalUpdateAlgorithm(g, h);
   
	for (int v = 1; v <= vertexnum; v++) g.node[v].kbound = 0;  //清空core计算时的标记
	double bound = 0.05;
	for(int k = 1; k <= 6; k++)
	{
		int count = vertexnum * bound * k;
		int maxhpeak = 0;
		VVI bin(g.node.size()); 
		for (int v = 1; v <= vertexnum; v++)  //找当前存在最大的peakness
		{
			bin[g.node[v].hpeakness].push_back(v);
			if (maxhpeak < g.node[v].hpeakness)  maxhpeak = g.node[v].hpeakness;
		}
		for (int i = maxhpeak; i >= 0; i--)
		{
			for (int j = 0; j < bin[i].size(); j++) 
			{
				g.node[bin[i][j]].kbound = 1;
				count--;
				if(count == 0) break;
			}
			if(count == 0) break;
		}
	
		int inedgenum = 0;
		int indegreesum = 0;
		for (int v = 1; v <= vertexnum; v++)
		{
			if (g.node[v].kbound == 1)
			{
				indegreesum += g.node[v].neighbor.size();
				for (auto& u : g.node[v].neighbor)
				{
					if (u > v && g.node[v].kbound == 1)  inedgenum++;  //边不重复计算
				}
			}
		}

		double peakresult = 0;
		double m = edgenum;
		peakresult = (inedgenum / m) - (indegreesum / (2 * m)) * (indegreesum / (2 * m));
		if(k == 1) cout << " kh-peak 模块度：" << '\n' << endl;
		cout << peakresult << '\n' << endl;
	}

}


int main(int argc, char* argv[])
{

	//int h = 4;
	for (int h = 2; h <= 3; h++) { //start from h = 2 
		// cout << " -----   H = " << h << '\n' << endl;
		// Graph g;
		// /* Open the input file Yeast facebook Vidal sister_cities Gnutella09 caida collaboration douban amazon*/
		// ifstream inputFile("data/Yeast.txt");
		// if (!inputFile) {
		// 	cerr << "Input file could not be opened!" << endl;
		// 	exit(EXIT_FAILURE);
		// }

		// inputFile >> vertexnum >> edgenum;
		// g.n = vertexnum;
		// g.m = edgenum;
		// g.node.resize(vertexnum + 1);

		// for (int i = 1; i <= edgenum; i++)  //数据集从1开始编号
		// {
		// 	int u, v;
		// 	inputFile >> u >> v;
		// 	if (u == v) continue;           //遇到loop需跳过，eg:数据集Yeast有loop
		// 	g.node[u].neighbor.push_back(v);
		// 	g.node[v].neighbor.push_back(u);

		// }
		// inputFile.close();
		// cout << " 图g加载成功!" << '\n' << endl;


		cout << " -----   H = " << h << '\n' << endl;
		Graph g;

		ifstream inputFile("data/Vidal.txt");
		if (!inputFile) {
			cerr << "Input file could not be opened!" << endl;
			exit(EXIT_FAILURE);
		}

		int maxVertexIndex = 0;
		int maxEdgeIndex = 0;

		inputFile >> vertexnum >> edgenum;
		for (int i = 1; i <= edgenum; i++)  {
			int u, v;
			inputFile >> u >> v;
			if (u == v) continue;
			maxEdgeIndex ++;
			maxVertexIndex = max(maxVertexIndex, max(u, v));
		}
		inputFile.close();

		g.n = maxVertexIndex;
		g.m = edgenum;
		g.node.resize(maxVertexIndex + 1);

		inputFile.open("data/Vidal.txt");
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
		cout << " Graph g loaded successfully!" << '\n' << endl;
		
		/* 计算聚类系数 */
		//Coefficient(g, h);

		/* 计算模块度 */
		//CoreModularity(g, h);
		PeakModularity(g, h);
		


		// Baseline算法：剥离法 
				//迭代get_khcore算法
		// struct timeval start, end;
		// double timeuse;
		// gettimeofday(&start, NULL);

		// BaselineAlgorithm(g, h);

		// gettimeofday(&end, NULL);
		// timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;
		// cout << " baseline run time: " << timeuse << '\n' << endl;
		





		// //  迭代para_get_khcore算法
		// struct timeval start, end;
		// // double timeuse;
		// gettimeofday(&start, NULL);

		// // All Paralle算法：全并行法
		// AllParalleAlgorithm(g, h);

		// gettimeofday(&end, NULL);
		// timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;
		// cout << " parallel run time: " << timeuse << '\n' << endl;
		// // 





		// //  // (para)get_khcore + local update算法
		// struct timeval start, end;
		// double timeuse;
		// gettimeofday(&start, NULL);

		// // Local Update算法：局部并行更新法 
		// LocalUpdateAlgorithm(g, h);

		// gettimeofday(&end, NULL);
		// timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;
		// cout << " local run time: " << timeuse << '\n' << endl;
		// // 


		//  ////////(提前收敛)(para)get_khcore + local update算法
		// vector<int> ispeak(vertexnum + 1, 1);  
		// para_get_khcore(g, h);
		// //get_khcore(g, h);
		// EarlyCovergeOpt(g, h, ispeak);
		// int earlynum = 0;
		// for (int v = 1; v <= vertexnum; ++v)
		// {
		//     if (ispeak[v] == 1)  
		//     {
		//         g.node[v].hpeakness = g.node[v].hcoreness;
		//         earlynum++;
		//     }  
		// }
		// cout << " early pruning number: " << earlynum << '\n' << endl;
		// //exit(0);
	
		// struct timeval start, end;
		// double timeuse;
		// gettimeofday(&start, NULL);

		// PlusLocalAlgorithm(g, h, ispeak);

		// gettimeofday(&end, NULL);
		// timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;
		// cout << " local_plus run time: " << timeuse << '\n' << endl;
		

		/* Output the Result */
		//ofstream outputFile("BaselineResult.txt");
		//ofstream outputFile("ParallelResult.txt");
		//ofstream outputFile("LocalResult.txt");
		//ofstream outputFile("PlusLocalResult.txt");

		
		// if (!outputFile) {
		//     cerr << "Output file could not be opened!" << endl;
		//     exit(EXIT_FAILURE);
		// }
		// for (int i = 1; i <= vertexnum; i++) {
		//     outputFile << i << ' ' << g.node[i].hpeakness << '\n';
		// }
		// outputFile.close();
		// 
	}


		
	return 0;
}
