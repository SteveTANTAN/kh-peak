//g++ localupdate.cpp -fopenmp
#include<omp.h> 
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


using namespace std;

int vertexnum = 0, edgenum = 0;
int compute_num = 0;
int compute_num_real = 0;
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
	for (auto u : g.node[v].neighbor)
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
	for (auto u : g.node[v].neighbor)
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


/* 并行的higher-order core decompsition */
void para_get_khcore(Graph& g, int h)
{
	vector<int> n0_order(g.node.size() + 1, 0);
	vector<int> n_order(g.node.size() + 1, 0);
	#pragma omp parallel for schedule(dynamic) num_threads(40)
	for (int v = 1; v <= g.node.size(); ++v)
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
		#pragma omp parallel for schedule(dynamic) num_threads(40)
		for (int v = 1; v <= g.node.size(); v++)
		{
			
			/* 原句：n_order[v] = ComputeH(g, v, h); */  //OpenMP内不使用复杂函数调用，要直接写
			vector<int> Reach1(vertexnum + 1, 0);  
			vector<int> Reach2(vertexnum + 1, 0);
			vector<bool> pathvisited(vertexnum + 1, false);
			queue<int> Q1;
			//printf("test1111 %c\n", v);
			//cout << v << ' '  << '\n' << endl;


			for (auto u : g.node[v].neighbor)
			{
				Reach1[u] = n0_order[u];
				Reach2[u] = n0_order[u];
				Q1.push(u);
				pathvisited[u] = true;
			}
			//printf("test1112 \n");

			int d = 2;
			while (d <= h)
			{
				queue<int> Q2;  //Q2每轮循环初始为空
				while (!Q1.empty())
				{
					int u = Q1.front();
					Q1.pop();
					//pathvisited[u] = true;
					//printf("test2.10 \n");

					for (auto w : g.node[u].neighbor)
					{

					
						if (pathvisited[w] == false || Reach2[w] < std::min(n0_order[w], Reach2[u]))
						{
							Reach1[w] = std::min(n0_order[w], Reach2[u]);
							Q2.push(w);
							pathvisited[w] = true;
						}
						//printf("test2.13 \n");

					}
					//printf("test2.2 \n");

					for (auto w : g.node[v].hneighbor)
					{
						Reach2[w] = Reach1[w];
					}
					//printf("test2.3 \n");

				}
				d++;
				Q1.swap(Q2);
				
			}
			/* 原句：int h_result = hIndex(v, g, Reach1); */
			//printf("test1113 \n");

			int h_index = 0;
			int neighborSize = g.node[v].hneighbor.size();
			if (neighborSize < 1) h_index = 0;
			else
			{
				vector<int> record(neighborSize + 1, 0); //下标从0开始，所以要开“+1”数组
				//printf("test1114 \n");
				
				for (auto u : g.node[v].hneighbor)
				{
					if (Reach1[u] >= neighborSize)
						++record[neighborSize];
					else
						++record[Reach1[u]];

				}
				int h = neighborSize, sum = record[neighborSize];
				//printf("test1115 \n");

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
		// gaixiancheng
        #pragma omp parallel for schedule(dynamic) num_threads(8)
        //#pragma omp parallel for schedule(static) num_threads(20)
		for (int i = 0; i < cand.size(); i++)
		{
            //printf("Thread: %d, iteration: %d\n", omp_get_thread_num(),i);
            
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


							if (pathvisited[w] == false || Reach2[w] < std::min(g.node[w].hcoreness, Reach2[q]))
							{
								Reach1[w] = std::min(g.node[w].hcoreness, Reach2[q]);
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
						compute_num ++;
						if (ispeak[w] == 1)  continue;  //k-shell优化
						if (new_cnv > g.node[w].hcoreness) continue;
						else
						{
							if (pre_cnv < g.node[w].hcoreness)  continue;
							else {
								Q.push(w); compute_num_real ++;
							}
						}

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
	// #pragma omp parallel for schedule(dynamic) num_threads(40)
	for (auto u : g.node[v].neighbor)
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
	printf("test1\n");
	gettimeofday(&t2, NULL);
	cost = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
	//if (g.n == vertexnum)  cout << " Core run time: " << cost << '\n' << endl;

	for (int v = 1; v <= vertexnum; ++v)
	{
		get_hneighbor(g, v, h);  //刚开始需要算一遍hneighbor，后续在此基础上更新即可
		g.node[v].hcoreness = g_copy.node[v].hcoreness;  //仅需要将算得的hcoreness记录下
	}
	printf("test2\n");
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
						compute_num ++;
						if (new_cnv > g.node[w].hcoreness) continue;
						else
						{
							if (pre_cnv < g.node[w].hcoreness)  continue;
							else {
								Q.push(w);
								compute_num_real ++;
							}
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


void EarlyCovergeOpt(Graph& g, int h, vector<int>& ispeak)
{
    // Calculate the sdegree in parallel
    #pragma omp parallel for num_threads(40)
    for (int v = 1; v <= vertexnum; v++)
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

    queue<int> Q;
	printf("333\n");

    // Fill Q in parallel
    #pragma omp parallel for num_threads(40)
    for (int v = 1; v <= vertexnum; v++)
    {
        if (g.node[v].sdegree < g.node[v].hcoreness)
        {
            #pragma omp critical
            {
                Q.push(v);
            }
        }
    }
	printf("444\n");

    // Process the queue
    while (!Q.empty())
    {
        int u;
        #pragma omp critical
        {
            u = Q.front();
            Q.pop();
        }

        ispeak[u] = 0;

        #pragma omp parallel for num_threads(40)
        for (auto v : g.node[u].shell)
        {
            get_shell(g, v, h, ispeak);
            if (g.node[v].sdegree < g.node[v].hcoreness)  
            {
                #pragma omp critical
                {
                    Q.push(v);
                }
            }
        }       
    }
}

// void EarlyCovergeOpt(Graph& g, int h, vector<int>& ispeak)
// {
// 	//此前已获取全局顶点的hcoreness
// 	for (int v = 1; v <= vertexnum; v++)  //初始计算k-shell子图中每个顶点的shell度数
// 	{
// 		for (auto u : g.node[v].hneighbor)
// 		{
// 			if (g.node[u].hcoreness == g.node[v].hcoreness)
// 			{
// 				g.node[v].shell.push_back(u);
// 				g.node[v].sdegree++;
// 			}
// 		}
// 	}
// 	printf("333\n");
// 	queue<int> Q;  // Q存储不满足shell内sdegree小于coreness的顶点
// 	for (int v = 1; v <= vertexnum; v++)
// 	{
// 		if (g.node[v].sdegree < g.node[v].hcoreness)
// 		{
// 			Q.push(v);
// 		}
// 	}
// 	printf("444\n");

// 	while (Q.empty() == false)
// 	{
// 		int u = Q.front();
// 		Q.pop();
// 		ispeak[u] = 0;  //先把u删除，再重新计算它shell层内每个点的h-hop邻居
		
// 		for (auto v : g.node[u].shell)
// 		{
// 			get_shell(g, v, h, ispeak);  //重新计算它shell层内每个点的h-hop邻居
// 			if (g.node[v].sdegree < g.node[v].hcoreness)  Q.push(v);
// 		}       
// 	}
// 	 /*
// 	for (int v = 1; v <= vertexnum; v++)
// 	{
// 		if (ispeak[v] == 1)  cout << v << '\n' << endl;
// 	}
// 	// */

// }





int main(int argc, char* argv[])
{


    //for (int h = 2; h <= 4; h++) {
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
	//vertexnum = maxVertexIndex;
    g.n = maxVertexIndex;
    g.m = edgenum;
    g.node.resize(maxVertexIndex + 1);

    inputFile.open(filePath);
    if (!inputFile) {
        cerr << "Input file could not be opened!" << endl;
        exit(EXIT_FAILURE);
    }
	printf("1111\n");
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
    
    vector<int> ispeak(vertexnum + 1, 1);  
	printf("222\n");

    para_get_khcore(g, h);
	printf("222\n");

    //get_khcore(g, h);
    EarlyCovergeOpt(g, h, ispeak);
    int earlynum = 0;
    for (int v = 1; v <= vertexnum; ++v)
    {
        if (ispeak[v] == 1)  
        {
            g.node[v].hpeakness = g.node[v].hcoreness;
            earlynum++;
        }  
    }
    cout << " early pruning number: " << earlynum << '\n' << endl;
    //exit(0);



    struct timeval start, end;
    double timeuse;
    gettimeofday(&start, NULL);

	//LocalUpdateAlgorithm(g, h);

    PlusLocalAlgorithm(g, h, ispeak);

    gettimeofday(&end, NULL);
    timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;
    cout << " local_plus run time: " << timeuse << '\n' << endl;

    cout << " compute_num_real: " << compute_num_real << '\n' << endl;
    //cout << " compute_prunning: " << compute_num-compute_num_real << '\n' << endl;
	cout << " compute_num: " << compute_num << '\n' << endl;
    
    ofstream outputFile("PlusLocalResult.txt");

    
    if (!outputFile) {
        cerr << "Output file could not be opened!" << endl;
        exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= vertexnum; i++) {
        outputFile << i << ' ' << g.node[i].hpeakness << '\n';
    }
    outputFile.close();
		
    //}
	
	return 0;
}

