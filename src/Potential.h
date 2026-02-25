#pragma once

#include <queue>
#include <tuple>

#include "TransportationData.h"

#define INF 100000007

using namespace std;

class Potential
{
private:

    static bool isOptimal(ProblemData pb_data, Approximation  apr,vector<int>& u,vector<int>& v){
        size_t rows = pb_data.numSuppliers;
        size_t columns = pb_data.numConsumers;

        for(int i = 0; i < rows; i++){
            for(int j = 0; j < columns; j++){
                if(apr.isEmpty[i][j]){
                    if (v[j] - u[i] > pb_data.costs[i][j]){
                        return false; //решение, увы, неоптимальное
                    }
                }
            }
        }

        return true; //а вот здесь - порядок
    };



    static bool calculatePotentialsBFS(
        ProblemData pb_data, Approximation apr,
        vector<int>& u,
        vector<int>& v) {
        
        int m = apr.numSuppliers; 
        int n = apr.numConsumers;

        u.assign(m, INF);
        v.assign(n, INF);
        
        vector<vector<pair<int, int>>> graph(m + n);
        
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (apr.isEmpty[i][j] == false) {
                    graph[i].push_back({m + j, pb_data.costs[i][j]});
                    graph[m + j].push_back({i, pb_data.costs[i][j]});
                }
            }
        }
        

        queue<int> q;
        vector<bool> visited(m + n, false);
        
        u[0] = 0;
        q.push(0);
        visited[0] = true;
        
        while (!q.empty()) {
            int current = q.front();
            q.pop();
            
            for (auto& edge : graph[current]) {
                int next = edge.first;
                int edgeCost = edge.second;
                
                if (!visited[next]) {
                    visited[next] = true;
                    
                    if (current < m) {
                        int j = next - m;
                        v[j] = edgeCost - u[current];
                        q.push(next);
                    } else {
        
                        int j = current - m;
                        int i = next;
                        u[i] = edgeCost - v[j];
                        q.push(next);
                    }
                }
            }
        }
        
        for (int i = 0; i < m; i++) {
            if (!visited[i]) return false;
        }
        for (int j = 0; j < n; j++) {
            if (!visited[m + j]) return false;
        }
        
        return true;
    }

public:

    static void getOptimal(ProblemData pb_data, Approximation  apr){

        vector<int> v;
        vector<int> u;

        calculatePotentialsBFS(pb_data,apr,u,v);

        if(isOptimal(pb_data,apr,u,v)){
            cout << "Успех" << endl;
        }
    }


};


