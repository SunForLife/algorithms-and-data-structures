#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

const int INF = 1e9, MAX_SIZE = 1e5;

class Graph {
    struct Edge {
        int to, f, c, p, rev, id;
    };

    int n, s, t;
    std::vector <std::vector <Edge> > g;
    std::vector <int> p;
    std::vector <std::pair <int, int> > pfrom;
    std::vector <std::vector <int> > dec;

    inline void ford_bellman() {
        for (size_t step = 0; step < n; ++step) {
            for (size_t from = 1; from <= n; ++from) {
                for (size_t i = 0; i < g[from].size(); ++i){
                    int to = g[from][i].to, ep = g[from][i].p;
                    if (g[from][i].c - g[from][i].f > 0 && p[to] > p[from] + ep) {
                        p[to] = p[from] + ep;
                        pfrom[to] = {from, i};
                    }
                }
            }
        }
    }

    std::vector <int> d;
    std::set <std::pair <int, int> > ms;
    inline void dij(int x) {
        if (!ms.empty())
            ms.erase(ms.begin());

        for (size_t i = 0; i < g[x].size(); ++i) {
            auto e = g[x][i];
            int to = e.to, ep = e.p;
            if (e.f < e.c && d[to] > p[x] + d[x] + ep - p[to]) {
                ms.erase({d[to], to});
                d[to] = p[x] + d[x] + ep - p[to];
                pfrom[to] = {x, i};
                ms.insert({d[to], to});
            }
        }

        if (!ms.empty())
            dij(ms.begin()->second);
    }

    std::vector <int> used;
    inline bool dfs(int x, int step) {
        used[x] = step + 1;

        if (x == t)
            return true;

        for (auto& e : g[x]) {
            int to = e.to;
            if (e.f > 0 && used[to] <= step) {
                if (dfs(to, step)) {
                    e.f--;
                    dec[step].push_back(e.id);
                    return true;
                }
            }
        }

        return false;
    }

public:

    Graph() {};
    Graph(int _n, int _s, int _t) {
        n = _n, s = _s, t = _t;
        g.resize(n + 1);
    }

    inline void make_edge(int from, int to, int p, int id) {
        g[from].push_back({to, 0, 1, p, g[to].size(), id});
        g[to].push_back({from, 0, 0, -p, g[from].size() - 1, id});
    }

    inline int build_flow(int k) {
        int flow = 0, cnt = 0;
        p.clear(), pfrom.clear();
        p.resize(n + 1, INF);
        p[s] = 0;
        pfrom.resize(n + 1);

        ford_bellman();

        while (p[t] < INF && cnt < k) {
            flow += p[t];
            cnt++;

            int it = t;
            while (it != s) {
                auto& e = g[pfrom[it].first][pfrom[it].second];
                e.f++;
                g[e.to][e.rev].f--;
                it = pfrom[it].first;
            }

            d.clear();
            d.resize(n + 1, 2 * INF);
            pfrom.clear();
            pfrom.resize(n + 1);
            d[s] = 0;
            dij(s);
            for (int i = 1; i <= n; ++i)
                p[i] += d[i];
        }

        return cnt == k ? flow : -1;
    }

    inline void clear_roads() {
        std::vector <int> cnt(MAX_SIZE + 1);
        for (int i = 1; i <= n; ++i) {
            for (auto& e : g[i]) {
                if (e.f == 1)
                    cnt[e.id]++;
            }
        }
        for (int i = 1; i <= n; ++i) {
            for (auto& e : g[i]) {
                if (e.f == 1 && cnt[e.id] == 2)
                    e.f--;
            }
        }
    }

    inline std::vector <std::vector <int> > decomposition(int k) {
        dec.resize(k);
        used.resize(n + 1);
        for (int i = 0; i < k; ++i)
            dfs(s, i);
        return dec;
    }
};

int main() {
    int n, m, k;
    std::cin >> n >> m >> k;

    Graph g(n, 1, n);

    for (size_t i = 0; i < m; ++i) {
        int from, to, p;
        std::cin >> from >> to >> p;

        g.make_edge(from, to, p, i + 1);
        g.make_edge(to, from, p, i + 1);
    }

    double flow = g.build_flow(k);
    if (flow == -1)
        return std::cout << -1, 0;
    else
        std::cout << flow / k << '\n';

    g.clear_roads();
    auto ans = g.decomposition(k);
    for (size_t i = 0; i < k; ++i) {
        std::cout << ans[i].size() << ' ';
        reverse(ans[i].begin(), ans[i].end());
        for (int x : ans[i])
            std::cout << x << ' ';
        std::cout << '\n';
    }
}
