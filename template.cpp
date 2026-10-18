#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef pair<int,int> pii;
typedef pair<ll,int> pli;
const ll INF = 1e18;
const int N = 200005;

/* ==================================================
   1. DO THI: DFS / BFS (danh sach ke, dinh 1..n)
   ================================================== */
int n, m;
vector<int> adj[N];
bool vis[N];

// DFS de quy (n lon ~1e5 tren duong thang co the tran stack -> dung ban lap)
void dfs(int u, int p = 0) {
    vis[u] = true;
    for (int v : adj[u]) if (v != p && !vis[v]) dfs(v, u);
}

// DFS lap bang stack
void dfs_iter(int s) {
    stack<int> st; st.push(s); vis[s] = true;
    while (!st.empty()) {
        int u = st.top(); st.pop();
        for (int v : adj[u]) if (!vis[v]) { vis[v] = true; st.push(v); }
    }
}

// BFS: khoang cach ngan nhat tren do thi khong trong so
int d[N];
void bfs(int s) {
    fill(d, d + n + 1, -1);
    queue<int> q; q.push(s); d[s] = 0;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : adj[u]) if (d[v] == -1) { d[v] = d[u] + 1; q.push(v); }
    }
}

// Dem so thanh phan lien thong
int components() {
    int c = 0;
    fill(vis, vis + n + 1, false);
    for (int i = 1; i <= n; i++) if (!vis[i]) { c++; dfs_iter(i); }
    return c;
}

/* ==================================================
   2. LUOI / MA TRAN: BFS 4 huong (0-indexed)
   '#' la tuong, '.' la duong di
   ================================================== */
int R, C;
vector<string> g;
int dx[4] = {1, -1, 0, 0};
int dy[4] = {0, 0, 1, -1};
// 8 huong: dx[8]={1,1,1,0,0,-1,-1,-1}, dy[8]={1,0,-1,1,-1,1,0,-1}

bool inside(int x, int y) { return x >= 0 && x < R && y >= 0 && y < C; }

// BFS tu (sx,sy), tra ve mang khoang cach (-1 neu khong toi duoc)
vector<vector<int>> gridBFS(int sx, int sy) {
    vector<vector<int>> dist(R, vector<int>(C, -1));
    queue<pii> q; q.push({sx, sy}); dist[sx][sy] = 0;
    while (!q.empty()) {
        auto [x, y] = q.front(); q.pop();
        for (int k = 0; k < 4; k++) {
            int nx = x + dx[k], ny = y + dy[k];
            if (!inside(nx, ny) || g[nx][ny] == '#' || dist[nx][ny] != -1) continue;
            dist[nx][ny] = dist[x][y] + 1;
            q.push({nx, ny});
        }
    }
    return dist;
}
// Multi-source BFS: push TAT CA nguon vao queue voi dist = 0 truoc khi chay.

// Dem so vung lien thong tren luoi (flood fill)
void flood(int sx, int sy, vector<vector<bool>>& seen) {
    queue<pii> q; q.push({sx, sy}); seen[sx][sy] = true;
    while (!q.empty()) {
        auto [x, y] = q.front(); q.pop();
        for (int k = 0; k < 4; k++) {
            int nx = x + dx[k], ny = y + dy[k];
            if (inside(nx, ny) && g[nx][ny] != '#' && !seen[nx][ny]) {
                seen[nx][ny] = true; q.push({nx, ny});
            }
        }
    }
}

// 0-1 BFS: canh trong so 0 hoac 1 (deque, nhanh hon Dijkstra)
vector<int> bfs01(int s, vector<vector<pii>>& G) { // G[u] = {v, w}, w in {0,1}
    vector<int> dd(n + 1, INT_MAX);
    deque<int> dq; dd[s] = 0; dq.push_back(s);
    while (!dq.empty()) {
        int u = dq.front(); dq.pop_front();
        for (auto [v, w] : G[u])
            if (dd[u] + w < dd[v]) {
                dd[v] = dd[u] + w;
                if (w == 0) dq.push_front(v); else dq.push_back(v);
            }
    }
    return dd;
}

/* ==================================================
   3. DIJKSTRA (trong so khong am)
   ================================================== */
vector<pii> G[N]; // G[u] = {v, w}
ll dist_[N];
int par[N]; // truy vet duong di

void dijkstra(int s) {
    fill(dist_, dist_ + n + 1, INF);
    fill(par, par + n + 1, 0);
    priority_queue<pli, vector<pli>, greater<pli>> pq;
    dist_[s] = 0; pq.push({0, s});
    while (!pq.empty()) {
        auto [du, u] = pq.top(); pq.pop();
        if (du > dist_[u]) continue;
        for (auto [v, w] : G[u])
            if (dist_[u] + w < dist_[v]) {
                dist_[v] = dist_[u] + w;
                par[v] = u;
                pq.push({dist_[v], v});
            }
    }
}
vector<int> getPath(int t) { // goi sau dijkstra(s)
    vector<int> path;
    if (dist_[t] == INF) return path;
    for (int v = t; v; v = par[v]) path.push_back(v);
    reverse(path.begin(), path.end());
    return path;
}

// Dijkstra tren luoi: chi phi di vao o (x,y) la cost[x][y]
vector<vector<ll>> gridDijkstra(int sx, int sy, vector<vector<int>>& cost) {
    vector<vector<ll>> D(R, vector<ll>(C, INF));
    typedef tuple<ll,int,int> T;
    priority_queue<T, vector<T>, greater<T>> pq;
    D[sx][sy] = cost[sx][sy]; pq.push({D[sx][sy], sx, sy});
    while (!pq.empty()) {
        auto [du, x, y] = pq.top(); pq.pop();
        if (du > D[x][y]) continue;
        for (int k = 0; k < 4; k++) {
            int nx = x + dx[k], ny = y + dy[k];
            if (!inside(nx, ny) || g[nx][ny] == '#') continue;
            if (du + cost[nx][ny] < D[nx][ny]) {
                D[nx][ny] = du + cost[nx][ny];
                pq.push({D[nx][ny], nx, ny});
            }
        }
    }
    return D;
}

// Dijkstra co trang thai: dist[u][k] (vd: da dung k lan "ve mien phi")
// state = {dist, u, k}, canh thuong: k giu nguyen, canh dac biet: k+1

/* ==================================================
   4. LCA - BINARY LIFTING (cay, goc = 1)
   ================================================== */
const int LOG = 18; // 2^18 > 2e5
vector<int> tr[N];
int up[N][LOG], dep[N];

void buildLCA(int root = 1) {
    // BFS de tranh tran stack
    queue<int> q; q.push(root);
    dep[root] = 0;
    for (int j = 0; j < LOG; j++) up[root][j] = root;
    vector<bool> seen(n + 1, false); seen[root] = true;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : tr[u]) if (!seen[v]) {
            seen[v] = true;
            dep[v] = dep[u] + 1;
            up[v][0] = u;
            for (int j = 1; j < LOG; j++) up[v][j] = up[up[v][j-1]][j-1];
            q.push(v);
        }
    }
}
int lca(int u, int v) {
    if (dep[u] < dep[v]) swap(u, v);
    int diff = dep[u] - dep[v];
    for (int j = 0; j < LOG; j++) if (diff >> j & 1) u = up[u][j];
    if (u == v) return u;
    for (int j = LOG - 1; j >= 0; j--)
        if (up[u][j] != up[v][j]) { u = up[u][j]; v = up[v][j]; }
    return up[u][0];
}
int treeDist(int u, int v) { return dep[u] + dep[v] - 2 * dep[lca(u, v)]; }
int kthAncestor(int u, int k) {
    for (int j = 0; j < LOG; j++) if (k >> j & 1) u = up[u][j];
    return u;
}

/* ==================================================
   5. DIGIT DP
   Dem so x trong [L, R] thoa tinh chat nao do.
   Vi du ben duoi: tong cac chu so == S.
   Doi phan "state" va dieu kien o ham go() cho bai khac.
   ================================================== */
string num;
int S;
int memo[20][2][200];// [vi tri][state=tong][started]

// pos: dang xet chu so thu pos (tu trai), sum: state, tight: con bi chan boi num khong
// started: da gap chu so khac 0 chua (de xu ly so 0 dung dau)
int dfs(int pos, int tight, int sum){
    // if(check(sum) == 0) return 0;
    if(pos == digits.length()){
        if(check(sum) == 1) return 1;
        else return 0;
    }
    if(memo[pos][tight][sum]!=-1) return memo[pos][tight][sum];
    int ans = 0;
    int limit = tight ? (digits[pos]-'0') : 9;
    foi(d,0,limit){
        bool new_tight = tight && (d == limit);
        int new_sum = sum + d;
        ans+= dfs(pos+1,new_tight,new_sum);
    }
    memo[pos][tight][sum] = ans;
    return ans;
}
// Dem x trong [0, X] (X >= 0). Luu y so 0 duoc tinh la 1 so.
ll countUpTo(ll X) {
    if (X < 0) return 0;
    num = to_string(X);
    memset(memo, -1, sizeof(memo));
    return go(0, 0, true, false);
}
// dap an [L, R] = countUpTo(R) - countUpTo(L - 1)

/* Mot so state hay dung cho digit DP:
   - chia het cho K: state = so du (rem), nrem = (rem * 10 + dgt) % K
   - khong co 2 chu so lien tiep bang nhau: state = chu so truoc (prev), 10 = chua co
   - dem so chu so 'd': state = so lan xuat hien
   - cac chu so khong giam: state = chu so truoc, dgt >= prev
   Nho kich thuoc memo theo so state. */

/* ==================================================
   MAIN: vi du doc do thi vo huong
   ================================================== */

void nhiphan(){
    ll lo = 0, hi = 1e18, ans = -1;   // tìm x nhỏ nhất thỏa check
    while (lo <= hi) {
        ll mid = (lo + hi) / 2;
        if (check(mid)) ans = mid, hi = mid - 1;
        else lo = mid + 1;
    }
}

ll power(ll a, ll b, ll mod) {
    ll r = 1; a %= mod;
    for (; b; b >>= 1, a = a * a % mod) if (b & 1) r = r * a % mod;
    return r;
}

void cuphap(){
    vector<int> a;                    // rỗng
    vector<int> a(n);                 // n phần tử = 0
    vector<int> a(n, -1);             // n phần tử = -1
    vector<vector<int>> g(n, vector<int>(m, 0));  // ma trận n x m

    a.push_back(x);                   // thêm cuối
    a.pop_back();                     // xóa cuối
    a.size();  a.empty();  a.clear();
    a.front(); a.back();              // phần tử đầu / cuối
    a.insert(a.begin() + i, x);       // chèn tại i (O(n))
    a.erase(a.begin() + i);           // xóa tại i (O(n))
    a.resize(k);

    sort(a.begin(), a.end());                       // tăng dần
    sort(a.rbegin(), a.rend());                     // giảm dần
    sort(a.begin(), a.end(), greater<int>());       // giảm dần
    reverse(a.begin(), a.end());
    a.erase(unique(a.begin(), a.end()), a.end());   // xóa trùng (phải sort trước)

    *max_element(a.begin(), a.end());               // giá trị max
    min_element(a.begin(), a.end()) - a.begin();    // vị trí min
    accumulate(a.begin(), a.end(), 0LL);            // tổng (0LL để dùng long long)
    lower_bound(a.begin(), a.end(), x) - a.begin(); // vị trí đầu tiên >= x
    upper_bound(a.begin(), a.end(), x) - a.begin(); // vị trí đầu tiên > x
    set<int> st;                      // tự sắp tăng, không trùng
    st.insert(x);
    st.erase(x);                      // xóa theo giá trị
    st.erase(it);                     // xóa theo iterator
    st.count(x);                      // 0 hoặc 1 (dùng để kiểm tra tồn tại)
    st.find(x) != st.end();           // cũng là kiểm tra tồn tại
    st.size();  st.empty();  st.clear();

    *st.begin();                      // nhỏ nhất
    *st.rbegin();                     // lớn nhất
    st.lower_bound(x);                // iterator tới phần tử đầu tiên >= x
    st.upper_bound(x);                // iterator tới phần tử đầu tiên > x
    // KHÔNG dùng lower_bound(st.begin(), st.end(), x) vì chậm O(n)

    for (int x : st) cout << x;       // duyệt theo thứ tự tăng

    set<int, greater<int>> st2;       // sắp giảm dần
    auto it = st.upper_bound(x);
    if (it != st.begin()) { --it; /* *it là kết quả */ }
        map<string,int> mp;               // key tự sắp tăng
    mp[key]++;                        // tự tạo key với giá trị 0 nếu chưa có
    mp[key] = v;
    mp.count(key);                    // kiểm tra tồn tại
    mp.find(key) != mp.end();
    mp.erase(key);
    mp.size();  mp.empty();  mp.clear();

    for (auto &[k, v] : mp) cout << k << " " << v;   // duyệt theo key tăng
    for (auto &p : mp) p.first, p.second;             // cách viết cũ

    mp.begin()->first;                // key nhỏ nhất
    mp.rbegin()->first;               // key lớn nhất
    mp.lower_bound(x);  mp.upper_bound(x);

    map<pair<int,int>, int> mp2;      // key là cặp
    mp2[{1, 2}] = 5;

    pair<int,int> p = {1, 2};  p.first;  p.second;
vector<pair<int,int>> v;  v.push_back({1, 2});  v.emplace_back(1, 2);
sort(v.begin(), v.end());         // sort theo first rồi tới second
sort(v.begin(), v.end(), [](auto &a, auto &b){ return a.second < b.second; });  // comparator tự viết

queue<int> q;  q.push(x);  q.front();  q.pop();
stack<int> st; st.push(x);  st.top();  st.pop();
priority_queue<int> pq;                           // max-heap
priority_queue<int, vector<int>, greater<int>> pq2; // min-heap
deque<int> dq;  dq.push_front(x);  dq.push_back(x);
}


vector<vector<pair<int,int>>> adj(100005);
int visited[100005];
vi dist(100005,inf);


void dijisktra(int s, int n){
    dist[s] = 0;
    priority_queue<pair<int,int> , vector<pair<int,int>>, greater<pair<int,int>>> pq;
    pq.push({0,s});

    while(!pq.empty()){
        int u = pq.top().se;
        int d = pq.top().se;
        pq.pop();

        if(visited[u]) continue;
        visited[u] = 1;

        for(auto edge: adj[u]){
            int v = edge.fi;
            int weight = edge.se;

            if(dist[u] + weight < dist[v]){
                dist[v] = dist[u] + weight;
                pq.push({dist[v],v});
            }
        }
    }
    foi(i,1,n){
        cout<< dist[i] <<  ' ';
    }
}



int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // cin >> n >> m;
    // for (int i = 0; i < m; i++) {
    //     int u, v; cin >> u >> v;
    //     adj[u].push_back(v);
    //     adj[v].push_back(u);
    // }
    // Dijkstra: cin >> u >> v >> w; G[u].push_back({v,w}); G[v].push_back({u,w});
    // Luoi: cin >> R >> C; g.resize(R); for (auto &row : g) cin >> row;

    return 0;
}