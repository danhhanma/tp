#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef pair<int,int> pii;
const ll INF = 1e18;
const int N = 200005;

/* ==================================================
   1. DSU (co gop theo size) + Kruskal
   ================================================== */
struct DSU {
    vector<int> p, sz;
    DSU(int n = 0) { init(n); }
    void init(int n) { p.resize(n + 1); sz.assign(n + 1, 1); iota(p.begin(), p.end(), 0); }
    int find(int x) { return p[x] == x ? x : p[x] = find(p[x]); }
    bool unite(int a, int b) {            // true neu hai nhom khac nhau va da gop
        a = find(a); b = find(b);
        if (a == b) return false;
        if (sz[a] < sz[b]) swap(a, b);
        p[b] = a; sz[a] += sz[b];
        return true;
    }
    int size(int x) { return sz[find(x)]; }
};

// Kruskal: edges = {w, u, v}, tra ve tong trong so cay khung nho nhat (-1 neu khong lien thong)
ll kruskal(int n, vector<array<ll,3>>& edges) {
    sort(edges.begin(), edges.end());
    DSU d(n); ll total = 0; int cnt = 0;
    for (auto& [w, u, v] : edges)
        if (d.unite(u, v)) total += w, cnt++;
    return cnt == n - 1 ? total : -1;
}

/* ==================================================
   2. BIT (Fenwick): cap nhat diem, tong tien to (1-indexed)
   ================================================== */
struct BIT {
    int n; vector<ll> t;
    BIT(int n = 0) : n(n), t(n + 1, 0) {}
    void add(int i, ll v) { for (; i <= n; i += i & -i) t[i] += v; }
    ll sum(int i) { ll r = 0; for (; i > 0; i -= i & -i) r += t[i]; return r; }
    ll sum(int l, int r) { return sum(r) - sum(l - 1); }
};
// Dem nghich the: duyet i, cong (so phan tu da them > a[i]) = i-1 - bit.sum(a[i]), roi bit.add(a[i], 1) (nen nen toa do truoc)

/* ==================================================
   3. SEGMENT TREE (cap nhat diem, truy van doan, 0-indexed)
   Doi combine va gia tri trung tinh cho min/max/gcd...
   ================================================== */
struct SegTree {
    int n; vector<ll> t;
    static ll combine(ll a, ll b) { return a + b; }   // doi thanh min(a,b) / max(a,b)
    static const ll NEUTRAL = 0;                      // min: INF, max: -INF
    SegTree(int n = 0) : n(n), t(2 * n, NEUTRAL) {}
    void build(vector<ll>& a) {
        for (int i = 0; i < n; i++) t[n + i] = a[i];
        for (int i = n - 1; i >= 1; i--) t[i] = combine(t[2 * i], t[2 * i + 1]);
    }
    void update(int i, ll v) {
        for (t[i += n] = v; i > 1; i >>= 1) t[i >> 1] = combine(t[i & ~1], t[i | 1]);
    }
    ll query(int l, int r) {                          // [l, r] bao gom hai dau
        ll rl = NEUTRAL, rr = NEUTRAL;
        for (l += n, r += n + 1; l < r; l >>= 1, r >>= 1) {
            if (l & 1) rl = combine(rl, t[l++]);
            if (r & 1) rr = combine(t[--r], rr);
        }
        return combine(rl, rr);
    }
};

/* ==================================================
   4. SPARSE TABLE (min/max doan, khong cap nhat, O(1) truy van)
   ================================================== */
struct Sparse {
    vector<vector<int>> st; vector<int> lg;
    Sparse(vector<int>& a) {
        int n = a.size(), K = 1;
        while ((1 << K) <= n) K++;
        st.assign(K, vector<int>(n));
        st[0] = a;
        for (int j = 1; j < K; j++)
            for (int i = 0; i + (1 << j) <= n; i++)
                st[j][i] = min(st[j-1][i], st[j-1][i + (1 << (j-1))]);   // doi thanh max neu can
        lg.assign(n + 2, 0);
        for (int i = 2; i <= n + 1; i++) lg[i] = lg[i / 2] + 1;
    }
    int query(int l, int r) {                         // [l, r]
        int k = lg[r - l + 1];
        return min(st[k][l], st[k][r - (1 << k) + 1]);
    }
};

/* ==================================================
   5. STACK DON DIEU
   ================================================== */
// Phan tu lon hon dau tien ben phai (chi so), -1 neu khong co
vector<int> nextGreater(vector<ll>& a) {
    int n = a.size(); vector<int> res(n, -1); stack<int> st;
    for (int i = 0; i < n; i++) {
        while (!st.empty() && a[st.top()] < a[i]) res[st.top()] = i, st.pop();
        st.push(i);
    }
    return res;
}
// Hinh chu nhat lon nhat trong histogram
ll largestRectangle(vector<ll>& h) {
    int n = h.size(); ll best = 0; stack<int> st;
    for (int i = 0; i <= n; i++) {
        ll cur = (i == n ? 0 : h[i]);
        while (!st.empty() && h[st.top()] >= cur) {
            ll height = h[st.top()]; st.pop();
            int left = st.empty() ? -1 : st.top();
            best = max(best, height * (i - left - 1));
        }
        st.push(i);
    }
    return best;
}

/* ==================================================
   6. MANG HIEU, PREFIX SUM 2D, NEN TOA DO
   ================================================== */
// Cong v len doan [l, r] (0-indexed): diff[l] += v, diff[r+1] -= v, cuoi cung lay prefix sum
// vector<ll> diff(n + 1); diff[l] += v; diff[r+1] -= v; for i: diff[i] += diff[i-1];

// Prefix sum 2D (1-indexed): pre[i][j] = a[i][j] + pre[i-1][j] + pre[i][j-1] - pre[i-1][j-1]
// Tong hinh chu nhat (x1,y1)-(x2,y2) = pre[x2][y2] - pre[x1-1][y2] - pre[x2][y1-1] + pre[x1-1][y1-1]

// Nen toa do: chuyen gia tri lon ve 0..k-1 giu nguyen thu tu
vector<int> compress(vector<ll>& a) {
    vector<ll> s = a;
    sort(s.begin(), s.end()); s.erase(unique(s.begin(), s.end()), s.end());
    vector<int> r(a.size());
    for (size_t i = 0; i < a.size(); i++) r[i] = lower_bound(s.begin(), s.end(), a[i]) - s.begin();
    return r;
}

/* ==================================================
   7. DO THI: TOPO SORT, FLOYD, HAI PHIA, CAY
   ================================================== */
int n, m;
vector<int> adj[N];

// Topo sort Kahn: tra ve thu tu, rong neu co chu trinh
vector<int> topoSort() {
    vector<int> indeg(n + 1, 0), order;
    for (int u = 1; u <= n; u++) for (int v : adj[u]) indeg[v]++;
    queue<int> q;
    for (int u = 1; u <= n; u++) if (!indeg[u]) q.push(u);
    while (!q.empty()) {
        int u = q.front(); q.pop(); order.push_back(u);
        for (int v : adj[u]) if (--indeg[v] == 0) q.push(v);
    }
    if ((int)order.size() < n) order.clear();         // co chu trinh
    return order;
}
// DP tren DAG: duyet theo thu tu topo, dp[v] = max(dp[v], dp[u] + w)

// Floyd-Warshall (n <= 400). d[i][i] = 0, khong co canh = INF
// Nho kiem tra d[i][k] < INF && d[k][j] < INF khi co INF lon de tranh tran
void floyd(vector<vector<ll>>& d, int n) {
    for (int k = 1; k <= n; k++)
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                if (d[i][k] < INF && d[k][j] < INF)
                    d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
}

// Kiem tra do thi hai phia (to mau BFS)
bool isBipartite() {
    vector<int> color(n + 1, -1);
    for (int s = 1; s <= n; s++) {
        if (color[s] != -1) continue;
        queue<int> q; q.push(s); color[s] = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : adj[u]) {
                if (color[v] == -1) color[v] = color[u] ^ 1, q.push(v);
                else if (color[v] == color[u]) return false;
            }
        }
    }
    return true;
}

// Cay: kich thuoc cay con (goc 1), dung DFS de quy (cay sau ~1e5 co the can duyet lap)
int sub[N];
void calcSubtree(int u, int p) {
    sub[u] = 1;
    for (int v : adj[u]) if (v != p) { calcSubtree(v, u); sub[u] += sub[v]; }
}
// Duong kinh cay: BFS tu 1 tim dinh xa nhat A, BFS tu A tim dinh xa nhat B, dist(A,B) la duong kinh

/* ==================================================
   8. DP KINH DIEN
   ================================================== */
// Knapsack 0/1: moi mon dung toi da 1 lan, W la suc chua
ll knapsack01(vector<int>& wt, vector<ll>& val, int W) {
    vector<ll> dp(W + 1, 0);
    for (size_t i = 0; i < wt.size(); i++)
        for (int w = W; w >= wt[i]; w--)               // chay NGUOC
            dp[w] = max(dp[w], dp[w - wt[i]] + val[i]);
    return dp[W];
}
// Knapsack khong gioi han so lan: doi vong w thanh chay XUOI (w tu wt[i] len W)

// Doi tien: so dong xu it nhat de tao tong S (-1 neu khong duoc)
int coinChange(vector<int>& coins, int S) {
    vector<int> dp(S + 1, INT_MAX); dp[0] = 0;
    for (int s = 1; s <= S; s++)
        for (int c : coins)
            if (s >= c && dp[s - c] != INT_MAX) dp[s] = min(dp[s], dp[s - c] + 1);
    return dp[S] == INT_MAX ? -1 : dp[S];
}

// LIS O(n log n): do dai day con tang nghiem ngat dai nhat
int LIS(vector<ll>& a) {
    vector<ll> tail;
    for (ll x : a) {
        auto it = lower_bound(tail.begin(), tail.end(), x);   // khong nghiem ngat: dung upper_bound
        if (it == tail.end()) tail.push_back(x); else *it = x;
    }
    return tail.size();
}

// LCS hai chuoi: O(n * m)
int LCS(const string& s, const string& t) {
    int a = s.size(), b = t.size();
    vector<vector<int>> dp(a + 1, vector<int>(b + 1, 0));
    for (int i = 1; i <= a; i++)
        for (int j = 1; j <= b; j++)
            dp[i][j] = s[i-1] == t[j-1] ? dp[i-1][j-1] + 1 : max(dp[i-1][j], dp[i][j-1]);
    return dp[a][b];
}

// Khoang cach chinh sua (edit distance)
int editDistance(const string& s, const string& t) {
    int a = s.size(), b = t.size();
    vector<vector<int>> dp(a + 1, vector<int>(b + 1, 0));
    for (int i = 0; i <= a; i++) dp[i][0] = i;
    for (int j = 0; j <= b; j++) dp[0][j] = j;
    for (int i = 1; i <= a; i++)
        for (int j = 1; j <= b; j++)
            dp[i][j] = min({dp[i-1][j] + 1, dp[i][j-1] + 1, dp[i-1][j-1] + (s[i-1] != t[j-1])});
    return dp[a][b];
}

// Duong di tren luoi (chi phai / xuong), dem so cach, bo qua o '#'
// dp[i][j] = (g[i][j]=='#') ? 0 : dp[i-1][j] + dp[i][j-1], dp[0][0] = 1

/* ==================================================
   9. XAU: KMP (ham tien to), Z-function
   ================================================== */
vector<int> prefixFunction(const string& s) {        // pi[i] = do dai tien to = hau to dai nhat cua s[0..i]
    int n = s.size(); vector<int> pi(n, 0);
    for (int i = 1; i < n; i++) {
        int j = pi[i - 1];
        while (j > 0 && s[i] != s[j]) j = pi[j - 1];
        if (s[i] == s[j]) j++;
        pi[i] = j;
    }
    return pi;
}
// Tim mau p trong van ban t: tinh prefixFunction(p + '#' + t), vi tri i co pi[i] == |p| la mot lan khop
// Chu ky nho nhat cua s: k = n - pi[n-1], neu n % k == 0 thi chu ky la k, nguoc lai la n

vector<int> zFunction(const string& s) {             // z[i] = do dai tien to chung dai nhat cua s va s[i..]
    int n = s.size(); vector<int> z(n, 0);
    for (int i = 1, l = 0, r = 0; i < n; i++) {
        if (i < r) z[i] = min(r - i, z[i - l]);
        while (i + z[i] < n && s[z[i]] == s[i + z[i]]) z[i]++;
        if (i + z[i] > r) l = i, r = i + z[i];
    }
    return z;
}

// Hash chuoi (mod 2^61 hoac 2 mod): neu can, dung mod 1e9+7 va 998244353 cung luc de tranh va cham

/* ==================================================
   10. CHAT NHI PHAN, BITMASK
   ================================================== */
// Chat nhi phan tren dap an: tim x nho nhat thoa check(x), check don dieu (sai...sai dung...dung)
// ll lo = 0, hi = 1e18, ans = -1;
// while (lo <= hi) { ll mid = lo + (hi - lo) / 2; if (check(mid)) ans = mid, hi = mid - 1; else lo = mid + 1; }
// Chat nhi phan so thuc: lap 100 lan  mid = (lo + hi) / 2

// Duyet moi tap con
// for (int mask = 0; mask < (1 << n); mask++)
//     for (int i = 0; i < n; i++) if (mask >> i & 1) { /* phan tu i co trong tap */ }
// Duyet tap con cua mask: for (int sub = mask; sub > 0; sub = (sub - 1) & mask)
// __builtin_popcount(x) dem bit 1, __builtin_ctz(x) so bit 0 cuoi, __lg(x) vi tri bit cao nhat
// Dung ll thi dung 1LL << i va __builtin_popcountll

// Quay lui (n nho): sinh hoan vi bang next_permutation
// sort(a.begin(), a.end()); do { /* xu ly a */ } while (next_permutation(a.begin(), a.end()));

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    return 0;
}