#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
const ll MOD = 1e9 + 7;   // doi thanh 998244353 neu de yeu cau

/* ==================================================
   0. CO BAN: gcd, lcm, luy thua nhanh, nghich dao
   ================================================== */
ll gcdll(ll a, ll b) { return b ? gcdll(b, a % b) : a; }   // hoac __gcd(a,b), std::gcd
ll lcmll(ll a, ll b) { return a / gcdll(a, b) * b; }       // chia truoc nhan sau de do tran

ll power(ll a, ll b, ll mod = MOD) {
    ll r = 1; a %= mod; if (a < 0) a += mod;
    for (; b > 0; b >>= 1, a = a * a % mod)
        if (b & 1) r = r * a % mod;
    return r;
}
// Nghich dao khi mod NGUYEN TO: a^(mod-2)
ll inv(ll a, ll mod = MOD) { return power(a, mod - 2, mod); }

// Euclid mo rong: tim x, y sao cho a*x + b*y = gcd(a,b)
ll extgcd(ll a, ll b, ll &x, ll &y) {
    if (!b) { x = 1; y = 0; return a; }
    ll x1, y1; ll g = extgcd(b, a % b, x1, y1);
    x = y1; y = x1 - (a / b) * y1;
    return g;
}
// Nghich dao khi mod KHONG nguyen to (gcd(a,mod) phai = 1), tra -1 neu khong ton tai
ll invGeneral(ll a, ll mod) {
    ll x, y; ll g = extgcd(((a % mod) + mod) % mod, mod, x, y);
    if (g != 1) return -1;
    return ((x % mod) + mod) % mod;
}

/* ==================================================
   1. KIEM TRA SO NGUYEN TO
   ================================================== */
// n <= 1e12 : thu chia den sqrt, O(sqrt n)
bool isPrime(ll n) {
    if (n < 2) return false;
    if (n < 4) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (ll i = 5; i * i <= n; i += 6)
        if (n % i == 0 || n % (i + 2) == 0) return false;
    return true;
}

// n lon toi 1e18 : Miller-Rabin (xac dinh)
ll mulmod(ll a, ll b, ll m) { return (__int128)a * b % m; }
ll powmod128(ll a, ll b, ll m) {
    ll r = 1; a %= m;
    for (; b; b >>= 1, a = mulmod(a, a, m)) if (b & 1) r = mulmod(r, a, m);
    return r;
}
bool millerRabin(ll n) {
    if (n < 2) return false;
    for (ll p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        if (n % p == 0) return n == p;
    }
    ll d = n - 1; int s = 0;
    while (d % 2 == 0) d /= 2, s++;
    for (ll a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        ll x = powmod128(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int i = 1; i < s; i++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { comp = false; break; }
        }
        if (comp) return false;
    }
    return true;
}

/* ==================================================
   2. SANG NGUYEN TO
   ================================================== */
const int MX = 10000005;       // chinh theo gioi han de (1e7 la toi da hop ly)
bitset<MX> comp;               // comp[i] = true neu i la hop so
vector<int> primes;

// Sang Eratosthenes: O(n log log n)
void sieve(int n) {
    comp[0] = comp[1] = true;
    for (int i = 2; (ll)i * i <= n; i++)
        if (!comp[i])
            for (int j = i * i; j <= n; j += i) comp[j] = true;
    for (int i = 2; i <= n; i++) if (!comp[i]) primes.push_back(i);
}

// Sang tuyen tinh + uoc nguyen to nho nhat (spf): O(n)
// Dung khi can phan tich NHIEU so <= 1e6 ... 1e7
int spf[MX];                   // spf[x] = uoc nguyen to nho nhat cua x
void linearSieve(int n) {
    primes.clear();
    for (int i = 2; i <= n; i++) {
        if (!spf[i]) { spf[i] = i; primes.push_back(i); }
        for (int p : primes) {
            if (p > spf[i] || (ll)i * p > n) break;
            spf[i * p] = p;
        }
    }
}

// Sang doan: tim so nguyen to trong [L, R] voi R <= 1e12, R - L <= 1e6
// Can goi sieve(sqrt(R)) truoc
vector<ll> segmentedSieve(ll L, ll R) {
    vector<bool> isp(R - L + 1, true);
    for (ll p : primes) {
        if (p * p > R) break;
        for (ll j = max(p * p, (L + p - 1) / p * p); j <= R; j += p) isp[j - L] = false;
    }
    vector<ll> res;
    for (ll i = max(L, 2LL); i <= R; i++) if (isp[i - L]) res.push_back(i);
    return res;
}

/* ==================================================
   3. PHAN TICH THUA SO NGUYEN TO
   ================================================== */
// (a) Thu chia, 1 so n <= 1e12 : O(sqrt n)
vector<pair<ll,int>> factorize(ll n) {
    vector<pair<ll,int>> f;                // {nguyen to, so mu}
    for (ll p = 2; p * p <= n; p++) {
        if (n % p == 0) {
            int c = 0;
            while (n % p == 0) n /= p, c++;
            f.push_back({p, c});
        }
    }
    if (n > 1) f.push_back({n, 1});
    return f;
}

// (b) Dung spf, nhieu so <= MX : O(log n) moi so (goi linearSieve truoc)
vector<pair<int,int>> factorizeSPF(int n) {
    vector<pair<int,int>> f;
    while (n > 1) {
        int p = spf[n], c = 0;
        while (n % p == 0) n /= p, c++;
        f.push_back({p, c});
    }
    return f;
}

// Cac ham suy ra tu phan tich
ll numDivisors(ll n) {                    // so uoc
    ll r = 1;
    for (auto [p, c] : factorize(n)) r *= (c + 1);
    return r;
}
ll sumDivisors(ll n) {                    // tong uoc (khong mod)
    ll r = 1;
    for (auto [p, c] : factorize(n)) {
        ll s = 1, pw = 1;
        for (int i = 0; i < c; i++) pw *= p, s += pw;
        r *= s;
    }
    return r;
}
ll phi(ll n) {                            // phi Euler: so so <= n nguyen to cung nhau voi n
    ll r = n;
    for (auto [p, c] : factorize(n)) r = r / p * (p - 1);
    return r;
}
vector<ll> allDivisors(ll n) {            // liet ke moi uoc, O(sqrt n)
    vector<ll> d;
    for (ll i = 1; i * i <= n; i++)
        if (n % i == 0) { d.push_back(i); if (i != n / i) d.push_back(n / i); }
    sort(d.begin(), d.end());
    return d;
}

// So mu cua nguyen to p trong n! (cong thuc Legendre), vd dem so chu so 0 cuoi cua n! la min(v2, v5) = v5
ll legendre(ll n, ll p) {
    ll r = 0;
    while (n) n /= p, r += n;
    return r;
}

// Dem so uoc cua MOI so 1..n (O(n log n)); tuong tu co the tinh tong uoc
vector<int> divisorCountAll(int n) {
    vector<int> cnt(n + 1, 0);
    for (int i = 1; i <= n; i++)
        for (int j = i; j <= n; j += i) cnt[j]++;
    return cnt;
}

// Phi cho moi so 1..n
vector<int> phiAll(int n) {
    vector<int> ph(n + 1);
    iota(ph.begin(), ph.end(), 0);
    for (int i = 2; i <= n; i++)
        if (ph[i] == i)                   // i la nguyen to
            for (int j = i; j <= n; j += i) ph[j] -= ph[j] / i;
    return ph;
}

/* ==================================================
   4. TO HOP
   ================================================== */
const int NC = 1000005;                   // chinh theo n toi da
ll fact[NC], ifact[NC];

// Tien xu ly giai thua va nghich dao giai thua (mod NGUYEN TO)
void initComb(int n = NC - 1) {
    fact[0] = 1;
    for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * i % MOD;
    ifact[n] = inv(fact[n]);
    for (int i = n; i >= 1; i--) ifact[i - 1] = ifact[i] * i % MOD;
}
ll C(int n, int k) {                      // C(n, k) mod MOD
    if (k < 0 || k > n) return 0;
    return fact[n] * ifact[k] % MOD * ifact[n - k] % MOD;
}
ll P(int n, int k) {                      // chinh hop A(n, k)
    if (k < 0 || k > n) return 0;
    return fact[n] * ifact[n - k] % MOD;
}

// n nho (<= 2000) : tam giac Pascal, dung duoc voi MOD bat ky, khong can nguyen to
ll pas[2005][2005];
void initPascal(int n, ll mod = MOD) {
    for (int i = 0; i <= n; i++) {
        pas[i][0] = 1;
        for (int j = 1; j <= i; j++)
            pas[i][j] = (pas[i - 1][j - 1] + (j <= i - 1 ? pas[i - 1][j] : 0)) % mod;
    }
}

// Khong mod, ket qua vua long long (k nho, n <= 60 an toan): nhan roi chia lien tuc
ll Cexact(ll n, ll k) {
    k = min(k, n - k);
    ll r = 1;
    for (ll i = 1; i <= k; i++) r = r * (n - k + i) / i;   // luon chia het
    return r;
}

// n rat lon (toi 1e18) nhung mod p NHO NGUYEN TO (p <= 1e6): dinh ly Lucas
// Can fact/ifact tinh theo mod p (viet lai voi p thay MOD). Ban gon:
ll lucas(ll n, ll k, ll p, vector<ll>& f, vector<ll>& fi) {  // f, fi la giai thua, nghich dao theo mod p
    if (k < 0 || k > n) return 0;
    ll res = 1;
    while (n || k) {
        ll a = n % p, b = k % p;
        if (b > a) return 0;
        res = res * f[a] % p * fi[b] % p * fi[a - b] % p;
        n /= p; k /= p;
    }
    return res;
}

// Cac cong thuc dem hay gap (dung C o tren)
ll catalan(int n) { return C(2 * n, n) * inv(n + 1) % MOD; }              // so day ngoac hop le, ...
ll starsAndBars(int n, int k) { return C(n + k - 1, k - 1); }             // chia n vat giong nhau cho k nguoi (moi nguoi >= 0)
// Chia n vat cho k nguoi, moi nguoi >= 1: C(n - 1, k - 1)
// So duong di tren luoi (0,0)->(a,b) chi di phai/xuong: C(a + b, a)
// So hoan vi co lap: n! / (c1! * c2! * ...)
// Bao ham - loai tru: |A1 U A2 U ...| = tong |Ai| - tong |Ai ^ Aj| + ... (dung bitmask khi so tap nho)

/* ==================================================
   5. MOT SO TRICK HAY GAP
   ================================================== */
// Tong 1..n = n*(n+1)/2   (chia truoc khi nhan khi mod: dung inv(2) hoac chia chan truoc)
// Tong cap so cong: (a1 + an) * n / 2
// Do dong du: (a - b) % m co the am -> ((a - b) % m + m) % m
// So uoc lon nhat co the cua n <= 1e9 la 1344, n <= 1e18 la ~103680
// Moi so <= 1e9 co toi da 9 uoc nguyen to khac nhau
// a^b mod m voi b rat lon (b la chuoi): dung b mod (m-1) khi m nguyen to va gcd(a,m)=1

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // vi du su dung:
    // sieve(1000000);                       // roi dung comp[x], primes
    // linearSieve(1000000);                 // roi factorizeSPF(x)
    // initComb(); cout << C(10, 3) << "\n"; // 120
    // for (auto [p, c] : factorize(360)) cout << p << "^" << c << "\n";

    return 0;
}