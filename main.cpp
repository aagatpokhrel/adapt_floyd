#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace std;

typedef long long ll;
const ll INF = 1e16; // Use a safe INF to avoid overflow

typedef vector<vector<ll>> Matrix;

Matrix multiply(const Matrix &A, const Matrix &B, int N)
{
    Matrix C(N, vector<ll>(N, INF));
    for (int i = 0; i < N; ++i)
    {
        for (int k = 0; k < N; ++k)
        {
            if (A[i][k] == INF)
                continue;
            for (int j = 0; j < N; ++j)
            {
                if (B[k][j] == INF)
                    continue;
                C[i][j] = min(C[i][j], A[i][k] + B[k][j]);
            }
        }
    }
    return C;
}

Matrix power(Matrix A, int p, int N)
{
    Matrix res(N, vector<ll>(N, INF));
    for (int i = 0; i < N; ++i)
        res[i][i] = 0;
    while (p > 0)
    {
        if (p & 1)
            res = multiply(res, A, N);
        A = multiply(A, A, N);
        p >>= 1;
    }
    return res;
}

int main()
{
    int N;
    while (cin >> N)
    {
        string line;
        getline(cin, line); // consume newline

        auto read_vec = [&]()
        {
            getline(cin, line);
            stringstream ss(line);
            int val;
            vector<int> v;
            while (ss >> val)
                v.push_back(val);
            return v;
        };

        vector<int> from = read_vec();
        vector<int> to = read_vec();
        vector<int> weight = read_vec();
        ll charges;
        cin >> charges;

        Matrix S(N, vector<ll>(N, INF));
        Matrix L(N, vector<ll>(N, -1));

        for (int i = 0; i < N; ++i)
            S[i][i] = 0;
        for (size_t i = 0; i < from.size(); ++i)
        {
            int u = from[i] - 1;
            int v = to[i] - 1;
            S[u][v] = min(S[u][v], (ll)weight[i]);
            L[u][v] = max(L[u][v], (ll)weight[i]);
        }

        // Part 1: Floyd-Warshall
        for (int k = 0; k < N; ++k)
        {
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    S[i][j] = min(S[i][j], S[i][k] + S[k][j]);
                }
            }
        }

        if (charges == 0)
        {
            cout << S[0][N - 1] << endl;
            continue;
        }

        // Part 2: One charge matrix A
        Matrix A(N, vector<ll>(N, INF));
        for (int f = 0; f < N; ++f)
        {
            for (int t = 0; t < N; ++t)
            {
                for (int i = 0; i < N; ++i)
                {
                    for (int j = 0; j < N; ++j)
                    {
                        if (L[i][j] != -1 && S[f][i] != INF && S[j][t] != INF)
                        {
                            A[f][t] = min(A[f][t], S[f][i] - L[i][j] + S[j][t]);
                        }
                    }
                }
            }
        }

        // Allow for using 0 charges within the exponentiation step
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                A[i][j] = min(A[i][j], S[i][j]);
            }
        }

        // Part 3: Matrix Power
        Matrix finalA = power(A, charges, N);
        cout << finalA[0][N - 1] << endl;
    }
    return 0;
}