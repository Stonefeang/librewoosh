struct Matching {//numeracja od 1, skojarzenie znajduje sie w wektorze nxt
	int n; 
	vector<int> nxt, fa, S, Q, pre, Vis, head;
	int *Top, Tim = 0;
	struct edge {int to; int next;};
	vector<edge> e;
	Matching(int _n): n(_n), nxt(n + 1), fa(n + 1), S(n + 1), Q(n + 1), pre(n + 1), Vis(n + 1), head(n + 1, -1) {
		Top = Q.data();
	}
	inline void add_edge(int x, int y) {
		e.push_back(edge{y, head[x]});
		head[x] = e.size() - 1;
		e.push_back(edge{x, head[y]});
		head[y] = e.size() - 1;
	}
	int getfa(int x) {return (x == fa[x]) ? x : (fa[x] = getfa(fa[x]));}
	inline int lca(int x, int y) {
		for(++Tim, x = getfa(x), y = getfa(y); ; x^= y^= x^= y)
			if(x) {
				if(Vis[x] == Tim) return x;
				Vis[x] = Tim;
				x = getfa(pre[nxt[x]]);
			}
	}
	inline void blossom(int x, int y, int l) {
		while(getfa(x) != l) {
			pre[x] = y;
			if(S[nxt[x]] == 1) S[*Top = nxt[x]] = 0, *Top++;
			if(fa[x] == x) fa[x] = l;
			if(fa[nxt[x]] == nxt[x]) fa[nxt[x]] = l;
			y = nxt[x];
			x = pre[y];
		}
	}
	int match(int x) {
		for(int i = 1; i <= n; ++i) fa[i] = i;
		fill(S.begin(), S.end(), -1);
		S[*(Top = Q.data()) = x] = 0; ++Top;
		for(int *i = Q.data(); i != Top; *i++)
			for(int T = head[*i]; T >= 0; T = e[T].next) {
				int y = e[T].to;
				if(S[y] == -1) {
					pre[y] = *i; S[y] = 1;
					if(!nxt[y]) {
						for(int u = y, v = *i, lst; v; u = lst, v = pre[u])
							lst = nxt[v], nxt[v] = u, nxt[u] = v;
						return 1;
					}
					S[*Top = nxt[y]] = 0; *Top++;
				}
				else if(!S[y] && getfa(y) != getfa(*i)) {
					int l = lca(y, *i);
					blossom(y, *i, l);
					blossom(*i, y, l);
				}
			}
		return 0;
	}
	int matching_size() {
		int ret = 0;
		for (int i = 1; i <= n; i++)
			if (!nxt[i])
				ret += match(i);
		return ret;
	}
};
