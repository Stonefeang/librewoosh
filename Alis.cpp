vi alis(vi per) {
	auto invert=[](const vi& P, int n, int def) {
		vi ret(n, def);
		for (int i = 0; i < (int)P.size(); i++) if (P[i] != def) ret[P[i]] = i;
		return ret;
	};
	auto split=[](const vi& P, int k, vi& lo, vi& hi, vi& loInd, vi& hiInd) {
		int i = 0;
		for (int e : P) {
			if (e < k) lo.push_back(e), loInd.push_back(i++);
			else hi.push_back(e-k), hiInd.push_back(i++);
		}
	};
	auto expand=[](const vi& P, vi& ind1, vi& ind2, int n, int def) {
		vi ret(n, def);
		for (int i = 0; i < (int)P.size(); i++) if(P[i] != def) ret[ind1[i]] = ind2[P[i]];
		return ret;
	};
	function<vi(const vi&, const vi&)> comb=[&](const vi& P, const vi& invQ) {
		int n = P.size();
		if (n < 100) {
			vi ret = invert(P, n, -1);
			for (int i = 0; i < (int)invQ.size(); i++) {
				int from = invQ[i];
				for (int j = 0; j < i; j++) from += invQ[j] > invQ[i];
				for (int j = from; j > i; j--)
					if (ret[j-1] < ret[j]) swap(ret[j-1], ret[j]);
			}
			return invert(ret, n, -1);
		}
		vi p1, p2, q1, q2, i1, i2, j1, j2;
		split(P, n/2, p1, p2, i1, i2);
		split(invQ, n/2, q1, q2, j1, j2);
		p1 = expand(comb(p1, q1), i1, j1, n, -1);
		p2 = expand(comb(p2, q2), i2, j2, n, n);
		q1 = invert(p1, n, -1);
		q2 = invert(p2, n, n);
		vi ans(n, -1);
		int delta = 0, j = n;
		for (int i = 0; i < n; i++) {
			ans[i] = (p1[i] < 0 ? p2[i] : p1[i]);
			while (j > 0 && delta >= 0)
				delta -= (q2[--j] < i || q1[j] >= i);
			if (p2[i] < j || p1[i] >= j)
				if (delta++ < 0)
					if (q2[j] < i || q1[j] >= i)
						ans[i] = j;
		}
		return ans;
	};
	auto padPerm=[](const vi& P, vi& has, vi& pad, vi& ind, int n) {
		vector<bool> seen(n);
		for (int i = 0; i < (int)P.size(); i++) if(P[i] != -1) {
			ind.push_back(i);
			has.push_back(P[i]);
			seen[P[i]] = 1;
		}
		for (int i = 0; i < n; i++) if (!seen[i]) pad.push_back(i);
	};
	auto mongeMul=[&](const vi& P, const vi& Q, int n) {
		vi h1, p1, i1, h2, p2, i2;
		padPerm(P, h1, p1, i1, Q.size());
		padPerm(invert(Q, n, -1), h2, p2, i2, Q.size());
		h1.insert(h1.begin(), p1.begin(), p1.end());
		h2.insert(h2.end(), p2.begin(), p2.end());
		vi ans(P.size(), -1), tmp = comb(h1, h2);
		for (int i = 0; i < (int)i1.size(); i++) {
			int j = tmp[i+(int)p1.size()];
			if (j < (int)i2.size()) {
				ans[i1[i]] = i2[j];
			}
		}
		return ans;
	};
	function<vi(const vi&)> build=[&](const vi& seq) {
		int n = seq.size();
		if (!n) return vi{};
		int lo = *min_element(seq.begin(), seq.end());
		int hi = *max_element(seq.begin(), seq.end());
		if (lo == hi) {
			vi tmp(n);
			iota(tmp.begin(), tmp.end(), 1);
			tmp.back() = -1;
			return tmp;
		}
		int mid = (lo+hi+1)/2;
		vi p1, p2, i1, i2;
		split(seq, mid, p1, p2, i1, i2);
		p1 = expand(build(p1), i1, i1, n, -1);
		p2 = expand(build(p2), i2, i2, n, -1);
		for (int j : i1) p2[j] = j;
		for (int j : i2) p1[j] = j;
		return mongeMul(p1, p2, n);
	};
	vi ret(per.size(), -1);
	vi take = build(per);
	for (int i = 0; i < (int)per.size(); i++)
		if (take[i] != -1)
			ret[take[i]] = i;
	return ret;
}
