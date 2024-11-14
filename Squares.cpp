//wymaga przepisania funkcji Pref
vector<pair<pii, int>> squares(const char* s, int n) {
	vector<pair<pii, int>> ans;
	vector pos(n + 2, -1);
	for (int mid=1; mid<n; mid++) {
		int part = mid & ~(mid - 1), off = mid - part;
		int end = min(mid + part, n);
		string a(s + off, s + off + part),
			b(s + mid, s + end),
			ra(a.rbegin(), a.rend());
		for (int j=0; j<2; j++) {
			vi z1(part);
			Pref(a.data(), part, z1.data());
			string bha=b;
			bha.push_back('#');
			for(char x : a) bha.push_back(x);
			vi z2(bha.size());
			Pref(bha.data(), bha.size(), z2.data());
			for(auto *v : {&z1, &z2}) {
				v[0][0] = v[0].size();
				v->emplace_back(0);
			}
			for (int c=0; c<(int)a.size(); c++) {
				int l = (int)a.size() - c, x = c - min(l - 1, z1[l]),
					y = c - max(l - z2[(int)b.size() + c + 1], j),
					sb = (j ? end - y - l * 2 : off + x),
					se = (j ? end - x - l * 2 + 1 : off + y + 1),
					&p = pos[l];
				if (x > y) continue;
				if (p != -1 && ans[p].first.second + 1 == sb)
					ans[p].first.second = se - 1;
				else
					p = ans.size(), ans.push_back({{sb, se - 1}, l});
			}
			a = string(b.rbegin(), b.rend());
			b.swap(ra);
		 }
	}
	return ans;
}
