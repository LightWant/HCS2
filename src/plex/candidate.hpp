#ifndef CANDIDATE_HPP
#define CANDIDATE_HPP

#include <vector>
#include "../tools/types.hpp"
#include <iostream>



class candidate {
private:
    std::vector<ui> pre, nxt, heads, keys;
    ui n;
    ui k;
    ui minKey, maxKey;
    ui sz;

public:
    candidate() {}
    candidate(ui n, ui k):n(n), k(k) {
        pre.resize(n + 1);
        nxt.resize(n + 1);
        heads.resize(n + 1);
        // ids.resize(n + 1);
        keys.resize(n + 1);
        // idx.resize(n + 1);

        for(ui i = 0; i <= n; i++) pre[i] = n;
        for(ui i = 0; i <= n; i++) nxt[i] = n;
        for(ui i = 0; i <= n; i++) heads[i] = n;
        // for(ui i = 0; i <= n; i++) ids[i] = n;
        for(ui i = 0; i <= n; i++) keys[i] = n;
        // for(ui i = 0; i <= n; i++) idx[i] = n;

        minKey = n;
        maxKey = 0; 
        sz = 0; 
    }

    void resize(ui n, ui k) {
        this->n = n;
        this->k = k;

        pre.resize(n + 1);
        nxt.resize(n + 1);
        heads.resize(n + 1);
        // ids.resize(n + 1);
        keys.resize(n + 1);
        // idx.resize(n + 1);

        for(ui i = 0; i <= n; i++) pre[i] = n;
        for(ui i = 0; i <= n; i++) nxt[i] = n;
        for(ui i = 0; i <= n; i++) heads[i] = n;
        // for(ui i = 0; i <= n; i++) ids[i] = n;
        for(ui i = 0; i <= n; i++) keys[i] = n;
        // for(ui i = 0; i <= n; i++) idx[i] = n;

        minKey = n;
        maxKey = 0;  
        sz = 0;


    }

//insert at the head
    void insert(ui v, ui key) {
// printf("insert %u %u %u %u:", v, key, keys[v], n);
        if(keys[v] != n) return;
// printf("insert %u %u:", v, key);
        keys[v] = key;
        // pre[v] = pre[heads[key]]; error! pre[n] may be not n!
        pre[v] = n;
        nxt[v] = heads[key];
        pre[heads[key]] = v;
        heads[key] = v;

        minKey = std::min(minKey, key);
        maxKey = std::max(maxKey, key);

        ++sz;
    }

    void remove(ui v) {
// printf("remove %u, %u, %u, nxt %u\n", v, keys[v], heads[keys[v]], nxt[v]);fflush(stdout);
        if(keys[v] == n) return;
        if(heads[keys[v]] == v) {
            heads[keys[v]] = nxt[v];
        }

        nxt[pre[v]] = nxt[v];
        pre[nxt[v]] = pre[v];

        nxt[v] = pre[v] = n;
        keys[v] = n;

        --sz;
    }

    void increment(ui v) {
        ui newKey = keys[v] + 1;
        remove(v);
        insert(v, newKey);
    }

    void decrement(ui v) {
        ui newKey = keys[v] - 1;
        remove(v);
        insert(v, newKey);
    }

    bool isIn(ui v) { return keys[v] != n; }

    ui operator [](ui v) { return keys[v]; }

    void flush() {
        while(minKey <= maxKey && heads[minKey] == n) ++minKey;
        while(maxKey >= minKey && heads[maxKey] == n) --maxKey;
    }

    bool empty() { return sz == 0; }

    ui size() { return sz; }

    void clear() {
        for(ui key = minKey; key <= maxKey; key++) {
            while(heads[key] != n) {
// printf("min %u max %u, key %u, head %u, n %u\n", minKey, maxKey, key, heads[key], n);
                remove(heads[key]);
            }
        }
        // for(ui i = 0; i <= n; i++) pre[i] = n;
        // for(ui i = 0; i <= n; i++) nxt[i] = n;
        // for(ui i = 0; i <= n; i++) heads[i] = n;
        // // for(ui i = 0; i <= n; i++) ids[i] = n;
        // for(ui i = 0; i <= n; i++) keys[i] = n;
        // // for(ui i = 0; i <= n; i++) idx[i] = n;

        minKey = n;
        maxKey = 0;
        sz = 0;
    }

    const ui toNxt(ui v) const {
        if(nxt[v] != n) return nxt[v];

        ui nxtKey = keys[v] + 1;
        while(nxtKey <= maxKey && heads[nxtKey] == n) nxtKey++;

        return nxtKey > maxKey ? n : heads[nxtKey];
    }

    ui ibegin() { return heads[minKey]; }
    ui ibegin(ui startKey) { return heads[startKey]; }
    ui iend() { return n; }


    struct iter {
        ui v;
        const candidate * c;

        iter(const candidate & c, ui v) {
            this->c = &c;
            this->v = v;
        }
        ui operator++ () { return v = c->toNxt(v); }

        bool operator != (const iter & t) { return v != t.v; }

        ui operator * () { return v; }
    };
    iter begin() { return iter(*this, heads[minKey]); }
    iter begin(ui startKey) { return iter(*this, heads[startKey]); }
    iter end() { return iter(*this, n); }
    
};



#endif