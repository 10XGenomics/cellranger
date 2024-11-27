// The file contains reference-free code for finding paths in a graph and
// simplifying it.

#![allow(clippy::many_single_char_names)]

use graph_simple::GraphSimple;
use hyperbase::Hyper;
use std::cmp::max;
use std::iter::zip;
use std::mem::swap;
use vector_utils::{
    bin_member, bin_position1_2, count_instances, erase_if, make_freq, next_diff, next_diff1_2,
    next_diff1_3, reverse_sort, unique_sort,
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// UTILITIES USED BY PATH FINDING AND GRAPH SIMPLIFICATION FUNCTIONS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// umis1: append all UMIs supporting an edge, in order.

pub fn umis1(x: &Hyper, umi_id: &[i32], f: i32, u: &mut Vec<i32>) {
    u.extend(x.ids[f as usize].iter().map(|&id| umi_id[id as usize]));
}

pub fn umis1_u(x: &Hyper, umi_id: &[i32], f: i32, u: &mut Vec<usize>) {
    u.extend(
        x.ids[f as usize]
            .iter()
            .map(|&id| umi_id[id as usize] as usize),
    );
}

// umis2: append all UMIs supporting an edge or its involution.

pub fn umis2(x: &Hyper, umi_id: &[i32], f: i32, u: &mut Vec<i32>) {
    for pass in 0..2 {
        let g = if pass == 1 {
            x.inv[f as usize] as i32
        } else {
            f
        };
        umis1(x, umi_id, g, u);
    }
}

// Determine if at a particular vertex, the majority of reads appear to be fw.

pub fn looks_fw(x: &mut Hyper, v: usize) -> bool {
    let mut fw = 0;
    let mut rc = 0;
    let n = x.h.g.n_from(v);
    for j in 0..n {
        let f = x.h.g.e_from(v, j);
        fw += x.ids[f].len();
        rc += x.ids[x.inv[f] as usize].len();
    }
    fw >= rc
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// PATH FINDING FUNCTIONS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Find strongest edge for each UMI, which is defined to be the edge that has the
// most read support.  In the event of a tie, a choice is made in a manner that is
// independent of graph numbering.

pub fn find_strongest_edges(
    x: &Hyper,
    umi_id: &[i32],
    best: &mut Vec<i32>,
    best_val: &mut Vec<i32>,
) {
    let numi = umi_id.iter().copied().max().map_or(0, |u| u + 1);
    best.resize(numi as usize, -1);
    best_val.resize(numi as usize, 0);
    let mut best_len: Vec<i32> = vec![1000000; numi as usize];
    for e in 0..x.h.g.edge_count() {
        let mut u = Vec::<usize>::new();
        umis1_u(x, umi_id, e as i32, &mut u);
        let mut j = 0;
        loop {
            if j >= u.len() {
                break;
            }
            let k = next_diff(&u, j);
            let n = k - j;
            if best[u[j]] < 0
                || n > best_val[u[j]] as usize
                || (n == best_val[u[j]] as usize && x.h.kmers(e as u32) < best_len[u[j]] as usize)
                || (n == best_val[u[j]] as usize
                    && x.h.kmers(e as u32) == best_len[u[j]] as usize
                    && x.h.g.edge_obj(e as u32).to_string()
                        < x.h.g.edge_obj(best[u[j]] as u32).to_string())
            {
                best[u[j]] = e as i32;
                best_val[u[j]] = n as i32;
                best_len[u[j]] = x.h.kmers(e as u32) as i32;
            }
            j = k;
        }
    }
}

// Find a path for each UMI u:
// 1. Find an edge having maximal support by u.  If the support is less than 50,
//    give up (reporting an empty path for u).
// 2. Extend the edge in both directions, stopping if we encounter an edge having
//    less than 10 support or a branch that cannot be resolving using a 10-to-1
//    support ratio requirement.
// This is a very conservative path definition.

pub fn strong_paths(x: &mut Hyper, umi_id: &[i32], strong: &mut Vec<(i32, Vec<i32>)>) {
    let numi = umi_id.iter().copied().max().map_or(0, |u| u + 1);
    let mut best = Vec::<i32>::new();
    let mut best_val = Vec::<i32>::new();
    find_strongest_edges(x, umi_id, &mut best, &mut best_val);
    const MIN_RAT: i32 = 10;
    const MIN_BEST: i32 = 50;
    strong.clear();
    for u in 0..numi {
        let e = best[u as usize];
        if e < 0 || best_val[u as usize] < MIN_BEST {
            let m = Vec::<i32>::new();
            strong.push((u, m));
            continue;
        }
        let mut z: Vec<i32> = vec![e];
        loop {
            let v = x.h.g.to_right(z[z.len() - 1] as u32);
            let n = x.h.g.n_from(v as usize);
            if n == 0 {
                break;
            }
            let mut supp = Vec::<(i32, i32)>::new();
            for j in 0..n {
                let mut s = 0;
                let e = x.h.g.e_from(v as usize, j);
                for pass in 0..2 {
                    let mut f = e;
                    if pass == 1 {
                        f = x.inv[e] as usize;
                    }
                    // INEFFICIENT!
                    for l in 0..x.ids[f].len() {
                        let id = x.ids[f][l];
                        if umi_id[id as usize] == u {
                            s += 1;
                        }
                    }
                }
                supp.push((s, j as i32));
            }
            reverse_sort(&mut supp);
            if supp[0].0 < MIN_RAT {
                break;
            }
            if n == 1 || MIN_RAT * supp[1].0 <= supp[0].0 {
                let e = x.h.g.e_from(v as usize, supp[0].1 as usize);
                if z.contains(&(e as i32)) {
                    break;
                }
                z.push(e as i32);
            } else {
                break;
            }
        }
        loop {
            let v = x.h.g.to_left(z[0] as u32);
            let n = x.h.g.n_to(v as usize);
            if n == 0 {
                break;
            }
            let mut supp = Vec::<(i32, i32)>::new();
            for j in 0..n {
                let mut s = 0;
                let e = x.h.g.e_to(v as usize, j);
                for pass in 0..2 {
                    let mut f = e;
                    if pass == 1 {
                        f = x.inv[e] as usize;
                    }
                    // INEFFICIENT!
                    for l in 0..x.ids[f].len() {
                        let id = x.ids[f][l];
                        if umi_id[id as usize] == u {
                            s += 1;
                        }
                    }
                }
                supp.push((s, j as i32));
            }
            reverse_sort(&mut supp);
            if supp[0].0 < MIN_RAT {
                break;
            }
            if n == 1 || MIN_RAT * supp[1].0 <= supp[0].0 {
                let e = x.h.g.e_to(v as usize, supp[0].1 as usize);
                if z.contains(&(e as i32)) {
                    break;
                }
                z.insert(0, e as i32);
            } else {
                break;
            }
        }
        strong.push((u, z));
    }
}

// Find strong paths, in a way appropriate for plasma cells where there are very large
// numbers of umis, each having a tiny number of reads.

pub fn alt_strong_paths(x: &mut Hyper, umi_id: &[i32], alt_strong: &mut Vec<Vec<i32>>) {
    alt_strong.clear();
    const MIN_WIN: usize = 10;

    // Go through the edges.

    let mut us = Vec::<i32>::new();
    let numi = umi_id.iter().copied().max().map_or(0, |u| u + 1);
    let mut using = vec![false; numi as usize];
    let mut freq = Vec::<(usize, i32)>::new();
    for e in 0..x.h.g.edge_count() {
        // Find edges supported by many umis.

        us.clear();
        for m in 0..x.ids[e].len() {
            us.push(umi_id[x.ids[e][m] as usize]);
        }
        unique_sort(&mut us);
        if us.len() < MIN_WIN {
            continue;
        }
        for i in 0..us.len() {
            using[us[i] as usize] = true;
        }

        // Start a path using the edge.  Then extend forward using the given umis,
        // always requiring strong evidence at a branch.

        let mut p = vec![e as i32];
        loop {
            let e = p[p.len() - 1];
            let v = x.h.g.to_right(e as u32);
            freq.clear();
            for l in 0..x.h.g.n_from(v as usize) {
                let f = x.h.g.e_from(v as usize, l);
                let mut s = 0;
                for m in 0..x.ids[f].len() {
                    let u = umi_id[x.ids[f][m] as usize];
                    if using[u as usize] {
                        s += 1;
                    }
                }
                freq.push((s, f as i32));
            }
            reverse_sort(&mut freq);
            if freq.is_empty() || freq[0].0 < MIN_WIN {
                break;
            }
            if freq.len() >= 2 && freq[0].0 < MIN_WIN * freq[1].0 {
                break;
            }
            if p.contains(&freq[0].1) {
                break;
            }
            p.push(freq[0].1);
        }

        // Now extend backward in the same way.  Save the path.

        loop {
            let e = p[0];
            let v = x.h.g.to_left(e as u32);
            freq.clear();
            for l in 0..x.h.g.n_to(v as usize) {
                let f = x.h.g.e_to(v as usize, l);
                let mut s = 0;
                for m in 0..x.ids[f].len() {
                    let u = umi_id[x.ids[f][m] as usize];
                    if using[u as usize] {
                        s += 1;
                    }
                }
                freq.push((s, f as i32));
            }
            reverse_sort(&mut freq);
            if freq.is_empty() || freq[0].0 < MIN_WIN {
                break;
            }
            if freq.len() >= 2 && freq[0].0 < MIN_WIN * freq[1].0 {
                break;
            }
            if p.contains(&freq[0].1) {
                break;
            }
            p.insert(0, freq[0].1);
        }
        alt_strong.push(p.clone());
        for i in 0..us.len() {
            using[us[i] as usize] = false;
        }
    }
    unique_sort(alt_strong);
}

// The next function is similar to strong_paths but uses a much more aggressive
// extension rule.  The rule itself is laden with a bunch of heuristics.

pub fn uber_strong_paths(x: &mut Hyper, umi_id: &[i32], strong: &mut Vec<(i32, Vec<i32>)>) {
    let numi = umi_id.iter().copied().max().map_or(0, |u| u + 1);
    let mut best = Vec::<i32>::new();
    let mut best_val = Vec::<i32>::new();
    find_strongest_edges(x, umi_id, &mut best, &mut best_val);
    strong.clear();

    // Create us, showing umi support.  For every edge e, us[e] is {(u,n} where
    // u is a umi and n is the number of reads from that umi supporting the edge.

    let us = {
        let mut us = vec![Vec::<(i32, i32)>::new(); x.h.g.edge_count()];
        let mut ulist = Vec::<i32>::new();
        for (e, umi_e) in us.iter_mut().enumerate().take(x.h.g.edge_count()) {
            ulist.clear();
            ulist.extend(x.ids[e].iter().map(|&m| umi_id[m as usize]));
            ulist.sort_unstable();
            let mut i = 0;
            while i < ulist.len() {
                let j = next_diff(&ulist, i);
                umi_e.push((ulist[i], (j - i) as i32));
                i = j;
            }
        }
        us
    };

    // Main loop.

    for u in 0..numi {
        let mut e = best[u as usize];
        if e < 0 {
            continue;
        }
        let mut z: Vec<i32> = vec![e];
        for _ in 0..2 {
            e = *z.last().unwrap();
            loop {
                let v = x.h.g.to_right(e as u32);
                if x.h.g.n_from(v as usize) == 0 {
                    break;
                }
                let mut freq = Vec::<(usize, usize, usize, usize)>::new();
                let mut xfreq = Vec::<(usize, usize, usize, usize)>::new();
                for l in 0..x.h.g.n_from(v as usize) {
                    let f = x.h.g.e_from(v as usize, l);
                    let mut s = 0_i32;
                    let xs = 0;
                    let mut nonu = 0_i32;
                    let mut g = f;
                    for pass in 0..2 {
                        if pass == 1 {
                            g = x.inv[f] as usize;
                        }
                        nonu += x.ids[g].len() as i32;
                        let p = bin_position1_2(&us[g], &u);
                        if p >= 0 {
                            s += us[g][p as usize].1;
                            nonu -= us[g][p as usize].1;
                        }
                    }
                    freq.push((s as usize, nonu as usize, x.h.kmers(f as u32), f));
                    xfreq.push((xs, s as usize, x.h.kmers(f as u32), f));
                }
                reverse_sort(&mut freq);
                reverse_sort(&mut xfreq);
                let accept = freq.len() == 1
                    || (freq[0].1 > 0 && freq[1].0 + freq[1].1 == 1 && freq[0].3 == xfreq[0].3)
                    || (freq[0].0 >= 4 * max(1, freq[1].0))
                    || (freq[1].0 <= 1 && freq[0].0 >= 2)
                    || (freq[0].0 >= 2 * freq[1].0 && freq[1].1 == 0);
                if !accept {
                    if xfreq[1].0 > 2 || xfreq[0].0 < 4 * max(1, xfreq[1].0) {
                        break;
                    }
                    e = xfreq[0].3 as i32;
                    if z.contains(&e) {
                        break;
                    }
                    z.push(e);
                    continue;
                }
                e = freq[0].3 as i32;
                if z.contains(&e) {
                    break;
                }
                z.push(e);
            }
            z.reverse();
            for v in &mut z {
                *v = x.inv[*v as usize] as i32;
            }
        }
        strong.push((u, z));
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// GRAPH SIMPLIFICATION FUNCTIONS THAT REMOVE READS FROM THE GRAPH
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Simplification based on incompatible edges from the same UMI.  For each UMI u:
// 1. Find its strongpath.
// 2. If any edge e with support by u does not lie in a path that includes the
//    strongpath, delete the support by u of e.

pub fn incompat_clean(x: &mut Hyper, umi_id: &[i32]) {
    // Find strongpath for each UMI.

    let mut strong = Vec::<(i32, Vec<i32>)>::new();
    strong_paths(x, umi_id, &mut strong);

    // Find predecessors and successors for strongpaths.

    let mut lefts = Vec::<i32>::new();
    let mut rights = Vec::<i32>::new();
    for (_, j) in &strong {
        if j.is_empty() {
            continue;
        }
        lefts.push(x.h.g.to_left(*j.first().unwrap() as u32) as i32);
        rights.push(x.h.g.to_right(*j.last().unwrap() as u32) as i32);
    }
    unique_sort(&mut lefts);
    unique_sort(&mut rights);
    let mut p: Vec<Vec<u32>> = Vec::with_capacity(lefts.len());
    let mut s: Vec<Vec<u32>> = Vec::with_capacity(rights.len());
    for &left in &lefts {
        let mut pj = Vec::<u32>::new();
        x.h.g.get_predecessors1(left, &mut pj);
        p.push(pj);
    }
    for &right in &rights {
        let mut sj = Vec::<u32>::new();
        x.h.g.get_successors1(right, &mut sj);
        s.push(sj);
    }

    // Delete support.

    for e in 0..x.h.g.edge_count() {
        let mut to_delete: Vec<bool> = vec![false; x.ids[e].len()];
        let mut j = 0;
        loop {
            if j >= x.ids[e].len() {
                break;
            }
            let id = x.ids[e][j] as usize;
            let u = umi_id[id];
            let mut k = j + 1;
            loop {
                if k == x.ids[e].len() {
                    break;
                }
                if umi_id[x.ids[e][k] as usize] != u {
                    break;
                }
                k += 1;
            }
            if !strong[u as usize].1.is_empty() {
                let y = &strong[u as usize].1;
                let v = x.h.g.to_left(*y.first().unwrap() as u32);
                let w = x.h.g.to_right(*y.last().unwrap() as u32);
                let pre = &p[lefts.binary_search(&(v as i32)).unwrap()];
                let suc = &s[rights.binary_search(&(w as i32)).unwrap()];
                if !y.contains(&(e as i32))
                    && pre.binary_search(&x.h.g.to_right(e as u32)).is_err()
                    && suc.binary_search(&x.h.g.to_left(e as u32)).is_err()
                {
                    for v in &mut to_delete[j..k] {
                        *v = true;
                    }
                }
            }
            j = k;
        }
        erase_if(&mut x.ids[e], &to_delete);
    }
}

// Drop UMIs in the bottom 1% for a given component: for each component, find the
// support by each UMI, and delete the support for that UMI in that component if
// the support lies in the bottom 1% of total support for that component.

pub fn drop_bottom(x: &mut Hyper, umi_id: &[i32]) {
    let mut comp = Vec::<Vec<u32>>::new();
    x.h.g.components_e(&mut comp);
    for c in comp {
        let mut c2: Vec<u32> = c.clone();
        c2.extend(c.iter().map(|&e| x.inv[e as usize]));
        unique_sort(&mut c2);
        let mut u = c2
            .iter()
            .flat_map(|&e| x.ids[e as usize].iter().map(|&id| umi_id[id as usize]))
            .collect::<Vec<_>>();
        u.sort_unstable();
        let mut freq = Vec::<(u32, i32)>::new();
        make_freq(&u, &mut freq);
        let mut delu = Vec::<i32>::new();
        let mut tot = 0;
        for (j, &f) in freq.iter().enumerate() {
            tot += f.0;
            if tot as f64 >= 0.99_f64 * u.len() as f64 {
                delu.extend(freq.iter().skip(j + 1).map(|f| f.1));
                break;
            }
        }
        delu.sort_unstable();
        if delu.is_empty() {
            continue;
        }
        for e in c2 {
            let to_delete: Vec<bool> = x.ids[e as usize]
                .iter()
                .map(|&xid| bin_member(&delu, &umi_id[xid as usize]))
                .collect();
            erase_if(&mut x.ids[e as usize], &to_delete);
        }
    }
}

// For each branch, and each UMI, if one branch has ten times more reads
// supporting it, delete the read support on the weak branch.

pub fn branch_clean(x: &mut Hyper, umi_id: &[i32]) {
    let mut dels = vec![Vec::<i32>::new(); x.h.g.edge_count()];
    for v in 0..x.h.g.node_count() {
        let n = x.h.g.n_from(v);
        if n <= 1 {
            continue;
        }
        let mut xx = Vec::<(i32, i32, i32)>::new();
        let mut yy = Vec::<(i32, i32, i32)>::new();
        for j in 0..n {
            let e = x.h.g.e_from(v, j);
            for pass in 0..2 {
                let mut f = e;
                if pass == 1 {
                    f = x.inv[e] as usize;
                }
                let mut l = 0;
                loop {
                    if l == x.ids[f].len() {
                        break;
                    }
                    let id = x.ids[f][l];
                    let u = umi_id[id as usize];
                    let mut m = l + 1;
                    loop {
                        if m == x.ids[f].len() {
                            break;
                        }
                        if umi_id[x.ids[f][m] as usize] != u {
                            break;
                        }
                        m += 1;
                    }
                    xx.push((u, j as i32, (m - l) as i32));
                    l = m;
                }
            }
        }
        xx.sort_unstable();
        let mut j = 0;
        loop {
            if j == xx.len() {
                break;
            }
            let mut k = j + 1;
            let mut n = 0;
            loop {
                if k == xx.len() {
                    break;
                }
                if xx[k].0 != xx[j].0 || xx[k].1 != xx[j].1 {
                    break;
                }
                k += 1;
            }
            for &l in &xx[j..k] {
                n += l.2;
            }
            yy.push((xx[j].0, n, xx[j].1));
            j = k;
        }
        reverse_sort(&mut yy);
        let mut z = Vec::<(i32, i32)>::new();
        let mut j = 0;
        loop {
            if j == yy.len() {
                break;
            }
            let k = next_diff1_3(&yy, j);
            for l in j + 1..k {
                if yy[j].1 >= 10 * max(1, yy[l].1) {
                    z.push((yy[l].2, yy[l].0));
                }
            }
            j = k;
        }
        z.sort_unstable(); // z = { (branchid, u ) } to kill
        let mut j = 0;
        loop {
            if j == z.len() {
                break;
            }
            let k = next_diff1_2(&z, j);
            let us = z[j..k].iter().map(|(_, l)| *l).collect::<Vec<_>>();
            let b = z[j].0;
            let e = x.h.g.e_from(v, b as usize);
            dels[e].append(&mut us.clone());
            dels[x.inv[e] as usize].append(&mut us.clone());
            j = k;
        }
    }
    for e in 0..x.h.g.edge_count() {
        unique_sort(&mut dels[e]);
        let mut to_delete: Vec<bool> = vec![false; x.ids[e].len()];
        for l in 0..x.ids[e].len() {
            if bin_member(&dels[e], &umi_id[x.ids[e][l] as usize]) {
                to_delete[l] = true;
            }
        }
        erase_if(&mut x.ids[e], &to_delete);
    }
}

// For each UMI, if one component has ten times more reads supporting it from
// that UMI than another component, delete support in the weak component.
// ◼ Unnecessarily O( numis x ncomponents ).

pub fn comp_clean(x: &mut Hyper, umi_id: &[i32]) {
    const MIN_RATIO: i32 = 10;
    let mut comp = Vec::<Vec<u32>>::new();
    x.h.g.components_e(&mut comp);
    let mut to_comp: Vec<i32> = vec![-1; x.h.g.edge_count()];
    for (c, co) in comp.iter().enumerate() {
        for &x in co {
            to_comp[x as usize] = c as i32;
        }
    }
    let numi = umi_id.iter().copied().max().map_or(0, |u| u + 1);
    let mut ccomp = vec![Vec::<i32>::new(); numi as usize];
    // Reserve initial capacity assuming even distribution.
    let mut have = Vec::<u32>::with_capacity(if comp.is_empty() {
        0
    } else {
        umi_id.len() / comp.len()
    });
    let mut havex: Vec<bool> = vec![false; umi_id.len()];
    for (c, co) in comp.iter().enumerate() {
        have.clear();
        for &e in co {
            for pass in 0..2 {
                let mut f = e;
                if pass == 1 {
                    f = x.inv[e as usize];
                }
                for &id in &x.ids[f as usize] {
                    if !havex[id as usize] {
                        have.push(id);
                        havex[id as usize] = true;
                    }
                }
            }
        }
        for &h in &have {
            let id = h as usize;
            ccomp[umi_id[id] as usize].push(c as i32);
            havex[id] = false;
        }
    }
    let mut delu = vec![Vec::<i32>::new(); comp.len()];
    for (u, mut c) in ccomp.into_iter().enumerate().take(numi as usize) {
        if c.is_empty() {
            continue;
        }
        c.sort_unstable();
        let mut supp = Vec::<(i32, i32)>::new();
        let mut i = 0;
        while i < c.len() {
            let j = next_diff(&c, i);
            supp.push(((j - i) as i32, c[i]));
            i = j;
        }
        reverse_sort(&mut supp);
        let s0 = supp[0].0;
        for s in supp.into_iter().skip(1) {
            if s0 >= MIN_RATIO * max(1, s.0) {
                delu[s.1 as usize].push(u as i32);
            }
        }
    }
    for (e, c) in zip(&mut x.ids, to_comp) {
        let c = c as usize;
        if delu[c].is_empty() {
            continue;
        }
        let to_delete: Vec<bool> = e
            .iter()
            .map(|id| bin_member(&delu[c], &umi_id[*id as usize]))
            .collect();
        erase_if(e, &to_delete);
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// GRAPH SIMPLIFICATION FUNCTIONS THAT REMOVE EDGES FROM THE GRAPH
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// If a branch has ten times as many UMIs and ten times as many reads as another
// branch, delete the weak branch.

pub fn power_clean(x: &mut Hyper, umi_id: &[i32]) {
    let mut dels = Vec::<u32>::new();
    for v in 0..x.h.g.node_count() {
        let n = x.h.g.n_from(v);
        if n <= 1 {
            continue;
        }
        const MIN_RATIO: usize = 10;
        for j1 in 0..n {
            for j2 in 0..n {
                if j1 == j2 {
                    continue;
                }
                let e1 = x.h.g.e_from(v, j1);
                let e2 = x.h.g.e_from(v, j2);
                let mut u1 = Vec::<i32>::new();
                let mut u2 = Vec::<i32>::new();
                umis2(x, umi_id, e1 as i32, &mut u1);
                umis2(x, umi_id, e2 as i32, &mut u2);
                let (n1, n2) = (u1.len(), u2.len());
                unique_sort(&mut u1);
                unique_sort(&mut u2);
                if u1.len() >= MIN_RATIO * max(1, u2.len()) && n1 >= MIN_RATIO * max(1, n2) {
                    dels.push(e2 as u32);
                    break;
                }
            }
        }
    }
    x.kill_edges(&dels);
}

// Remove components that have less than 'cap' kmers or less than 5 reads.

pub fn tiny_comp_clean(x: &mut Hyper, cap: i32) {
    let mut dels = Vec::<u32>::new();
    let mut comp = Vec::<Vec<u32>>::new();
    x.h.g.components_e(&mut comp);
    for compj in comp {
        let n: usize = compj.iter().map(|&e| x.h.kmers(e)).sum();
        let mut ids = Vec::<u32>::new();
        for &e in &compj {
            let xid = &x.ids[e as usize];
            let xidr = &x.ids[x.inv[e as usize] as usize];
            ids.reserve(xid.len() + xidr.len());
            ids.extend(xid);
            ids.extend(xidr);
        }
        unique_sort(&mut ids);
        if n < cap as usize || ids.len() < 5_usize {
            dels.extend(compj);
        }
    }
    x.kill_edges(&dels);
}

// Look for branches which have 'essentially single UMI support'
// at one weak branch (see code), then delete these weak branches.
// Note: could only look at branches that are at a SNP.  We could also require that
// branches rejoin eventually, which might be safer.

pub fn solo_clean(x: &mut Hyper, umi_id: &[i32]) {
    let mut dels = Vec::<u32>::new();
    for v in 0..x.h.g.node_count() {
        let n = x.h.g.n_from(v);
        if n <= 1 {
            continue;
        }
        const MIN_RATIO: i32 = 5;
        const MAX_KILL: i32 = 6;
        const MIN_READS_PER_BC: i32 = 10;
        let mut supp: Vec<i32> = Vec::with_capacity(n);
        let mut e: Vec<i32> = Vec::with_capacity(n);
        for j in 0..n {
            let f = x.h.g.e_from(v, j);
            e.push(f as i32);
            supp.push(x.supp(f) as i32);
        }
        let mut ux: Vec<Vec<i32>> = Vec::with_capacity(n);
        let mut good: Vec<bool> = Vec::with_capacity(n);
        let mut count: Vec<i32> = Vec::with_capacity(n);
        for &f in &e {
            let mut u = Vec::<i32>::new();
            umis2(x, umi_id, f, &mut u);
            u.sort_unstable();
            count.push(u.len() as i32);
            let mut uu = Vec::<i32>::new();
            let mut l = 0;
            loop {
                if l == u.len() {
                    break;
                }
                let m = next_diff(&u, l);
                uu.push((m - l) as i32);
                l = m;
            }
            reverse_sort(&mut uu);
            ux.push(uu.clone());
            good.push(uu.len() >= 2 && uu[0] >= MIN_READS_PER_BC && uu[1] >= MIN_READS_PER_BC);
        }
        let mut uall = Vec::<i32>::new();
        for &f in &e {
            umis2(x, umi_id, f, &mut uall);
        }
        uall.sort_unstable();
        for (f, &c) in zip(e, &count) {
            let mut uf = Vec::<i32>::new();
            umis2(x, umi_id, f, &mut uf);
            uf.sort_unstable();
            let mut ucount = 0;
            let mut l = 0;
            loop {
                if l == uf.len() {
                    break;
                }
                let m = next_diff(&uf, l);
                if m - l > MAX_KILL as usize {
                    ucount += 1;
                    l = m;
                    continue;
                }
                let count = count_instances(&uall, &uf[l]);
                if count < MIN_RATIO * (m - l) as i32 {
                    ucount += 1;
                    l = m;
                    continue;
                }
                l = m;
            }
            if ucount <= 1 {
                let ok = zip(&count, &good).any(|(&ck, &goodk)| goodk && ck >= c);
                if ok {
                    dels.push(f as u32);
                }
            }
        }
    }
    x.kill_edges(&dels);
}

// Pop certain simple bubbles.  Suppose:
// 1. One branch is supported by at most one UMI, whereas the second branch
//    is supported by at least two UMIs.
// 2. The total reads supporting the strong branch is at least that of the weak
//    branch.
// 3. The length of the branches is at most K+1 kmers.
// Then we delete the weak branch.

pub fn pop_bubbles(x: &mut Hyper, umi_id: &[i32]) {
    let mut dels = Vec::<u32>::new();
    for v in 0..x.h.g.node_count() {
        if x.h.g.n_to(v) != 1 || x.h.g.n_from(v) != 2 {
            continue;
        }
        let w = x.h.g.v_from(v, 0);
        if x.h.g.v_from(v, 1) != w || x.h.g.n_to(w) != 2 || x.h.g.n_from(w) != 1 {
            continue;
        }
        let mut e1 = x.h.g.e_from(v, 0);
        let mut e2 = x.h.g.e_from(v, 1);
        for pass in 0..2 {
            if pass == 1 {
                swap(&mut e1, &mut e2);
            }
            let mut u1 = Vec::<i32>::new();
            let mut u2 = Vec::<i32>::new();
            umis2(x, umi_id, e1 as i32, &mut u1);
            umis2(x, umi_id, e2 as i32, &mut u2);
            unique_sort(&mut u1);
            unique_sort(&mut u2);
            if u1.len() >= 2
                && u2.len() <= 1
                && x.supp(e1) >= x.supp(e2)
                && x.h.kmers(e1 as u32) <= (x.h.k + 1) as usize
                && x.h.kmers(e2 as u32) <= (x.h.k + 1) as usize
            {
                dels.push(e2 as u32);
                break;
            }
        }
    }
    x.kill_edges(&dels);
}

// For each vertex, first determine whether the edges emanating from it have more
// or less reads than their involutions.  Fix this strong orientation.  Now for
// each emanating edge e, and for each UMI, considering only the strong orientation,
// if e has at most max_kill reads and the total number of reads from all
// the edges is at least min_ratio times the number of reads on e, and this holds
// for all UMIs present on e, kill e.

pub fn simple_simp(x: &mut Hyper, umi_id: &[i32], min_ratio: i32, max_kill: i32) {
    let mut dels = Vec::<u32>::new();
    let numi = umi_id.iter().copied().max().map_or(0, |u| u + 1);
    for v in 0..x.h.g.node_count() {
        let n = x.h.g.n_from(v);
        if n <= 1 {
            continue;
        }
        let use_fw = looks_fw(x, v);
        let mut uallc: Vec<i32> = vec![0; numi as usize];
        for j in 0..n {
            let f = x.h.g.e_from(v, j);
            let mut g = f;
            if !use_fw {
                g = x.inv[f] as usize;
            }
            for l in 0..x.ids[g].len() {
                let id = x.ids[g][l];
                let u = umi_id[id as usize];
                uallc[u as usize] += 1;
            }
        }
        for j in 0..n {
            let f = x.h.g.e_from(v, j);
            let mut uf = Vec::<i32>::new();
            let mut g = f;
            if !use_fw {
                g = x.inv[f] as usize;
            }
            umis1(x, umi_id, g as i32, &mut uf);
            let mut ok: bool = true;
            let mut l = 0;
            loop {
                if l == uf.len() {
                    break;
                }
                let m = next_diff(&uf, l);
                if m - l > max_kill as usize {
                    ok = false;
                    break;
                }
                let count = uallc[uf[l] as usize];
                if count < min_ratio * (m - l) as i32 {
                    ok = false;
                    break;
                }
                l = m;
            }
            if ok {
                dels.push(f as u32);
            }
        }
    }
    x.kill_edges(&dels);
}

// simple_simp_type: very similar to simple_simp, differences flagged
// Given two branches at a vertex, this deletes
// the weak branch if all of the following hold:
// (a) there are at least twice as many reads on the strong branch;
// (b) the weak branch does not have more than 8 reads for any UMI;
// (c) for every UMI, the strong branch has at least as many reads,
//     with at most one possible exception, where it can have just one read.
// This could probably be strengthened.

pub fn simple_simp_type(x: &mut Hyper, umi_id: &[i32]) {
    let mut dels = Vec::<u32>::new();
    for v in 0..x.h.g.node_count() {
        let n = x.h.g.n_from(v);
        if n != 2 {
            continue;
        } // *************************************************
        let use_fw = looks_fw(x, v);
        let mut uall = Vec::<i32>::new();
        for j in 0..n {
            let f = x.h.g.e_from(v, j);
            let mut g = f;
            if !use_fw {
                g = x.inv[f] as usize;
            }
            umis1(x, umi_id, g as i32, &mut uall);
        }
        uall.sort_unstable();
        for j in 0..n {
            let f = x.h.g.e_from(v, j);
            let mut uf = Vec::<i32>::new();
            let mut g = f;
            if !use_fw {
                g = x.inv[f] as usize;
            }
            umis1(x, umi_id, g as i32, &mut uf);
            let mut ok: bool = true;
            // requiring twice as many total ***************************************
            if uall.len() < 3 * uf.len() {
                ok = false;
            } // **********************
            let mut exceptions = 0; // *********************************************
            let mut l = 0;
            loop {
                if l == uf.len() {
                    break;
                }
                let m = next_diff(&uf, l);
                if m - l > 8 {
                    // **************************************************
                    ok = false;
                    break;
                }
                let count = count_instances(&uall, &uf[l]);
                if count < 2 * (m - l) as i32 {
                    // ***********************************
                    exceptions += m - l; // ******************************************
                } // ***************************************************************
                l = m;
            }
            if exceptions > 1 {
                ok = false;
            } // ONLY ONE EXCEPTION ****************
            if ok {
                dels.push(f as u32);
            }
        }
    }
    x.kill_edges(&dels);
}

// Kill edges that are not on a UMI path, as defined by uber_strong_paths
// and alt_strong_paths.

pub fn path_clean(x: &mut Hyper, umi_id: &[i32]) {
    let mut strong = Vec::<(i32, Vec<i32>)>::new();
    uber_strong_paths(x, umi_id, &mut strong);
    let mut alt_strong = Vec::<Vec<i32>>::new();
    alt_strong_paths(x, umi_id, &mut alt_strong);
    let all = strong.into_iter().map(|(_, s)| s).chain(alt_strong);
    let mut keep = Vec::<u32>::new();
    for mut z in all {
        keep.reserve(z.len() * 4);
        for pass in 0..2 {
            if pass == 1 {
                z.reverse();
                for zj in &mut z {
                    *zj = x.inv[*zj as usize] as i32;
                }
            }
            for &zl in &z {
                keep.push(zl as u32);
                keep.push(x.inv[zl as usize]);
            }
        }
    }
    unique_sort(&mut keep);
    let dels: Vec<u32> = (0..x.h.g.edge_count() as u32)
        .filter(|e| !bin_member(&keep, e))
        .collect();
    x.kill_edges(&dels);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// MASTER GRAPH CLEANER
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Combo cleaner: the master reference-free cleaner.  The only justification for
// the particular order of operations and arguments is that on our test set of a
// dozen or so samples, and relative to the metrics we used to evaluate them, this
// was the best choice amongst those tested.

pub fn simplify_without_ref(x: &mut Hyper, umi_id: &[i32]) {
    let verbose = false;
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }

    power_clean(x, umi_id);

    comp_clean(x, umi_id);

    if verbose {
        println!("\ncalling simple_simp.1");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    simple_simp(x, umi_id, 8, 8);

    if verbose {
        println!("\ncalling comp_clean.2");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    comp_clean(x, umi_id);

    if verbose {
        println!("\ncalling simple_simp.2");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    simple_simp(x, umi_id, 8, 8);

    if verbose {
        println!("\ncalling branch_clean.1");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    branch_clean(x, umi_id);

    if verbose {
        println!("\ncalling simple_simp.3");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    simple_simp(x, umi_id, 8, 8);

    if verbose {
        println!("\ncalling solo_clean");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    solo_clean(x, umi_id);

    simple_simp(x, umi_id, 5, 50);
    simple_simp(x, umi_id, 3, 8);
    if verbose {
        println!("\ncalling branch_clean.2");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    branch_clean(x, umi_id);

    if verbose {
        println!("\nbefore incompat_clean");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    incompat_clean(x, umi_id);

    if verbose {
        println!("\nafter incompat_clean");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    simple_simp(x, umi_id, 5, 50);

    if verbose {
        println!("\ncalling drop_bottom");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    drop_bottom(x, umi_id);
    if verbose {
        println!("after drop_bottom");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }

    simple_simp(x, umi_id, 2, 8);

    simple_simp_type(x, umi_id);

    if verbose {
        println!("\ncalling path_clean");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    path_clean(x, umi_id);

    if verbose {
        println!("\ncalling pop_bubbles");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
    pop_bubbles(x, umi_id);

    power_clean(x, umi_id);

    tiny_comp_clean(x, 150);

    if verbose {
        println!("\nnonref-cleaning complete");
    }
    if verbose {
        printme!(x.h.g.edge_count(), x.checksum());
    }
}
