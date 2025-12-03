#!/usr/bin/env python3
"""Benchmark scipy.cKDTree vs sklearn.NearestNeighbors for 100k points in 3D."""

import numpy as np
import time
from scipy.spatial import cKDTree
from sklearn.neighbors import NearestNeighbors

# Setup: 100k reference points, 10k query points, 3D, k=5 neighbors
np.random.seed(42)
n_ref = 100000
n_query = 10000
dims = 3
k = 5

X = np.random.rand(n_ref, dims).astype(np.float64)
Y = np.random.rand(n_query, dims).astype(np.float64)

print(f"Benchmark: {n_ref} reference points, {n_query} queries, {dims}D, k={k}")
print("=" * 70)

# Warm-up (compile any JIT, populate caches)
_ = cKDTree(X[:100])
_ = NearestNeighbors(n_neighbors=k).fit(X[:100])

# ========== scipy.spatial.cKDTree ==========
t0 = time.perf_counter()
tree_scipy = cKDTree(X)
t_build_scipy = time.perf_counter() - t0

t0 = time.perf_counter()
D_scipy, idx_scipy = tree_scipy.query(Y, k=k)
t_query_scipy = time.perf_counter() - t0

print(f"scipy.cKDTree:")
print(f"  Build time:  {t_build_scipy*1000:.2f} ms")
print(f"  Query time:  {t_query_scipy*1000:.2f} ms")
print(f"  Total time:  {(t_build_scipy + t_query_scipy)*1000:.2f} ms")
print()

# ========== sklearn (kd_tree algorithm) ==========
t0 = time.perf_counter()
nbrs_kd = NearestNeighbors(n_neighbors=k, algorithm='kd_tree').fit(X)
t_build_kd = time.perf_counter() - t0

t0 = time.perf_counter()
D_kd, idx_kd = nbrs_kd.kneighbors(Y)
t_query_kd = time.perf_counter() - t0

print(f"sklearn.NearestNeighbors (algorithm='kd_tree'):")
print(f"  Build time:  {t_build_kd*1000:.2f} ms")
print(f"  Query time:  {t_query_kd*1000:.2f} ms")
print(f"  Total time:  {(t_build_kd + t_query_kd)*1000:.2f} ms")
print()

# ========== sklearn (ball_tree algorithm) ==========
t0 = time.perf_counter()
nbrs_ball = NearestNeighbors(n_neighbors=k, algorithm='ball_tree').fit(X)
t_build_ball = time.perf_counter() - t0

t0 = time.perf_counter()
D_ball, idx_ball = nbrs_ball.kneighbors(Y)
t_query_ball = time.perf_counter() - t0

print(f"sklearn.NearestNeighbors (algorithm='ball_tree'):")
print(f"  Build time:  {t_build_ball*1000:.2f} ms")
print(f"  Query time:  {t_query_ball*1000:.2f} ms")
print(f"  Total time:  {(t_build_ball + t_query_ball)*1000:.2f} ms")
print()

# ========== Speedup summary ==========
total_scipy = t_build_scipy + t_query_scipy
total_kd = t_build_kd + t_query_kd
total_ball = t_build_ball + t_query_ball

print("=" * 70)
print("Speedup vs scipy.cKDTree:")
print(f"  sklearn kd_tree:   {total_kd/total_scipy:.2f}× slower")
print(f"  sklearn ball_tree: {total_ball/total_scipy:.2f}× slower")
print()
print(f"Winner: scipy.cKDTree ({total_scipy*1000:.2f} ms total)")
