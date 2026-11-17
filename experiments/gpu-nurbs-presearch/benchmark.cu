// Focused CUDA feasibility study for Monstertruck's cached 51x51 NURBS
// inverse-projection seed search. No CUDA is used by production code.
// Each query scans the same 2,601 precomputed f64 positions, taking the
// lexicographically first (distance, original row-major index) minimum.
#include <cuda_runtime.h>
#include <cfloat>
#include <functional>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

constexpr int kGridWidth = 51;
constexpr int kGrid = kGridWidth * kGridWidth;
constexpr int kThreads = 256;
constexpr int kMaxQueries = 8192;

struct Point {
    double x, y, z;
};

struct Point32 {
    float x, y, z;
};

#define CUDA_OK(expr) do { cudaError_t e = (expr); if (e != cudaSuccess) { \
  std::cerr << #expr << ": " << cudaGetErrorString(e) << "\n"; std::exit(3); } } while (false)

__host__ __device__ inline double distance2(Point a, Point b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

// Spatial-index control experiment: this is the *same* FP64 argmin,
// including earliest-index ties, with subtree AABB pruning. Build once
// per cached grid, not once per query.
class GridKdTree {
    struct Node {
        int sample, axis, left, right;
        Point lo, hi;
    };
    const Point *grid_;
    std::vector<Node> nodes_;
    int root_;

    static double axis(Point point, int direction) {
        if (direction == 0) return point.x;
        if (direction == 1) return point.y;
        return point.z;
    }
    static void extend(Point &lo, Point &hi, Point p) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    int build(std::vector<int> &ids, int begin, int end) {
        if (begin == end) return -1;
        Point lo = grid_[ids[begin]];
        Point hi = lo;
        for (int i = begin + 1; i < end; ++i)
            extend(lo, hi, grid_[ids[i]]);
        const Point span{hi.x - lo.x, hi.y - lo.y, hi.z - lo.z};
        const int split_axis = span.y > span.x && span.y >= span.z
            ? 1 : span.z > span.x && span.z > span.y ? 2 : 0;
        const int mid = begin + (end - begin) / 2;
        std::nth_element(ids.begin() + begin, ids.begin() + mid,
                         ids.begin() + end, [&](int a, int b) {
            const double x = axis(grid_[a], split_axis);
            const double y = axis(grid_[b], split_axis);
            return x < y || (x == y && a < b);
        });
        const int node_id = static_cast<int>(nodes_.size());
        nodes_.push_back(Node{ids[mid], split_axis, -1, -1, lo, hi});
        const int left = build(ids, begin, mid);
        const int right = build(ids, mid + 1, end);
        nodes_[node_id].left = left;
        nodes_[node_id].right = right;
        return node_id;
    }
    static double lower_bound_dist2(Point p, const Node &node) {
        auto dist = [](double q, double low, double high) {
            return q < low ? low - q : q > high ? q - high : 0.0;
        };
        const double x = dist(p.x, node.lo.x, node.hi.x);
        const double y = dist(p.y, node.lo.y, node.hi.y);
        const double z = dist(p.z, node.lo.z, node.hi.z);
        return x * x + y * y + z * z;
    }
    void search(int node_index, Point query, double &best, int &result) const {
        if (node_index < 0) return;
        const Node &node = nodes_[node_index];
        const double distance = distance2(grid_[node.sample], query);
        if (distance < best || (distance == best && node.sample < result)) {
            best = distance;
            result = node.sample;
        }
        const bool left_first =
            axis(query, node.axis) < axis(grid_[node.sample], node.axis);
        const int near_child = left_first ? node.left : node.right;
        const int far_child = left_first ? node.right : node.left;
        search(near_child, query, best, result);
        if (far_child >= 0) {
            // The AABB is an exact geometric lower bound. Keep a generous
            // *extra search* guard against floating-point rounding when
            // deciding if pruning is valid; this can only cost extra work.
            const double safe_bound =
                best + std::max(1e-12, 64.0 * DBL_EPSILON * best);
            if (lower_bound_dist2(query, nodes_[far_child]) <= safe_bound)
                search(far_child, query, best, result);
        }
    }
public:
    explicit GridKdTree(const Point *grid) : grid_(grid), root_(-1) {
        std::vector<int> ids(kGrid);
        std::iota(ids.begin(), ids.end(), 0);
        nodes_.reserve(kGrid);
        root_ = build(ids, 0, kGrid);
    }
    int nearest(Point query) const {
        double best = INFINITY;
        int result = 0;
        search(root_, query, best, result);
        return result;
    }
};

// Original CPU scan uses strict-<, so exact ties pick the earliest sample.
int cpu_argmin(const Point *grid, Point query) {
    double best = INFINITY;
    int selected = 0;
    for (int i = 0; i < kGrid; ++i) {
        const double dist = distance2(grid[i], query);
        if (dist < best) {
            best = dist;
            selected = i;
        }
    }
    return selected;
}

__device__ inline bool better(double a, int ia, double b, int ib) {
    return a < b || (a == b && ia < ib);
}

// One CTA per query. Grid stays resident on the GPU between calls;
// query coordinates and winning row-major indices cross PCIe.
__global__ void nearest_grid(const Point *__restrict__ grid,
                             const Point *__restrict__ queries,
                             int *__restrict__ winners, int count) {
    const int q = blockIdx.x;
    if (q >= count) return;

    const Point point = queries[q];
    double best = INFINITY;
    int best_index = kGrid;
    for (int i = threadIdx.x; i < kGrid; i += blockDim.x) {
        const double dist = distance2(grid[i], point);
        if (dist < best || (dist == best && i < best_index)) {
            best = dist;
            best_index = i;
        }
    }

    __shared__ double costs[kThreads];
    __shared__ int indices[kThreads];
    costs[threadIdx.x] = best;
    indices[threadIdx.x] = best_index;
    __syncthreads();

    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
        if (threadIdx.x < stride) {
            const double other = costs[threadIdx.x + stride];
            const int other_index = indices[threadIdx.x + stride];
            if (better(other, other_index, costs[threadIdx.x], indices[threadIdx.x])) {
                costs[threadIdx.x] = other;
                indices[threadIdx.x] = other_index;
            }
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        winners[q] = indices[0] == kGrid ? 0 : indices[0];
    }
}


// Upper-bound experiment: float32 *local coordinates* on the same grid.
// A CUDA float32 winner is NOT an exact substitute for the float64 result.
__global__ void nearest_grid_float(const Point32 *__restrict__ grid,
                                   const Point32 *__restrict__ queries,
                                   int *__restrict__ winners, int count) {
    const int q = blockIdx.x;
    if (q >= count) return;
    const Point32 point = queries[q];
    float best = INFINITY;
    int best_index = kGrid;
    for (int i = threadIdx.x; i < kGrid; i += blockDim.x) {
        const Point32 sample = grid[i];
        const float dx = sample.x - point.x;
        const float dy = sample.y - point.y;
        const float dz = sample.z - point.z;
        const float dist = dx * dx + dy * dy + dz * dz;
        if (dist < best || (dist == best && i < best_index)) {
            best = dist;
            best_index = i;
        }
    }
    __shared__ float costs[kThreads];
    __shared__ int indices[kThreads];
    costs[threadIdx.x] = best;
    indices[threadIdx.x] = best_index;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
        if (threadIdx.x < stride) {
            const float other = costs[threadIdx.x + stride];
            const int other_index = indices[threadIdx.x + stride];
            if (other < costs[threadIdx.x] ||
                (other == costs[threadIdx.x] && other_index < indices[threadIdx.x])) {
                costs[threadIdx.x] = other;
                indices[threadIdx.x] = other_index;
            }
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
        winners[q] = indices[0] == kGrid ? 0 : indices[0];
}

uint64_t next_random(uint64_t &state) {
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * 2685821657736338717ULL;
}

double unit_random(uint64_t &state) {
    return static_cast<double>(next_random(state) >> 11) * (1.0 / 9007199254740992.0);
}

void populate(Point *grid, Point *queries) {
    // An extruded curved surface around the actual source-part coordinate
    // scale (approximately x=67, y=5.3, z=2). This is a *representative*
    // smooth grid, not a trace of the source STEP or NURBS basis evaluator.
    for (int u = 0; u < kGridWidth; ++u) {
        for (int v = 0; v < kGridWidth; ++v) {
            const double a = 0.12 + 1.5 * (double(v) / 50.0);
            grid[u * kGridWidth + v] = {
                67.1195 + 1.27 * (double(u) / 50.0),
                5.36457203235 + 0.6731 * std::cos(a),
                1.67943412145 + 0.6731 * std::sin(a),
            };
        }
    }
    uint64_t state = 0xF14C6356BEEFCAFEULL;
    for (int q = 0; q < kMaxQueries; ++q) {
        const int index = static_cast<int>(next_random(state) % kGrid);
        const Point p = grid[index];
        // Nominal probe within a few microns to a few hundredths of a mm.
        queries[q] = {
            p.x + (unit_random(state) - 0.5) * 0.04,
            p.y + (unit_random(state) - 0.5) * 0.04,
            p.z + (unit_random(state) - 0.5) * 0.04,
        };
    }
    // Exact ties and exact samples exercise deterministic first-index logic.
    queries[0] = grid[0];
    queries[1] = grid[kGrid - 1];
    // Deliberately probe near exact Voronoi boundaries. At the 0.0254 mm
    // grid spacing, a 1e-10 mm perturbation remains visible in f64 but
    // generally disappears when f32 receives local-space coordinates.
    for (int q = kMaxQueries / 2; q < kMaxQueries; ++q) {
        const int u = (q - kMaxQueries / 2) % 50;
        const int v = (q / 50) % 51;
        const Point first = grid[u * kGridWidth + v];
        const Point second = grid[(u + 1) * kGridWidth + v];
        queries[q] = {
            (first.x + second.x) / 2.0 +
                ((q & 1) ? 1e-10 : -1e-10),
            first.y, first.z
        };
    }
}

using Clock = std::chrono::steady_clock;
double milliseconds(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
double median(std::vector<double> runs) {
    std::sort(runs.begin(), runs.end());
    return runs[runs.size() / 2];
}

double cpu_time(const Point *grid, const Point *queries, int *out,
                int count, int threads, int repeats) {
    std::vector<double> samples;
    for (int repeat = 0; repeat < repeats; ++repeat) {
        const auto begin = Clock::now();
        if (threads == 1) {
            for (int i = 0; i < count; ++i)
                out[i] = cpu_argmin(grid, queries[i]);
        } else {
            // OpenMP retains its worker pool, unlike spawning eight OS threads
            // on every call. Include the real batch-dispatch/barrier cost.
            #pragma omp parallel for num_threads(threads) schedule(static)
            for (int i = 0; i < count; ++i)
                out[i] = cpu_argmin(grid, queries[i]);
        }
        samples.push_back(milliseconds(begin));
    }
    return median(samples);
}

int main() {
    cudaDeviceProp props{};
    CUDA_OK(cudaGetDeviceProperties(&props, 0));
    std::cout << "GPU: " << props.name << ", SM " << props.major << "."
              << props.minor << ", grid=" << kGrid << " f64 points"
              << ", CPU threads=1/8\n";

    Point *host_grid = nullptr, *host_queries = nullptr;
    int *host_winners = nullptr;
    CUDA_OK(cudaMallocHost(&host_grid, kGrid * sizeof(Point)));
    CUDA_OK(cudaMallocHost(&host_queries, kMaxQueries * sizeof(Point)));
    CUDA_OK(cudaMallocHost(&host_winners, kMaxQueries * sizeof(int)));
    populate(host_grid, host_queries);
    Point32 *host_grid32 = nullptr, *host_queries32 = nullptr;
    int *host_winners32 = nullptr;
    CUDA_OK(cudaMallocHost(&host_grid32, kGrid * sizeof(Point32)));
    CUDA_OK(cudaMallocHost(&host_queries32, kMaxQueries * sizeof(Point32)));
    CUDA_OK(cudaMallocHost(&host_winners32, kMaxQueries * sizeof(int)));
    const Point origin = host_grid[0];
    for (int i = 0; i < kGrid; ++i)
        host_grid32[i] = {
            float(host_grid[i].x - origin.x),
            float(host_grid[i].y - origin.y),
            float(host_grid[i].z - origin.z)
        };
    for (int i = 0; i < kMaxQueries; ++i)
        host_queries32[i] = {
            float(host_queries[i].x - origin.x),
            float(host_queries[i].y - origin.y),
            float(host_queries[i].z - origin.z)
        };
    Point *device_grid = nullptr, *device_queries = nullptr;
    int *device_winners = nullptr;
    CUDA_OK(cudaMalloc(&device_grid, kGrid * sizeof(Point)));
    CUDA_OK(cudaMalloc(&device_queries, kMaxQueries * sizeof(Point)));
    CUDA_OK(cudaMalloc(&device_winners, kMaxQueries * sizeof(int)));
    Point32 *device_grid32 = nullptr, *device_queries32 = nullptr;
    int *device_winners32 = nullptr;
    CUDA_OK(cudaMalloc(&device_grid32, kGrid * sizeof(Point32)));
    CUDA_OK(cudaMalloc(&device_queries32, kMaxQueries * sizeof(Point32)));
    CUDA_OK(cudaMalloc(&device_winners32, kMaxQueries * sizeof(int)));

    cudaStream_t stream{};
    CUDA_OK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    CUDA_OK(cudaMemcpyAsync(device_grid, host_grid, kGrid * sizeof(Point),
                            cudaMemcpyHostToDevice, stream));
    CUDA_OK(cudaMemcpyAsync(device_grid32, host_grid32, kGrid * sizeof(Point32),
                            cudaMemcpyHostToDevice, stream));
    CUDA_OK(cudaStreamSynchronize(stream));
    nearest_grid<<<1, kThreads, 0, stream>>>(device_grid, device_queries, device_winners, 1);
    CUDA_OK(cudaGetLastError());
    CUDA_OK(cudaStreamSynchronize(stream));

    std::vector<int> cpu(kMaxQueries);
    const auto tree_build_start = Clock::now();
    const GridKdTree tree(host_grid);
    const double tree_build_us = milliseconds(tree_build_start) * 1000.0;
    std::cout << "kd_tree_build_us=" << std::fixed << std::setprecision(2)
              << tree_build_us << "\n";
    omp_set_dynamic(0);
    cpu_time(host_grid, host_queries, cpu.data(), 128, 8, 2);
    std::cout << "batch,cpu1_us,cpu8_us,kd1_us,kd8_us,gpu64_e2e_us,gpu64_resident_us,vs_cpu1,vs_cpu8,mismatches64,gpu32_e2e_us,gpu32_resident_us,mismatches32,mismatches_kd\n";
    for (int count : {1, 8, 32, 128, 512, 2048, 8192}) {
        // The CPU reference and GPU consume exactly the same query coordinates.
        const int cpu_repeats = count <= 512 ? 9 : 5;
        const double cpu1 = cpu_time(host_grid, host_queries, cpu.data(),
                                     count, 1, cpu_repeats);
        const double cpu8 = count == 1
            ? cpu1
            : cpu_time(host_grid, host_queries, cpu.data(),
                       count, 8, cpu_repeats);
        std::vector<double> kd1_samples, kd8_samples;
        std::vector<int> kd_winners(count);
        for (int repeat = 0; repeat < cpu_repeats; ++repeat) {
            const auto start = Clock::now();
            for (int i = 0; i < count; ++i)
                kd_winners[i] = tree.nearest(host_queries[i]);
            kd1_samples.push_back(milliseconds(start));
        }
        for (int repeat = 0; repeat < cpu_repeats; ++repeat) {
            const auto start = Clock::now();
            #pragma omp parallel for num_threads(8) schedule(static)
            for (int i = 0; i < count; ++i)
                kd_winners[i] = tree.nearest(host_queries[i]);
            kd8_samples.push_back(milliseconds(start));
        }
        int mismatches_kd = 0;
        for (int i = 0; i < count; ++i)
            mismatches_kd += cpu[i] != kd_winners[i];
        std::vector<double> e2e;
        std::vector<double> resident;
        const int repetitions = count <= 32 ? 80 : count <= 512 ? 40 : 12;
        for (int repetition = 0; repetition < repetitions; ++repetition) {
            const auto begin = Clock::now();
            CUDA_OK(cudaMemcpyAsync(device_queries, host_queries,
                                    count * sizeof(Point),
                                    cudaMemcpyHostToDevice, stream));
            nearest_grid<<<count, kThreads, 0, stream>>>(
                device_grid, device_queries, device_winners, count);
            CUDA_OK(cudaGetLastError());
            CUDA_OK(cudaMemcpyAsync(host_winners, device_winners,
                                    count * sizeof(int),
                                    cudaMemcpyDeviceToHost, stream));
            CUDA_OK(cudaStreamSynchronize(stream));
            e2e.push_back(milliseconds(begin));
        }
        for (int repetition = 0; repetition < repetitions; ++repetition) {
            const auto begin = Clock::now();
            nearest_grid<<<count, kThreads, 0, stream>>>(
                device_grid, device_queries, device_winners, count);
            CUDA_OK(cudaGetLastError());
            CUDA_OK(cudaStreamSynchronize(stream));
            resident.push_back(milliseconds(begin));
        }
        CUDA_OK(cudaMemcpyAsync(host_winners, device_winners, count * sizeof(int),
                                cudaMemcpyDeviceToHost, stream));
        CUDA_OK(cudaStreamSynchronize(stream));
        std::vector<double> e2e32;
        std::vector<double> resident32;
        for (int repetition = 0; repetition < repetitions; ++repetition) {
            const auto begin = Clock::now();
            CUDA_OK(cudaMemcpyAsync(device_queries32, host_queries32,
                                    count * sizeof(Point32),
                                    cudaMemcpyHostToDevice, stream));
            nearest_grid_float<<<count, kThreads, 0, stream>>>(
                device_grid32, device_queries32, device_winners32, count);
            CUDA_OK(cudaGetLastError());
            CUDA_OK(cudaMemcpyAsync(host_winners32, device_winners32,
                                    count * sizeof(int),
                                    cudaMemcpyDeviceToHost, stream));
            CUDA_OK(cudaStreamSynchronize(stream));
            e2e32.push_back(milliseconds(begin));
        }
        for (int repetition = 0; repetition < repetitions; ++repetition) {
            const auto begin = Clock::now();
            nearest_grid_float<<<count, kThreads, 0, stream>>>(
                device_grid32, device_queries32, device_winners32, count);
            CUDA_OK(cudaGetLastError());
            CUDA_OK(cudaStreamSynchronize(stream));
            resident32.push_back(milliseconds(begin));
        }
        CUDA_OK(cudaMemcpyAsync(host_winners32, device_winners32, count * sizeof(int),
                                cudaMemcpyDeviceToHost, stream));
        CUDA_OK(cudaStreamSynchronize(stream));
        int mismatch32 = 0;
        for (int i = 0; i < count; ++i)
            mismatch32 += cpu[i] != host_winners32[i];
        int mismatch = 0;
        for (int i = 0; i < count; ++i) {
            mismatch += cpu[i] != host_winners[i];
        }
        const double us_cpu1 = cpu1 * 1000.0;
        const double us_cpu8 = cpu8 * 1000.0;
        const double us_gpu = median(e2e) * 1000.0;
        std::cout << count << "," << std::fixed << std::setprecision(2)
                  << us_cpu1 << "," << us_cpu8 << ","
                  << median(kd1_samples) * 1000.0 << ","
                  << median(kd8_samples) * 1000.0 << "," << us_gpu << ","
                  << median(resident) * 1000.0 << ","
                  << us_cpu1 / us_gpu << "," << us_cpu8 / us_gpu << ","
                  << mismatch << ","
                  << median(e2e32) * 1000.0 << ","
                  << median(resident32) * 1000.0 << ","
                  << mismatch32 << "," << mismatches_kd << "\n";
        if (mismatch || mismatches_kd) {
            std::cerr << "MISMATCH: " << mismatch << " queries at batch=" << count << "\n";
            return 2;
        }
    }
    CUDA_OK(cudaStreamDestroy(stream));
    CUDA_OK(cudaFree(device_grid));
    CUDA_OK(cudaFree(device_queries));
    CUDA_OK(cudaFree(device_winners));
    CUDA_OK(cudaFree(device_grid32));
    CUDA_OK(cudaFree(device_queries32));
    CUDA_OK(cudaFree(device_winners32));
    CUDA_OK(cudaFreeHost(host_grid32));
    CUDA_OK(cudaFreeHost(host_queries32));
    CUDA_OK(cudaFreeHost(host_winners32));
    CUDA_OK(cudaFreeHost(host_grid));
    CUDA_OK(cudaFreeHost(host_queries));
    CUDA_OK(cudaFreeHost(host_winners));
    return 0;
}
