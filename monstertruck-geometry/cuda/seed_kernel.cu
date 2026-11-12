// Optional FP64 CUDA acceleration experiment, intentionally not used by
// the default portable geometry kernel. Owned per cached surface grid.
#include <cuda_runtime.h>
#include <cstdint>
#include <cstdlib>
#include <new>

struct GpuSeed {
    uint32_t count;
    uint32_t dims;
    double *grid = nullptr;
    double *query = nullptr;
    uint32_t *result = nullptr;
    double *host_query = nullptr;
    uint32_t *host_result = nullptr;
    cudaStream_t stream = nullptr;
};

__global__ void nearest_fp64(const double *grid, const double *query,
                              uint32_t *result, uint32_t count, uint32_t dims) {
    constexpr int THREADS = 256;
    int index = int(threadIdx.x);
    double winner_cost = INFINITY;
    uint32_t winner = count;
    for (uint32_t i = uint32_t(index); i < count; i += THREADS) {
        const double *coordinates = grid + 4 * i;
        double distance = 0.0;
        for (uint32_t axis = 0; axis < dims; ++axis) {
            const double delta = coordinates[axis] - query[axis];
            distance += delta * delta;
        }
        if (distance < winner_cost || (distance == winner_cost && i < winner)) {
            winner_cost = distance;
            winner = i;
        }
    }
    __shared__ double costs[THREADS];
    __shared__ uint32_t indices[THREADS];
    costs[index] = winner_cost;
    indices[index] = winner;
    __syncthreads();
    for (int stride = THREADS / 2; stride > 0; stride /= 2) {
        if (index < stride) {
            const double other = costs[index + stride];
            const uint32_t other_id = indices[index + stride];
            if (other < costs[index] || (other == costs[index] && other_id < indices[index])) {
                costs[index] = other;
                indices[index] = other_id;
            }
        }
        __syncthreads();
    }
    if (index == 0) result[0] = indices[0];
}

extern "C" void mt_gpu_grid_release(void *handle) {
    auto *grid = static_cast<GpuSeed *>(handle);
    if (!grid) return;
    if (grid->stream) cudaStreamSynchronize(grid->stream);
    if (grid->result) cudaFree(grid->result);
    if (grid->query) cudaFree(grid->query);
    if (grid->grid) cudaFree(grid->grid);
    if (grid->host_query) cudaFreeHost(grid->host_query);
    if (grid->host_result) cudaFreeHost(grid->host_result);
    if (grid->stream) cudaStreamDestroy(grid->stream);
    delete grid;
}

extern "C" int mt_gpu_grid_create(const double *data, uint32_t count,
                                    uint32_t dims, void **output) {
    if (!data || !output || count == 0 || dims == 0 || dims > 4) return -1;
    *output = nullptr;
    auto *grid = new (std::nothrow) GpuSeed;
    if (!grid) return -2;
    grid->count = count;
    grid->dims = dims;
    auto prepare = [&]() {
        if (cudaStreamCreateWithFlags(&grid->stream, cudaStreamNonBlocking) != cudaSuccess) return false;
        if (cudaMalloc(reinterpret_cast<void **>(&grid->grid),
                       size_t(count) * 4 * sizeof(double)) != cudaSuccess) return false;
        if (cudaMalloc(reinterpret_cast<void **>(&grid->query), 4 * sizeof(double)) != cudaSuccess) return false;
        if (cudaMalloc(reinterpret_cast<void **>(&grid->result), sizeof(uint32_t)) != cudaSuccess) return false;
        if (cudaHostAlloc(reinterpret_cast<void **>(&grid->host_query), 4 * sizeof(double),
                          cudaHostAllocDefault) != cudaSuccess) return false;
        if (cudaHostAlloc(reinterpret_cast<void **>(&grid->host_result), sizeof(uint32_t),
                          cudaHostAllocDefault) != cudaSuccess) return false;
        if (cudaMemcpyAsync(grid->grid, data, size_t(count) * 4 * sizeof(double),
                            cudaMemcpyHostToDevice, grid->stream) != cudaSuccess) return false;
        return cudaStreamSynchronize(grid->stream) == cudaSuccess;
    };
    if (!prepare()) {
        mt_gpu_grid_release(grid);
        return -3;
    }
    *output = grid;
    return 0;
}

extern "C" int mt_gpu_grid_nearest(void *handle, const double *target,
                                     uint32_t *winner) {
    auto *grid = static_cast<GpuSeed *>(handle);
    if (!grid || !target || !winner) return -1;
    for (uint32_t i = 0; i < 4; ++i) grid->host_query[i] = target[i];
    if (cudaMemcpyAsync(grid->query, grid->host_query, sizeof(double) * 4,
                        cudaMemcpyHostToDevice, grid->stream) != cudaSuccess) return -2;
    nearest_fp64<<<1, 256, 0, grid->stream>>>(
        grid->grid, grid->query, grid->result, grid->count, grid->dims);
    if (cudaGetLastError() != cudaSuccess) return -3;
    if (cudaMemcpyAsync(grid->host_result, grid->result, sizeof(uint32_t),
                        cudaMemcpyDeviceToHost, grid->stream) != cudaSuccess) return -4;
    if (cudaStreamSynchronize(grid->stream) != cudaSuccess) return -5;
    *winner = *grid->host_result;
    return 0;
}
