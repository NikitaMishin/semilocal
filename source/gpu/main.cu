

#include "memory_management.h"
#include "first_iteration/transposition_network_4symbol_gpu.cu"
using namespace memory_management;

int main() {
    auto c_host = allocate_1D_array_aligned_on_cpu<unsigned int>(100);
    auto c_device = allocate_1D_array_aligned_on_gpu<unsigned int>(100);

    copy_from_cpu_to_gpu_sync(c_host,c_device,100);
    copy_from_gpu_to_cpu_sync(c_device,c_host,100);
    four_symbol_gpu_runner_fully_gpu(c_device,100,100,c_device,100,100,c_device,c_device,1);
    synchronize_with_gpu();

    free_1D_array_cpu(c_host);
    free_1D_array_gpu(c_device);
    return 0;
}