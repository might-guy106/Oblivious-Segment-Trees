#include <iostream>
#include <functional>
#include "types.hpp"
#include "duoram.hpp"
#include "shapes.hpp"
#include "test_custompath.hpp"


// Test program to demonstrate CustomPath usage
void test_custompath(MPCIO &mpcio, const PRACOptions &opts) {
    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);
    
    run_coroutines(tio, [&tio, &mpcio](yield_t &yield) {
        int player = tio.player();
        size_t tree_size = 16;  // Example: tree with 16 nodes
        
        std::cout << "[Player " << player << "] Starting CustomPath test" << std::endl;
        
        // Create a Duoram of additively shared values
        Duoram<RegAS> oram(player, tree_size);
        
        std::cout << "[Player " << player << "] Created Duoram" << std::endl;
        
        // Get a Flat view
        auto flat = oram.flat(tio, yield);
        
        std::cout << "[Player " << player << "] Got Flat view" << std::endl;
        
        // Initialize with some values (for testing)
        flat.init([](size_t i) -> size_t {
            return i * 100;  // 0, 100, 200, 300, ...
        });

        // Reconstruct ONCE, then iterate over the result
        auto reconstructed_flat = flat.reconstruct();
        for(size_t i = 0; i < tree_size; ++i) {
            std::cout << "[Player " << player << "] flat[" << i << "] = " << reconstructed_flat[i].share() << std::endl;
        }
        
        std::cout << "[Player " << player << "] Initialized flat array" << std::endl;
        
        if (player == 0) {
            std::cout << "Testing CustomPath Shape" << std::endl;
            std::cout << "=========================" << std::endl;
        }
        
        // Test different start indices
        std::vector<size_t> test_indices = {10, 11, 12};
        
        for (size_t start_idx : test_indices) {
            if (player == 0) {
                std::cout << "\nCustomPath from node " << start_idx << ":" << std::endl;
            }
            
            std::cout << "[Player " << player << "] Creating CustomPath from node " << start_idx << std::endl;
            
            // Create a CustomPath starting from start_idx
            typename Duoram<RegAS>::CustomPath cpath(flat, tio, yield, start_idx);
            
            std::cout << "[Player " << player << "] CustomPath created, size: " << cpath.size() << std::endl;
                        
            // Read all values first
            for (size_t i = 0; i < cpath.size(); ++i) {
                auto val = cpath.path_indices[i];
                std::cout << "[Player " << player << "]  Path[" << i << "] = " << val << std::endl;
            }

            // store and print path values
            std::vector<RegAS> path_values(cpath.size());
            for (size_t i = 0; i < cpath.size(); ++i) {
                path_values[i] = cpath[i];
                auto curVal = mpc_reconstruct(tio, yield, path_values[i]);
                std::cout << "[Player " << player << "]  Value[" << i << "] = " << curVal << std::endl;
            }

            
            std::cout << "[Player " << player << "] Completed iteration for start_idx " << start_idx << std::endl;
        }
    });
}
