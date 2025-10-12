// #include <iostream>
// #include <functional>
// #include "types.hpp"
// #include "duoram.hpp"
// #include "shapes.hpp"
// #include "test_custompath.hpp"


// // Test program to demonstrate CustomPath usage
// void test_custompath(MPCIO &mpcio, const PRACOptions &opts) {
//     MPCTIO tio(mpcio, 0, opts.num_cpu_threads);

//     run_coroutines(tio, [&tio, &mpcio](yield_t &yield) {
//         int player = tio.player();
//         size_t tree_size = 16;  // Example: tree with 16 nodes

//         std::cout << "[Player " << player << "] Starting CustomPath test" << std::endl;

//         // ========== TEST 1: Duoram with 2-element Tuples ==========
//         std::cout << "\n[Player " << player << "] === Testing Duoram<std::tuple<RegAS, RegAS>> ===" << std::endl;

//         // Create a Duoram of tuples of additively shared values
//         Duoram<std::tuple<RegAS, RegAS>> oram_tuple2(player, tree_size);

//         std::cout << "[Player " << player << "] Created Duoram with 2-element tuples" << std::endl;

//         // Get a Flat view
//         auto flat_tuple2 = oram_tuple2.flat(tio, yield);

//         std::cout << "[Player " << player << "] Got Flat view for 2-element tuples" << std::endl;

//         // Initialize with tuples (i*100, i*100+1)
//         flat_tuple2.init([](size_t i) -> std::tuple<RegAS, RegAS> {
//             RegAS first, second;
//             first.set(i * 100);
//             second.set(i * 100 + 1);
//             return std::make_tuple(first, second);
//         });

//         std::cout << "[Player " << player << "] Initialized flat array with 2-element tuples" << std::endl;

//         // Reconstruct and print tuples
//         auto reconstructed_tuples2 = flat_tuple2.reconstruct();
//         for(size_t i = 0; i < tree_size; ++i) {
//             std::cout << "[Player " << player << "] flat_tuple2[" << i << "] = ("
//                       << std::get<0>(reconstructed_tuples2[i]).share() << ", "
//                       << std::get<1>(reconstructed_tuples2[i]).share() << ")" << std::endl;
//         }

//         // ========== TEST 1b: Duoram with 3-element Tuples ==========
//         std::cout << "\n[Player " << player << "] === Testing Duoram<std::tuple<RegAS, RegAS, RegAS>> ===" << std::endl;

//         // Create a Duoram of 3-element tuples
//         Duoram<std::tuple<RegAS, RegAS, RegAS>> oram_tuple3(player, tree_size);

//         std::cout << "[Player " << player << "] Created Duoram with 3-element tuples" << std::endl;

//         // Get a Flat view
//         auto flat_tuple3 = oram_tuple3.flat(tio, yield);

//         std::cout << "[Player " << player << "] Got Flat view for 3-element tuples" << std::endl;

//         // Initialize with tuples (i*100, i*100+1, i*100+2)
//         flat_tuple3.init([](size_t i) -> std::tuple<RegAS, RegAS, RegAS> {
//             RegAS first, second, third;
//             first.set(i * 100);
//             second.set(i * 100 + 1);
//             third.set(i * 100 + 2);
//             return std::make_tuple(first, second, third);
//         });

//         std::cout << "[Player " << player << "] Initialized flat array with 3-element tuples" << std::endl;

//         // Reconstruct and print tuples
//         auto reconstructed_tuples3 = flat_tuple3.reconstruct();
//         for(size_t i = 0; i < tree_size; ++i) {
//             std::cout << "[Player " << player << "] flat_tuple3[" << i << "] = ("
//                       << std::get<0>(reconstructed_tuples3[i]).share() << ", "
//                       << std::get<1>(reconstructed_tuples3[i]).share() << ", "
//                       << std::get<2>(reconstructed_tuples3[i]).share() << ")" << std::endl;
//         }
//         // ========== TEST 2: Original Duoram with single RegAS ==========
//         std::cout << "\n[Player " << player << "] === Testing Duoram<RegAS> ===" << std::endl;

//         // Create a Duoram of additively shared values
//         Duoram<RegAS> oram(player, tree_size);

//         std::cout << "[Player " << player << "] Created Duoram" << std::endl;

//         // Get a Flat view
//         auto flat = oram.flat(tio, yield);

//         std::cout << "[Player " << player << "] Got Flat view" << std::endl;

//         // Initialize with some values (for testing)
//         flat.init([](size_t i) -> size_t {
//             return i * 100;  // 0, 100, 200, 300, ...
//         });

//         // Reconstruct ONCE, then iterate over the result
//         auto reconstructed_flat = flat.reconstruct();
//         for(size_t i = 0; i < tree_size; ++i) {
//             std::cout << "[Player " << player << "] flat[" << i << "] = " << reconstructed_flat[i].share() << std::endl;
//         }

//         std::cout << "[Player " << player << "] Initialized flat array" << std::endl;

//         if (player == 0) {
//             std::cout << "Testing CustomPath Shape" << std::endl;
//             std::cout << "=========================" << std::endl;
//         }

//         // Define custom parent function
//         // The parent of node at index ind is:
//         // - (ind+1)/2 - 1 if ind is one less than a power of two
//         // - (ind+1)/2 otherwise
//         auto my_compute_parent = [](size_t ind) -> size_t {
//             if (ind <= 1) return 0;  // Root or invalid

//             // Check if ind is one less than a power of two
//             // i.e., check if (ind + 1) is a power of two
//             if (((ind + 1) & ind) == 0) {
//                 // ind is one less than a power of two
//                 return ((ind + 1) / 2) - 1;
//             } else {
//                 return (ind + 1) / 2;
//             }
//         };

//         auto my_compute_parent2 = [](size_t ind) -> size_t {
//             if (ind <= 1) return 0;  // Root or invalid

//             // Check if ind is a power of two
//             if ((ind & (ind - 1)) == 0) {
//                 return (ind / 2);
//             } else {
//                 return (ind - 1) / 2;
//             }
//         };

//         // Test different start indices
//         std::vector<size_t> test_indices = {10, 11, 12};

//         for (size_t start_idx : test_indices) {
//             if (player == 0) {
//                 std::cout << "\nCustomPath from node " << start_idx << ":" << std::endl;
//             }

//             std::cout << "[Player " << player << "] Creating CustomPath from node " << start_idx << std::endl;

//             // Create a CustomPath starting from start_idx with custom parent function
//             typename Duoram<RegAS>::CustomPath cpath(flat, tio, yield, start_idx, my_compute_parent2);

//             std::cout << "[Player " << player << "] CustomPath created, size: " << cpath.size() << std::endl;

//             // Read all values first
//             for (size_t i = 0; i < cpath.size(); ++i) {
//                 auto val = cpath.path_indices[i];
//                 std::cout << "[Player " << player << "]  Path[" << i << "] = " << val << std::endl;
//             }

//             // store and print path values
//             std::vector<RegAS> path_values(cpath.size());
//             for (size_t i = 0; i < cpath.size(); ++i) {
//                 path_values[i] = cpath[i];
//                 auto curVal = mpc_reconstruct(tio, yield, path_values[i]);
//                 std::cout << "[Player " << player << "]  Value[" << i << "] = " << curVal << std::endl;
//             }


//             std::cout << "[Player " << player << "] Completed iteration for start_idx " << start_idx << std::endl;
//         }
//     });
// }
