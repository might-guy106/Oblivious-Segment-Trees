# Oblivious Segment Tree Implementation

This document describes the oblivious segment tree implementation in the PRAC framework.

### Key Features

- **Range Sum Queries**: Compute the sum of elements in any range [left, right] in O(log n) time
- **Point Updates**: Update individual array elements and propagate changes through the tree


## Usage

### Command Format

```bash
./prac -o -t <threads> <player> [addresses] segmenttree -d <depth> -u <updates> -q <queries>
```

### Parameters

- `-d <depth>`: Depth of the segment tree (default: 5)
  - Creates a tree with 2^depth total nodes
  - Has 2^(depth-1) leaf nodes representing array elements
  - Example: depth 4 creates 16 total nodes with 8 leaf elements

- `-u <updates>`: Number of point updates to perform (default: 1)
  - Updates cycle through array indices 0, 1, 2, ..., 2^(depth-1)-1
  - Update values increment by 50 each time: 50, 100, 150, ...

- `-q <queries>`: Number of range sum queries to perform (default: 1)
  - Queries test different ranges to demonstrate functionality

### Example Usage

Run in three separate terminals (for 3-party computation):

**Terminal 1 (Player 0):**
```bash
./prac -o -t 8 0 segmenttree -d 4 -u 3 -q 5
```

**Terminal 2 (Player 1):**
```bash
./prac -o -t 8 1 localhost segmenttree -d 4 -u 3 -q 5
```

**Terminal 3 (Player 2 - Server):**
```bash
./prac -o -t 8 2 localhost localhost segmenttree -d 4 -u 3 -q 5
```

## Segment Tree Structure

### For Depth 4 Example:

```
Tree Structure (16 total nodes):
                    [1]: 2800
                   /          \
               [2]: 600      [3]: 2200
              /       \      /         \
         [4]: 100  [5]: 500 [6]: 900  [7]: 1300
        /   \     /   \     /   \      /     \
    [8]:0 [9]:100 [10]:200 [11]:300 [12]:400 [13]:500 [14]:600 [15]:700
```

- **Leaf nodes** (indices 8-15): Original array values [0, 100, 200, 300, 400, 500, 600, 700]
- **Internal nodes** (indices 1-7): Contain range sums of their children
- **Root** (index 1): Contains sum of entire array = 2800


## Sample Execution Output

```
(base) ➜  Oblivious-Segment-Trees git:(main) ✗ ./prac -o -t 8 0 segmenttree -d 4 -u 3 -q 5
Segment Tree of depth 4 with 16 nodes created
===== Segment Tree Initialized =====
Depth: , Size: 16
Updates: 3, Queries: 5

===== Initial Segment Tree =====
SegTreeArray[1] = 2800
SegTreeArray[2] = 600
SegTreeArray[3] = 2200
SegTreeArray[4] = 100
SegTreeArray[5] = 500
SegTreeArray[6] = 900
SegTreeArray[7] = 1300
SegTreeArray[8] = 0
SegTreeArray[9] = 100
SegTreeArray[10] = 200
SegTreeArray[11] = 300
SegTreeArray[12] = 400
SegTreeArray[13] = 500
SegTreeArray[14] = 600
SegTreeArray[15] = 700

===== Update 1 begins =====
Index to be updated in the original array = 0
Index to be updated in the segment tree array = 8
Current Value at index = 0
New Value to be updated = 50
Diff = 50
Updated Parent Index = 4 with value = 150
Updated Parent Index = 2 with value = 650
Updated Parent Index = 1 with value = 2850
Update 1 ends

===== Update 2 begins =====
Index to be updated in the original array = 1
Index to be updated in the segment tree array = 9
Current Value at index = 100
New Value to be updated = 100
Diff = 0
Updated Parent Index = 4 with value = 150
Updated Parent Index = 2 with value = 650
Updated Parent Index = 1 with value = 2850
Update 2 ends

===== Update 3 begins =====
Index to be updated in the original array = 2
Index to be updated in the segment tree array = 10
Current Value at index = 200
New Value to be updated = 150
Diff = 18446744073709551566
Updated Parent Index = 5 with value = 450
Updated Parent Index = 2 with value = 600
Updated Parent Index = 1 with value = 2800
Update 3 ends

===== Updated Segment Tree =====
SegTreeArray[1] = 2800
SegTreeArray[2] = 600
SegTreeArray[3] = 2200
SegTreeArray[4] = 150
SegTreeArray[5] = 450
SegTreeArray[6] = 900
SegTreeArray[7] = 1300
SegTreeArray[8] = 50
SegTreeArray[9] = 100
SegTreeArray[10] = 150
SegTreeArray[11] = 300
SegTreeArray[12] = 400
SegTreeArray[13] = 500
SegTreeArray[14] = 600
SegTreeArray[15] = 700

===== Range Sum Query 1 begins =====
Range Sum Query [0, 0]
Level: 3 Left Node [0] (bitVec[ 8]) isincluded: 1 Right Node [0] (bitVec[8]) isincluded: 1 valid: 1
Level: 2 Left Node [0] (bitVec[ 4]) isincluded: 0 Right Node [0] (bitVec[4]) isincluded: 0 valid: 1
Level: 1 Left Node [0] (bitVec[ 2]) isincluded: 0 Right Node [0] (bitVec[2]) isincluded: 0 valid: 1
Level: 0 Left Node [0] (bitVec[ 1]) isincluded: 0 Right Node [0] (bitVec[1]) isincluded: 0 valid: 1
Sum = 50
Range Sum Query 1 ends

===== Range Sum Query 2 begins =====
Range Sum Query [1, 2]
Level: 3 Left Node [1] (bitVec[ 9]) isincluded: 1 Right Node [2] (bitVec[10]) isincluded: 1 valid: 1
Level: 2 Left Node [1] (bitVec[ 5]) isincluded: 0 Right Node [0] (bitVec[4]) isincluded: 0 valid: 0
Level: 1 Left Node [1] (bitVec[ 3]) isincluded: 0 Right Node [0] (bitVec[2]) isincluded: 0 valid: 0
Level: 0 Left Node [0] (bitVec[ 1]) isincluded: 0 Right Node [0] (bitVec[1]) isincluded: 0 valid: 1
Sum = 250
Range Sum Query 2 ends

===== Range Sum Query 3 begins =====
Range Sum Query [2, 4]
Level: 3 Left Node [2] (bitVec[ 10]) isincluded: 0 Right Node [4] (bitVec[12]) isincluded: 1 valid: 1
Level: 2 Left Node [1] (bitVec[ 5]) isincluded: 1 Right Node [1] (bitVec[5]) isincluded: 1 valid: 1
Level: 1 Left Node [1] (bitVec[ 3]) isincluded: 0 Right Node [0] (bitVec[2]) isincluded: 0 valid: 0
Level: 0 Left Node [0] (bitVec[ 1]) isincluded: 0 Right Node [0] (bitVec[1]) isincluded: 0 valid: 1
Sum = 850
Range Sum Query 3 ends

===== Range Sum Query 4 begins =====
Range Sum Query [3, 6]
Level: 3 Left Node [3] (bitVec[ 11]) isincluded: 1 Right Node [6] (bitVec[14]) isincluded: 1 valid: 1
Level: 2 Left Node [2] (bitVec[ 6]) isincluded: 1 Right Node [2] (bitVec[6]) isincluded: 1 valid: 1
Level: 1 Left Node [1] (bitVec[ 3]) isincluded: 0 Right Node [0] (bitVec[2]) isincluded: 0 valid: 0
Level: 0 Left Node [0] (bitVec[ 1]) isincluded: 0 Right Node [0] (bitVec[1]) isincluded: 0 valid: 1
Sum = 1800
Range Sum Query 4 ends

===== Range Sum Query 5 begins =====
Range Sum Query [4, 5]
Level: 3 Left Node [4] (bitVec[ 12]) isincluded: 0 Right Node [5] (bitVec[13]) isincluded: 0 valid: 1
Level: 2 Left Node [2] (bitVec[ 6]) isincluded: 1 Right Node [2] (bitVec[6]) isincluded: 1 valid: 1
Level: 1 Left Node [1] (bitVec[ 3]) isincluded: 0 Right Node [0] (bitVec[2]) isincluded: 0 valid: 0
Level: 0 Left Node [0] (bitVec[ 1]) isincluded: 0 Right Node [0] (bitVec[1]) isincluded: 0 valid: 1
Sum = 900
Range Sum Query 5 ends
2325 messages sent
77402 message bytes sent
2114 Lamport clock (latencies)
9570 local AES operations
14802 milliseconds wall clock time
{14800000000;200000000;500000000} nanoseconds {real;user;system}
Mem: 6396 KiB
Precomputed values used: T0 m:75 a:2 s:120 r1:40 r2:40 r3:80 r4:33 c:20
```
