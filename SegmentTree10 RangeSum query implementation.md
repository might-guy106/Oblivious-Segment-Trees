- its an iterative algorthm. it is as follows
```cpp
    long long query(int l, int r) {
        long long sum = 0;
        l += n;
        r += n;

        while (l <= r) {
            if (l & 1) {           // l is a right child
                sum += tree[l];
                l++;
            }
            if (!(r & 1)) {        // r is a left child
                sum += tree[r];
                r--;
            }
            l >>= 1;
            r >>= 1;
        }
        return sum;
    }
```

- to convert it into an oblivous algorithm we need to handle some points they are as follows.
	1. our setup is simliar to other Oblivous Segment tree implementations in Prac. we maintain a Doram for Segment Tree.
	2. we start with Xor shares of l and r (which are our endpoints of query)
	3. we divide furthen process into 2 steps. which are as follows
		1. you compute the path for l to root and r to root in first step. this is done y first converting Xor shares to Additive shares and then adding 1's shares to l and subtracting 1's share from r and then converting them back to Xor shares then we can apply >> 1 to those shares diectly. these become Xor shares of parent so we continue this process till we reach the root. as you can see this is a sequential process this cant be parallelised
		2. now you got the Xor shares and additive shares of the path from l to root and r to root. so this step can be done parallely for each level. In each level you check l <= r using a mpc_compare operation to get a isDone flag. then we get the shares of value present at index l and r using Doram (SegmentTree[l^{additive shares}] and SegmentTree[r^{additive shares}]). then to check if l is a right child and if r is a left child we use the last bit of their XOR shares, which then becomes boolean share of isOdd which denotes if it is a right child or not. so we use these 3 flags isDone, leftLastBit, rightLastBit wtih mpc_and and mpc_flagmult operation to get the additive shares of left contirubiton and right contribution which we directly add to the sum.