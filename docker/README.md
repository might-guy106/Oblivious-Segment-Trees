# PRAC Docker Setup for Oblivious Segment Trees

This directory contains the Docker setup to run the PRAC framework and Oblivious Segment Tree experiments in a simulated 3-party environment on a single machine.

## How It Works

The setup abstracts the complexity of running three separate parties by containerizing them.

### 1. The Environment (`Dockerfile`)
- Builds a custom Ubuntu 22.04 image.
- Installs all necessary dependencies (`libbsd-dev`, `libboost-all-dev`, etc.).
- **Crucially**: It copies the current source code into the container and compiles it (`make -j`) during the build process.

### 2. The Network (`start-docker`)
- Spins up **3 separate containers** (`prac_p0`, `prac_p1`, `prac_p2`) representing the three computation parties.
- Sets up hostnames and ensures they can communicate.
- Keeps them running in the background, ready to accept commands.

### 3. Orchestration (`run-experiment`)
- This is the main wrapper script you will use.
- **Automatic Discovery**: It automatically checks the internal IP addresses of the running containers.
- **Command Dispatch**: It takes your arguments (e.g., `-t 8 -o ...`) and runs the `prac` binary inside *all three* containers simultaneously.
- **Addressing**: It automatically handles the addressing arguments:
    - Party 0 receives: `./prac ... 0`
    - Party 1 receives: `./prac ... 1 <IP_of_P0>`
    - Party 2 receives: `./prac ... 2 <IP_of_P0> <IP_of_P1>`
- **Output**: It captures stdout/stderr from all three parties and displays them sequentially after execution finishes.

### 4. Network Simulation (`set-networking`)
- Uses Linux Traffic Control (`tc`) to inject artificial latency (default: 30ms) and bandwidth limits (default: 100mbit) on the container network interfaces. This allows you to test how the protocol performs under realistic network conditions without needing actual separate machines.

---

## How to Use

### Step 1: Build the Image
Build the docker image. You only need to do this once, or whenever you modify the C++ source code (since the code is baked into the image).

```bash
cd docker
./build-docker
```

> **Note**: If you change the code, you must rebuild the image (`./build-docker`) for changes to take effect inside the container!

### Step 2: Start the Containers
Launch the three party containers in the background.

```bash
./start-docker
```

### Step 3: Run Experiments (Segment Tree Examples)

You can now run commands using the `run-experiment` script. The arguments are exactly the same as the `./prac` binary, but you **don't** specify the player number or IP addresses.

#### A. Online-Only Mode (Easiest)
To run the Segment Tree demo without separate preprocessing:

```bash
# Format: ./run-experiment -o -t <threads> segmenttree9 -d <depth> -u <updates> -q <queries>
./run-experiment -o segmenttree9 -d 4 -u 3 -q 5
```

#### B. Two-Phase Mode (Preprocessing + Online)
**Phase 1: Preprocessing**
Generate the required resources. `run-experiment` passes arguments to all players; Player 2 handles the resource generation coordination.

```bash
# Example: Generate resources. You'll need the output string from a previous dry-run or calculate needed resources manually.
# For simplicity, you can often run online-only mode first to see what resources are needed.
./run-experiment -t 8 -p 2 <Resource_String> 
```
*Note: Since passing resource strings can be complex, Online-Only mode is often preferred for quick testing.*

### Step 4: Log Management and Analysis
Since logs are stored inside the container, use these helper scripts to manage them.

**1. Clean Logs**
To delete all logs from the container:
```bash
./clean-logs
```

**2. Analyze Results**
To run the performance analysis (inside the container) and copy plots to your host machine:
```bash
./analyze-logs
```
*Plots will be saved to `../logs/plots/` on your host machine.*

> **URGENT**: You MUST rebuild the docker image (`./build-docker`) first! The original image does not have the Python libraries (`pandas`, `matplotlib`, `seaborn`, `rich`) required for analysis.

### Step 5: Network Performance Study (Automated)
To automatically run experiments across different latencies (5ms, 10ms... 30ms) and bandwidths:

```bash
./run-network-study segmenttree9 -d 4 -u 3 -q 5
```
*   This script will clean valid logs, run the experiment for each setting, tag the logs, run the analysis, and save plots to `../logs/plots/network_study/`.
*   You can edit the `LATENCIES` and `BANDWIDTHS` arrays in the script to change the sweep parameters.

### Step 6: Simulate Network Latency (Manual)
To manually set network conditions:

```bash
./set-networking 30ms 100mbit
```

To reset networking:
```bash
./unset-networking
```

### Step 7: Stop Containers
When finished, stop and remove the containers:

```bash
./stop-docker
```

## Troubleshooting
- **Code Changes**: If you edit `../segmentTree9.cpp`, you **must** run `./build-docker` again. The container has its own copy of the code.
- **Permissions**: Ensure scripts have execute permissions (`chmod +x *`).
