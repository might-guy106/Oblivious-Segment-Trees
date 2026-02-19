# Oblivious Segment Tree Implementation

This document describes the oblivious segment tree implementation in the PRAC framework. The code is split into `segmentTreeBasic.cpp`/`segmentTreeBasic.hpp` (basic implementation) and `segmentTreeOpt.cpp`/`segmentTreeOpt.hpp` (optimized implementation).

## Setup
**Option A: Local Setup**
1.  **Clone the repository:**
    ```bash
    git clone git@github.com:might-guy106/Oblivious-Segment-Trees.git
    cd Oblivious-Segment-Trees
    ```

2.  **Build the project:**
    ```bash
    make
    ```

3.  **Python Environment Setup:**
    The performance analysis script requires Python packages (`pandas`, `matplotlib`, `seaborn`, `rich`).

    
    We recommend using `uv` for environment management.

    1.  Install uv (if not installed):
        ```bash
        curl -LsSf https://astral.sh/uv/install.sh | sh
        ```
    2.  Create a virtual environment and install dependencies:
        ```bash
        uv venv
        source .venv/bin/activate
        uv pip install pandas matplotlib seaborn rich
        ```

**Option B: Server Setup (172.27.21.172)**
If you are running on the specific server, the environment is already pre-configured.

1.  SSH into the server:
    ```bash
    ssh segmentrees@172.27.21.172
    ```
2.  Navigate to the project directory:
    ```bash
    cd /home/segmentrees/ugp/Oblivious-Segment-Trees
    ```
3.  **Build the project:**
    ```bash
    make
    ```
4.  Activate the Python environment:
    ```bash
    source myenv/bin/activate
    ```

## Manual Execution

You can run the experiment in two modes: **Online Only Mode** or **Two-Phase Mode** (Preprocessing + Online).

### Command Format

```bash
./prac -o -t <threads> <player> [addresses] segmenttreebasic -d <depth> -u <updates> -q <queries>
```
or for the optimized version:
```bash
./prac -o -t <threads> <player> [addresses] segmenttreeopt -d <depth> -u <updates> -q <queries>
```

### Parameters

-   `-d <depth>`: Depth of the segment tree (default: 5)
    -   Creates a tree with $2^{depth}$ total nodes.
    -   Has $2^{depth-1}$ leaf nodes representing array elements.
    -   Example: depth 4 creates 16 total nodes with 8 leaf elements.
-   `-u <updates>`: Number of point updates to perform (default: 1)
    -   Updates cycle through array indices 0, 1, 2, ..., $2^{depth-1}-1$.
    -   Update values increment by 50 each time: 50, 100, 150, ...
-   `-q <queries>`: Number of range sum queries to perform (default: 1)
    -   Queries test different ranges to demonstrate functionality.

### 1. Online Only Mode

Run in three separate terminals (for 3-party computation):

**Terminal 1 (Player 0):**
```bash
./prac -o -t 8 0 segmenttreebasic -d 4 -u 3 -q 5
```

**Terminal 2 (Player 1):**
```bash
./prac -o -t 8 1 localhost segmenttreebasic -d 4 -u 3 -q 5
```

**Terminal 3 (Player 2 - Server):**
```bash
./prac -o -t 8 2 localhost localhost segmenttreebasic -d 4 -u 3 -q 5
```

At the end of the run, PRAC will output the number of resources of each type created. This format is suitable for passing as arguments to the preprocessing mode.

### 2. Two-Phase Mode (Preprocessing + Online)

In this approach, you first generate resources using preprocessing mode and then use them in online mode.

**Step 1: Preprocessing**

Run in three terminals on the same host:

**Terminal 1 (Player 0):**
```bash
./prac -t 8 -p 0
```

**Terminal 2 (Player 1):**
```bash
./prac -t 8 -p 1 localhost
```

**Terminal 3 (Player 2):**
```bash
./prac -t 8 -p 2 localhost localhost <Precomputed values string from Online Only Mode output>
```

**Step 2: Online Execution**

After preprocessing is complete, run the online phase:

**Terminal 1 (Player 0):**
```bash
./prac -t 8 0 segmenttreebasic -d 20 -u 10 -q 10
```

**Terminal 2 (Player 1):**
```bash
./prac -t 8 1 localhost segmenttreebasic -d 20 -u 10 -q 10
```

**Terminal 3 (Player 2):**
```bash
./prac -t 8 2 localhost localhost segmenttreebasic -d 20 -u 10 -q 10
```

## Automated Testing (Recommended)

There is a script `run_experiment_tmux.sh` which automates the process using `tmux`. It creates a session with 3 panes for the 3 players.

### Usage

```bash
./run_experiment_tmux.sh <depth> <updates> <queries> <variant> <threads>
```

### Example

```bash
./run_experiment_tmux.sh 20 5 5 segmenttreebasic 8
```

This script automatically:
1.  Runs the experiment.
2.  Saves the results in the `logs` folder.

## Performance Analysis

You can use `analyze_performance.py` to get a summary of the runs and generate visualizations/plots. The results will be saved in the `logs/plots` folder.



### Running the Analysis

```bash
python3 analyze_performance.py logs/performance_*
```
