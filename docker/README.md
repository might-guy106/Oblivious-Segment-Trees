# Docker: Quickstart (3-Party PRAC + Oblivious Segment Trees)

This folder contains a Docker setup that runs the PRAC framework and Oblivious Segment Tree experiments as **three parties** (`prac_p0`, `prac_p1`, `prac_p2`) on a **single machine**.

At a high level:

- `Dockerfile` builds an Ubuntu image, installs dependencies, **copies the repo code**, and runs `make -j` during image build.
- `start-docker` launches **3 containers** (one per party) on a shared Docker network.
- `run-single-experiment` runs one experiment configuration in `onlineonly`, `preprocessing`, `online`, or `all` (and handles party IDs/IPs automatically).
- `run-batch-experiments` runs a sweep across multiple depths (runs both update and query configurations per depth).
- `run-network-study` runs batch experiments across sweeps of latency/bandwidth settings (via `set-networking`).
- `set-networking` / `unset-networking` can inject latency/bandwidth limits using `tc`.

---

## Prerequisites

- Docker installed and working (`docker ps` succeeds)
- On Linux, you may need to run scripts with `sudo` depending on your Docker setup and traffic control permissions.

---

## Quickstart

### 1) Build the image
You typically rebuild when C++ source changes (because the code is baked into the image during build).

```bash
cd docker
./build-docker
```

### 2) Start the 3 containers (parties)
```bash
./start-docker
```

This starts three containers representing the parties:
- `prac_p0`
- `prac_p1`
- `prac_p2`

### 3) Run a simple segment tree experiment
Use `run-single-experiment` (it runs inside all three containers and handles party IDs/IPs for you).

Example (online-only mode):
```bash
# Usage: ./run-single-experiment <mode> <depth> <updates> <queries> [variant] [threads]
./run-single-experiment onlineonly 4 3 5 segmenttreebasic
```

Run the full pipeline (onlineonly → preprocessing → online):
```bash
./run-single-experiment all 4 3 5 segmenttreebasic
```

---

## Common workflows

### A) Single experiment
Use `run-single-experiment` to run one configuration in one mode or all modes.

```bash
# Run all phases (onlineonly, preprocessing, online)
./run-single-experiment all <depth> <updates> <queries> [variant] [threads]

# Example:
./run-single-experiment all 8 5 0 segmenttreebasic
```

<!--Run individual phases:
```bash
./run-single-experiment onlineonly <depth> <updates> <queries> [variant] [threads]
./run-single-experiment preprocessing <depth> <updates> <queries> [variant] [threads]
./run-single-experiment online <depth> <updates> <queries> [variant] [threads]
```-->

### B) Batch experiments (multiple depths/configs)
Runs, for each depth, both:
- update experiment: `updates=<UPDATES>, queries=0`
- query experiment: `updates=0, queries=<QUERIES>`

```bash
./run-batch-experiments --depths "4 6 8" --updates "5" --queries "5" --variant segmenttreebasic
```

Optional threads control:
```bash
./run-batch-experiments --depths "4 6 8" --variant segmenttreebasic --threads 8
# or auto:
./run-batch-experiments --depths "4 6 8" --variant segmenttreebasic --threads auto
```

### C) Network study (vary latency/bandwidth)
Sweeps `--latencies` × `--bandwidths`. For each network configuration it:
1) cleans logs, 2) sets networking, 3) runs batch experiments, 4) tags logs, 5) unsets networking.

```bash
./run-network-study --latencies "0 20 40" --bandwidths "10 50 100" --depths "4 6 8" --variant segmenttreebasic
```

<!-----

## Network simulation (optional)

Set artificial latency and bandwidth limits:
```bash
./set-networking 30ms 100mbit
```

Reset to normal networking:
```bash
./unset-networking
```-->

---

## Logs analysis

### Bring results to your host (if available in this repo)
```bash
./bring-plots
```
- Logs copied to `../logs/`

### Local analysis (run on host)
From the repo root:
```bash
python3 python_scripts/analyze_batch_results.py
python3 python_scripts/analyze_network_study.py logs/
```

Plot locations :
- Batch experiments: `plots/local/default/`
- Network study (per config): `plots/local/network/<latency>ms_<bandwidth>mbit/`
- Network study (impact): `plots/local/network/`
---

## Cleanup and stopping

### Clean logs inside containers
```bash
./clean-logs
```

Stop and remove the containers :
```bash
./stop-docker
```

---





<!--## Notes & troubleshooting

- **Source changes require rebuild**: If you edit C++ sources, re-run:
  ```bash
  ./build-docker
  ```
- **Script permissions**: Ensure scripts are executable:
  ```bash
  chmod +x *
  ```
- **Docker permissions**: If you get permission errors, try running the scripts with `sudo` (common when using `tc` for networking).-->

## Simple Workflow

Build docker image
```bash
./build-docker
```
Start containers
```bash
./start-docker
```
Run batch experiments (below uses default parameters)
```bash
./run-batch-experiments
```
Run network study (below uses default parameters)
```bash
./run-network-study
```
Bring logs outside of docker containers
```bash
./bring-plots
```
once we bring logs outside of docker containers, we can analyze them locally using:
```bash
python3 python_scripts/analyze_batch_results.py
python3 python_scripts/analyze_network_study.py logs/
```
Clean the logs (optional)
```bash
./clean-logs
```
Stop the containers (optional)
```bash
./stop-docker
```
