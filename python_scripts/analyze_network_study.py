#!/usr/bin/env python3

"""
Network Study Analysis Script

This script performs comprehensive analysis of network study experiments.
It has two main parts:

Part A: Per Network Configuration Analysis
  - For each (latency, bandwidth) combination
  - Analyzes onlineonly vs (preprocessing + online) across depths
  - Generates plots similar to analyze_batch_results.py
  - Saves plots to: plots/local/network/<latency>ms_<bandwidth>mbit/

Part B: Network Parameter Impact Analysis
  - Fixed Bandwidth, Varying Latency: Shows how latency affects performance
  - Fixed Latency, Varying Bandwidth: Shows how bandwidth affects performance
  - Saves plots to: plots/local/network/
"""

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import glob
import os
import re
import sys

import matplotlib.pyplot as plt
import seaborn as sns
from rich import box
from rich.console import Console
from rich.table import Table

console = Console()


def parse_filename(filename):
    """
    Extracts parameters from filename.
    Expected format: performance_TYPE_st11_dD_uU_qQ_latLms_bwBmbit_TIMESTAMP.csv
    """
    basename = os.path.basename(filename)

    # Extract log type, depth, updates, queries
    match = re.search(r"performance_([a-z]+)_st11_d(\d+)_u(\d+)_q(\d+)", basename)
    if not match:
        return None, None, None, None, None, None

    log_type = match.group(1)  # onlineonly, preprocessing, online
    depth = int(match.group(2))
    updates = int(match.group(3))
    queries = int(match.group(4))

    # Extract latency and bandwidth
    lat_match = re.search(r"lat(\d+)ms", basename)
    bw_match = re.search(r"bw(\d+)mbit", basename)

    latency = int(lat_match.group(1)) if lat_match else None
    bandwidth = int(bw_match.group(1)) if bw_match else None

    return log_type, depth, updates, queries, latency, bandwidth


def get_avg_update_time(filepath):
    """Extract average update time from performance CSV"""
    try:
        df = pd.read_csv(filepath)
        update_times = df[
            (df["phase"] == "update") & (df["metric"] == "total_update_time")
        ]["value"]
        if not update_times.empty:
            return update_times.mean()
        return 0.0
    except Exception as e:
        console.print(f"[yellow]Warning: Error reading {filepath}: {e}[/yellow]")
        return 0.0


def get_avg_query_time(filepath):
    """Extract average query time from performance CSV"""
    try:
        df = pd.read_csv(filepath)
        query_times = df[
            (df["phase"] == "query") & (df["metric"] == "total_query_time")
        ]["value"]
        if not query_times.empty:
            return query_times.mean()
        return 0.0
    except Exception as e:
        console.print(f"[yellow]Warning: Error reading {filepath}: {e}[/yellow]")
        return 0.0


def get_preprocessing_time(filepath, num_ops):
    """Extract preprocessing time per operation from preprocessing CSV"""
    try:
        df = pd.read_csv(filepath)
        duration_row = df[df["Message"] == "execution_time"]
        if not duration_row.empty:
            duration_sec = float(duration_row["Duration(s)"].iloc[0])
            duration_ms = duration_sec * 1000.0
            if num_ops > 0:
                return duration_ms / num_ops
        return 0.0
    except Exception as e:
        console.print(f"[yellow]Warning: Error reading {filepath}: {e}[/yellow]")
        return 0.0


def load_all_data(log_dir):
    """
    Load all network study data from CSV files.
    Returns nested dict: data[(latency, bandwidth)][depth][config_type][log_type] = filepath
    """
    data = {}
    files = glob.glob(os.path.join(log_dir, "performance_*.csv"))

    for f in files:
        log_type, depth, updates, queries, latency, bandwidth = parse_filename(f)

        # Skip files without network tags
        if log_type is None or latency is None or bandwidth is None:
            continue

        # Initialize nested structure
        net_config = (latency, bandwidth)
        if net_config not in data:
            data[net_config] = {}
        if depth not in data[net_config]:
            data[net_config][depth] = {"update": {}, "query": {}}

        # Categorize by update/query configuration
        if updates > 0 and queries == 0:
            data[net_config][depth]["update"][log_type] = f
        elif updates == 0 and queries > 0:
            data[net_config][depth]["query"][log_type] = f

    return data


def annotate_points(x_vals, y_vals, ax=None):
    """Helper to annotate plot points with values"""
    target = ax if ax else plt
    for x, y in zip(x_vals, y_vals):
        target.annotate(
            f"{y:.1f}",
            (x, y),
            textcoords="offset points",
            xytext=(0, 5),
            ha="center",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.88),
        )


def analyze_single_network_config(net_data, latency, bandwidth, output_dir):
    """
    Part A: Analyze a single network configuration across depths.
    Similar to analyze_batch_results.py but for one network setting.
    """
    console.print(
        f"\n[cyan]Analyzing Network Configuration: Latency={latency}ms, Bandwidth={bandwidth}mbit[/cyan]"
    )

    depths = sorted(net_data.keys())

    # Helper to extract metrics
    def get_metrics_list(phase_key):
        vals_oo = []
        vals_pp = []
        vals_on = []
        vals_total = []
        valid_depths = []

        for d in depths:
            paths = net_data[d].get(phase_key, {})
            if not paths:
                continue

            # OnlineOnly
            t_oo = 0.0
            if "onlineonly" in paths:
                if phase_key == "update":
                    t_oo = get_avg_update_time(paths["onlineonly"])
                else:
                    t_oo = get_avg_query_time(paths["onlineonly"])

            # Preprocessing
            t_pp = 0.0
            if "preprocessing" in paths:
                t_pp = get_preprocessing_time(paths["preprocessing"], 5)

            # Online
            t_o = 0.0
            if "online" in paths:
                if phase_key == "update":
                    t_o = get_avg_update_time(paths["online"])
                else:
                    t_o = get_avg_query_time(paths["online"])

            vals_oo.append(t_oo)
            vals_pp.append(t_pp)
            vals_on.append(t_o)
            vals_total.append(t_pp + t_o)
            valid_depths.append(d)

        return valid_depths, vals_oo, vals_pp, vals_on, vals_total

    # Create plots directory for this network config
    config_dir = os.path.join(output_dir, f"{latency}ms_{bandwidth}mbit")
    os.makedirs(config_dir, exist_ok=True)

    # 1. Update Breakdown
    d_u, u_oo, u_pp, u_on, u_tot = get_metrics_list("update")
    if d_u:
        plt.figure(figsize=(10, 6))
        plt.plot(d_u, u_oo, marker="o", label="OnlineOnly (Single Phase)")
        annotate_points(d_u, u_oo)
        plt.plot(d_u, u_pp, marker="s", label="Preprocessing")
        annotate_points(d_u, u_pp)
        plt.plot(d_u, u_on, marker="^", label="Online")
        annotate_points(d_u, u_on)
        plt.plot(d_u, u_tot, marker="D", label="Total (Preproc+Online)")
        annotate_points(d_u, u_tot)
        plt.xlabel("Depth")
        plt.ylabel("Time (ms)")
        plt.title(f"Update Performance - Lat={latency}ms, BW={bandwidth}mbit")
        plt.grid(True)
        plt.legend()
        plt.savefig(os.path.join(config_dir, "update_breakdown.png"))
        plt.close()

    # 2. Query Breakdown
    d_q, q_oo, q_pp, q_on, q_tot = get_metrics_list("query")
    if d_q:
        plt.figure(figsize=(10, 6))
        plt.plot(d_q, q_oo, marker="o", label="OnlineOnly (Single Phase)")
        annotate_points(d_q, q_oo)
        plt.plot(d_q, q_pp, marker="s", label="Preprocessing")
        annotate_points(d_q, q_pp)
        plt.plot(d_q, q_on, marker="^", label="Online")
        annotate_points(d_q, q_on)
        plt.plot(d_q, q_tot, marker="D", label="Total (Preproc+Online)")
        annotate_points(d_q, q_tot)
        plt.xlabel("Depth")
        plt.ylabel("Time (ms)")
        plt.title(f"Query Performance - Lat={latency}ms, BW={bandwidth}mbit")
        plt.grid(True)
        plt.legend()
        plt.savefig(os.path.join(config_dir, "query_breakdown.png"))
        plt.close()

    # 3. Update Comparison
    if d_u:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        ax1.plot(d_u, u_oo, marker="o", color="blue")
        annotate_points(d_u, u_oo, ax=ax1)
        ax1.set_xlabel("Depth")
        ax1.set_ylabel("Time (ms)")
        ax1.set_title("Update: Single Phase (OnlineOnly)")
        ax1.grid(True)

        ax2.plot(d_u, u_tot, marker="D", color="red")
        annotate_points(d_u, u_tot, ax=ax2)
        ax2.set_xlabel("Depth")
        ax2.set_ylabel("Time (ms)")
        ax2.set_title("Update: Two Phase (Preproc + Online)")
        ax2.grid(True)

        plt.suptitle(f"Latency={latency}ms, Bandwidth={bandwidth}mbit")
        plt.tight_layout()
        plt.savefig(os.path.join(config_dir, "update_comparison.png"))
        plt.close()

    # 4. Query Comparison
    if d_q:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        ax1.plot(d_q, q_oo, marker="o", color="blue")
        annotate_points(d_q, q_oo, ax=ax1)
        ax1.set_xlabel("Depth")
        ax1.set_ylabel("Time (ms)")
        ax1.set_title("Query: Single Phase (OnlineOnly)")
        ax1.grid(True)

        ax2.plot(d_q, q_tot, marker="D", color="red")
        annotate_points(d_q, q_tot, ax=ax2)
        ax2.set_xlabel("Depth")
        ax2.set_ylabel("Time (ms)")
        ax2.set_title("Query: Two Phase (Preproc + Online)")
        ax2.grid(True)

        plt.suptitle(f"Latency={latency}ms, Bandwidth={bandwidth}mbit")
        plt.tight_layout()
        plt.savefig(os.path.join(config_dir, "query_comparison.png"))
        plt.close()

    console.print(f"[green]Plots saved to: {config_dir}/[/green]")


def analyze_network_impact(data, output_dir):
    """
    Part B: Analyze how network parameters affect performance.
    Creates plots showing impact of latency and bandwidth.
    """
    console.print("\n[cyan]Analyzing Network Parameter Impact...[/cyan]")

    # Collect all unique latencies, bandwidths, and depths
    latencies = sorted(set(lat for lat, bw in data.keys()))
    bandwidths = sorted(set(bw for lat, bw in data.keys()))

    # For a reference depth (choose middle depth if available)
    all_depths = set()
    for net_config in data.values():
        all_depths.update(net_config.keys())
    all_depths = sorted(all_depths)

    if not all_depths:
        console.print(
            "[yellow]No depth data found for network impact analysis[/yellow]"
        )
        return

    # Use middle depth as reference
    ref_depth = all_depths[len(all_depths) // 2]
    console.print(f"Using reference depth: {ref_depth}")

    # 1. Fixed Bandwidth, Varying Latency
    for bw in bandwidths:
        # Collect data for this bandwidth across all latencies
        lat_vals = []
        oo_update = []
        pp_update = []
        on_update = []
        oo_query = []
        pp_query = []
        on_query = []

        for lat in latencies:
            net_config = (lat, bw)
            if net_config not in data or ref_depth not in data[net_config]:
                continue

            lat_vals.append(lat)

            # Update metrics
            upd_paths = data[net_config][ref_depth].get("update", {})
            if "onlineonly" in upd_paths:
                oo_update.append(get_avg_update_time(upd_paths["onlineonly"]))
            else:
                oo_update.append(0)
            if "preprocessing" in upd_paths:
                pp_update.append(get_preprocessing_time(upd_paths["preprocessing"], 5))
            else:
                pp_update.append(0)
            if "online" in upd_paths:
                on_update.append(get_avg_update_time(upd_paths["online"]))
            else:
                on_update.append(0)

            # Query metrics
            qry_paths = data[net_config][ref_depth].get("query", {})
            if "onlineonly" in qry_paths:
                oo_query.append(get_avg_query_time(qry_paths["onlineonly"]))
            else:
                oo_query.append(0)
            if "preprocessing" in qry_paths:
                pp_query.append(get_preprocessing_time(qry_paths["preprocessing"], 5))
            else:
                pp_query.append(0)
            if "online" in qry_paths:
                on_query.append(get_avg_query_time(qry_paths["online"]))
            else:
                on_query.append(0)

        if not lat_vals:
            continue

        # Calculate two-phase totals
        two_phase_update = [pp + on for pp, on in zip(pp_update, on_update)]
        two_phase_query = [pp + on for pp, on in zip(pp_query, on_query)]

        # Plot: Update performance vs latency (separate file)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Single Phase - OnlineOnly
        ax1.plot(
            lat_vals, oo_update, marker="o", color="blue", label="OnlineOnly (Time)"
        )
        annotate_points(lat_vals, oo_update, ax=ax1)
        ax1.set_xlabel("Latency (ms)")
        ax1.set_ylabel("Time (ms)")
        ax1.set_title("Single Phase (OnlineOnly)")
        ax1.grid(True)
        ax1.legend()

        # Two Phase - Preprocessing + Online + Total
        ax2.plot(lat_vals, pp_update, marker="s", color="orange", label="Preprocessing")
        annotate_points(lat_vals, pp_update, ax=ax2)
        ax2.plot(lat_vals, on_update, marker="^", color="green", label="Online")
        annotate_points(lat_vals, on_update, ax=ax2)
        ax2.plot(
            lat_vals,
            two_phase_update,
            marker="D",
            color="red",
            label="Preproc + Online",
        )
        annotate_points(lat_vals, two_phase_update, ax=ax2)
        ax2.set_xlabel("Latency (ms)")
        ax2.set_ylabel("Time (ms)")
        ax2.set_title("Two Phase (Preproc + Online)")
        ax2.grid(True)
        ax2.legend()

        plt.suptitle(f"Update Performance vs Latency (BW={bw}mbit, Depth={ref_depth})")
        plt.tight_layout()
        plt.savefig(
            os.path.join(
                output_dir, f"latency_impact_update_bw{bw}mbit_d{ref_depth}.png"
            )
        )
        plt.close()
        console.print(
            f"[green]Saved: latency_impact_update_bw{bw}mbit_d{ref_depth}.png[/green]"
        )

        # Plot: Query performance vs latency (separate file)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Single Phase - OnlineOnly
        ax1.plot(
            lat_vals, oo_query, marker="o", color="blue", label="OnlineOnly (Time)"
        )
        annotate_points(lat_vals, oo_query, ax=ax1)
        ax1.set_xlabel("Latency (ms)")
        ax1.set_ylabel("Time (ms)")
        ax1.set_title("Single Phase (OnlineOnly)")
        ax1.grid(True)
        ax1.legend()

        # Two Phase - Preprocessing + Online + Total
        ax2.plot(lat_vals, pp_query, marker="s", color="orange", label="Preprocessing")
        annotate_points(lat_vals, pp_query, ax=ax2)
        ax2.plot(lat_vals, on_query, marker="^", color="green", label="Online")
        annotate_points(lat_vals, on_query, ax=ax2)
        ax2.plot(
            lat_vals, two_phase_query, marker="D", color="red", label="Preproc + Online"
        )
        annotate_points(lat_vals, two_phase_query, ax=ax2)
        ax2.set_xlabel("Latency (ms)")
        ax2.set_ylabel("Time (ms)")
        ax2.set_title("Two Phase (Preproc + Online)")
        ax2.grid(True)
        ax2.legend()

        plt.suptitle(f"Query Performance vs Latency (BW={bw}mbit, Depth={ref_depth})")
        plt.tight_layout()
        plt.savefig(
            os.path.join(
                output_dir, f"latency_impact_query_bw{bw}mbit_d{ref_depth}.png"
            )
        )
        plt.close()
        console.print(
            f"[green]Saved: latency_impact_query_bw{bw}mbit_d{ref_depth}.png[/green]"
        )

    # 2. Fixed Latency, Varying Bandwidth
    for lat in latencies:
        # Collect data for this latency across all bandwidths
        bw_vals = []
        oo_update = []
        pp_update = []
        on_update = []
        oo_query = []
        pp_query = []
        on_query = []

        for bw in bandwidths:
            net_config = (lat, bw)
            if net_config not in data or ref_depth not in data[net_config]:
                continue

            bw_vals.append(bw)

            # Update metrics
            upd_paths = data[net_config][ref_depth].get("update", {})
            if "onlineonly" in upd_paths:
                oo_update.append(get_avg_update_time(upd_paths["onlineonly"]))
            else:
                oo_update.append(0)
            if "preprocessing" in upd_paths:
                pp_update.append(get_preprocessing_time(upd_paths["preprocessing"], 5))
            else:
                pp_update.append(0)
            if "online" in upd_paths:
                on_update.append(get_avg_update_time(upd_paths["online"]))
            else:
                on_update.append(0)

            # Query metrics
            qry_paths = data[net_config][ref_depth].get("query", {})
            if "onlineonly" in qry_paths:
                oo_query.append(get_avg_query_time(qry_paths["onlineonly"]))
            else:
                oo_query.append(0)
            if "preprocessing" in qry_paths:
                pp_query.append(get_preprocessing_time(qry_paths["preprocessing"], 5))
            else:
                pp_query.append(0)
            if "online" in qry_paths:
                on_query.append(get_avg_query_time(qry_paths["online"]))
            else:
                on_query.append(0)

        if not bw_vals:
            continue

        # Calculate two-phase totals
        two_phase_update = [pp + on for pp, on in zip(pp_update, on_update)]
        two_phase_query = [pp + on for pp, on in zip(pp_query, on_query)]

        # Plot: Update performance vs bandwidth (separate file)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Single Phase - OnlineOnly
        ax1.plot(
            bw_vals, oo_update, marker="o", color="blue", label="OnlineOnly (Time)"
        )
        annotate_points(bw_vals, oo_update, ax=ax1)
        ax1.set_xlabel("Bandwidth (mbit)")
        ax1.set_ylabel("Time (ms)")
        ax1.set_title("Single Phase (OnlineOnly)")
        ax1.grid(True)
        ax1.legend()

        # Two Phase - Preprocessing + Online + Total
        ax2.plot(bw_vals, pp_update, marker="s", color="orange", label="Preprocessing")
        annotate_points(bw_vals, pp_update, ax=ax2)
        ax2.plot(bw_vals, on_update, marker="^", color="green", label="Online")
        annotate_points(bw_vals, on_update, ax=ax2)
        ax2.plot(
            bw_vals, two_phase_update, marker="D", color="red", label="Preproc + Online"
        )
        annotate_points(bw_vals, two_phase_update, ax=ax2)
        ax2.set_xlabel("Bandwidth (mbit)")
        ax2.set_ylabel("Time (ms)")
        ax2.set_title("Two Phase (Preproc + Online)")
        ax2.grid(True)
        ax2.legend()

        plt.suptitle(
            f"Update Performance vs Bandwidth (Lat={lat}ms, Depth={ref_depth})"
        )
        plt.tight_layout()
        plt.savefig(
            os.path.join(
                output_dir, f"bandwidth_impact_update_lat{lat}ms_d{ref_depth}.png"
            )
        )
        plt.close()
        console.print(
            f"[green]Saved: bandwidth_impact_update_lat{lat}ms_d{ref_depth}.png[/green]"
        )

        # Plot: Query performance vs bandwidth (separate file)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Single Phase - OnlineOnly
        ax1.plot(bw_vals, oo_query, marker="o", color="blue", label="OnlineOnly (Time)")
        annotate_points(bw_vals, oo_query, ax=ax1)
        ax1.set_xlabel("Bandwidth (mbit)")
        ax1.set_ylabel("Time (ms)")
        ax1.set_title("Single Phase (OnlineOnly)")
        ax1.grid(True)
        ax1.legend()

        # Two Phase - Preprocessing + Online + Total
        ax2.plot(bw_vals, pp_query, marker="s", color="orange", label="Preprocessing")
        annotate_points(bw_vals, pp_query, ax=ax2)
        ax2.plot(bw_vals, on_query, marker="^", color="green", label="Online")
        annotate_points(bw_vals, on_query, ax=ax2)
        ax2.plot(
            bw_vals, two_phase_query, marker="D", color="red", label="Preproc + Online"
        )
        annotate_points(bw_vals, two_phase_query, ax=ax2)
        ax2.set_xlabel("Bandwidth (mbit)")
        ax2.set_ylabel("Time (ms)")
        ax2.set_title("Two Phase (Preproc + Online)")
        ax2.grid(True)
        ax2.legend()

        plt.suptitle(f"Query Performance vs Bandwidth (Lat={lat}ms, Depth={ref_depth})")
        plt.tight_layout()
        plt.savefig(
            os.path.join(
                output_dir, f"bandwidth_impact_query_lat{lat}ms_d{ref_depth}.png"
            )
        )
        plt.close()
        console.print(
            f"[green]Saved: bandwidth_impact_query_lat{lat}ms_d{ref_depth}.png[/green]"
        )


def main():
    if len(sys.argv) < 2:
        console.print("[red]Usage: python3 analyze_network_study.py <logs_dir>[/red]")
        console.print("Example: python3 analyze_network_study.py logs/")
        sys.exit(1)

    log_dir = sys.argv[1]

    console.print("[bold cyan]=========================================[/bold cyan]")
    console.print("[bold cyan]Network Study Analysis[/bold cyan]")
    console.print("[bold cyan]=========================================[/bold cyan]")

    # Load all data
    console.print(f"\n[cyan]Loading data from: {log_dir}[/cyan]")
    data = load_all_data(log_dir)

    if not data:
        console.print("[red]No network study data found![/red]")
        console.print(
            "[yellow]Looking for files matching: performance_*_lat*ms_bw*mbit*.csv[/yellow]"
        )
        sys.exit(1)

    console.print(f"[green]Found {len(data)} network configurations[/green]")

    # Create output directories
    network_plots_dir = "plots/local/network"
    os.makedirs(network_plots_dir, exist_ok=True)

    # Part A: Analyze each network configuration
    console.print("\n[bold]Part A: Per-Configuration Analysis[/bold]")
    for (latency, bandwidth), net_data in sorted(data.items()):
        analyze_single_network_config(net_data, latency, bandwidth, network_plots_dir)

    # Part B: Network parameter impact analysis
    console.print("\n[bold]Part B: Network Parameter Impact Analysis[/bold]")
    analyze_network_impact(data, network_plots_dir)

    console.print(
        "\n[bold green]=========================================[/bold green]"
    )
    console.print("[bold green]Analysis Complete![/bold green]")
    console.print("[bold green]=========================================[/bold green]")
    console.print(f"\nPlots saved to: {network_plots_dir}/")
    console.print("  - Per-config plots: <latency>ms_<bandwidth>mbit/")
    console.print("  - Impact analysis:")
    console.print("    * latency_impact_update_*.png")
    console.print("    * latency_impact_query_*.png")
    console.print("    * bandwidth_impact_update_*.png")
    console.print("    * bandwidth_impact_query_*.png")


if __name__ == "__main__":
    main()
