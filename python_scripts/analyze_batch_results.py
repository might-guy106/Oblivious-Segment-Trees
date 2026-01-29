import glob
import os
import re

import pandas as pd
from rich import box
from rich.console import Console
from rich.table import Table

console = Console()

LOG_DIR = "logs"


def parse_filename(filename):
    basename = os.path.basename(filename)
    # Regex to extract type and parameters
    # Matches: performance_TYPE_st9_d(\d+)_u(\d+)_q(\d+)_...
    match = re.search(r"performance_([a-z]+)_st11_d(\d+)_u(\d+)_q(\d+)_", basename)
    if match:
        log_type = match.group(1)  # online, onlineonly, preprocessing
        depth = int(match.group(2))
        updates = int(match.group(3))
        queries = int(match.group(4))
        return log_type, depth, updates, queries
    return None, None, None, None


def get_avg_update_time(filepath):
    try:
        df = pd.read_csv(filepath)
        # Filter for total_update_time per operation
        update_times = df[
            (df["phase"] == "update") & (df["metric"] == "total_update_time")
        ]["value"]
        if not update_times.empty:
            return update_times.mean()
        return 0.0
    except Exception as e:
        console.print(f"[red]Error reading {filepath}: {e}[/red]")
        return 0.0


def get_avg_query_time(filepath):
    try:
        df = pd.read_csv(filepath)
        # Filter for total_query_time per operation
        # Note: In the viewed file, it seemed query metrics might be similar.
        # Accessing checking the analyze_performance.py, query metric is 'total_query_time'
        query_times = df[
            (df["phase"] == "query") & (df["metric"] == "total_query_time")
        ]["value"]
        if not query_times.empty:
            return query_times.mean()
        return 0.0
    except Exception as e:
        console.print(f"[red]Error reading {filepath}: {e}[/red]")
        return 0.0


def get_preprocessing_time(filepath, num_ops):
    try:
        df = pd.read_csv(filepath)
        # Format: Timestamp,Message,Duration(s)
        duration_row = df[df["Message"] == "execution_time"]
        if not duration_row.empty:
            duration_sec = float(duration_row["Duration(s)"].iloc[0])
            duration_ms = duration_sec * 1000.0
            if num_ops > 0:
                return duration_ms / num_ops
        return 0.0
    except Exception as e:
        console.print(f"[red]Error reading {filepath}: {e}[/red]")
        return 0.0


def analyze_results():
    # Data structure: data[depth][config_type] = { 'onlineonly': path, 'preprocessing': path, 'online': path }
    # config_type is 'update' (u=5, q=0) or 'query' (u=0, q=5)
    data = {}

    files = glob.glob(os.path.join(LOG_DIR, "performance_*.csv"))

    for f in files:
        log_type, depth, updates, queries = parse_filename(f)
        if not log_type:
            continue

        if depth not in data:
            data[depth] = {"update": {}, "query": {}}

        if updates > 0 and queries == 0:
            data[depth]["update"][log_type] = f
        elif updates == 0 and queries > 0:
            data[depth]["query"][log_type] = f

    # Create Tables
    depths = sorted(data.keys())

    # Update Performance Table
    upd_table = Table(title="Update Performance Analysis (ms)", box=box.SQUARE)
    upd_table.add_column("Depth", style="cyan", justify="right")
    upd_table.add_column("OnlineOnly (Avg)", style="green", justify="right")
    upd_table.add_column("Preproc (Avg)", style="yellow", justify="right")
    upd_table.add_column("Online (Avg)", style="blue", justify="right")
    upd_table.add_column("Single-Phase", style="magenta", justify="right")
    upd_table.add_column("Two-Phase", style="red", justify="right")

    for d in depths:
        paths = data[d]["update"]
        if not paths:
            continue

        # OnlineOnly
        t_oo = 0.0
        if "onlineonly" in paths:
            t_oo = get_avg_update_time(paths["onlineonly"])

        # Preprocessing
        t_pp = 0.0
        # Determine num_updates from filename or assume 5 based on script
        # Parsing filename again or passing it creates complexity.
        # Easier: open the file and read num_updates if needed, or rely on our convention.
        # Convention: We know for this batch script updates=5.
        if "preprocessing" in paths:
            t_pp = get_preprocessing_time(paths["preprocessing"], 5)

        # Online
        t_o = 0.0
        if "online" in paths:
            t_o = get_avg_update_time(paths["online"])

        single_phase = t_oo
        two_phase = t_pp + t_o

        upd_table.add_row(
            str(d),
            f"{t_oo:.2f}",
            f"{t_pp:.2f}",
            f"{t_o:.2f}",
            f"{single_phase:.2f}",
            f"{two_phase:.2f}",
        )

    console.print(upd_table)
    console.print()

    # Query Performance Table
    qry_table = Table(title="Query Performance Analysis (ms)", box=box.SQUARE)
    qry_table.add_column("Depth", style="cyan", justify="right")
    qry_table.add_column("OnlineOnly (Avg)", style="green", justify="right")
    qry_table.add_column("Preproc (Avg)", style="yellow", justify="right")
    qry_table.add_column("Online (Avg)", style="blue", justify="right")
    qry_table.add_column("Single-Phase", style="magenta", justify="right")
    qry_table.add_column("Two-Phase", style="red", justify="right")

    for d in depths:
        paths = data[d]["query"]
        if not paths:
            continue

        # OnlineOnly
        t_oo = 0.0
        if "onlineonly" in paths:
            t_oo = get_avg_query_time(paths["onlineonly"])

        # Preprocessing
        t_pp = 0.0
        # Convention: We know for this batch script queries=5.
        if "preprocessing" in paths:
            t_pp = get_preprocessing_time(paths["preprocessing"], 5)

        # Online
        t_o = 0.0
        if "online" in paths:
            t_o = get_avg_query_time(paths["online"])

        single_phase = t_oo
        two_phase = t_pp + t_o

        qry_table.add_row(
            str(d),
            f"{t_oo:.2f}",
            f"{t_pp:.2f}",
            f"{t_o:.2f}",
            f"{single_phase:.2f}",
            f"{two_phase:.2f}",
        )

    console.print(qry_table)

    # Plotting
    plot_results(data)


def plot_results(data):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        console.print("[red]matplotlib not found. Skipping plots.[/red]")
        return

    depths = sorted(data.keys())
    plots_dir = "plots/local/default"
    os.makedirs(plots_dir, exist_ok=True)

    # helper to extract metric list
    def get_metrics_list(metric_type, phase_key):
        vals_oo = []
        vals_pp = []
        vals_on = []
        vals_total = []
        valid_depths = []

        for d in depths:
            paths = data[d][phase_key]
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
            # Convention: updates=5, queries=5
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

    def annotate_points(x_vals, y_vals, ax=None):
        target = ax if ax else plt
        for x, y in zip(x_vals, y_vals):
            target.annotate(
                f"{y:.1f}",
                (x, y),
                textcoords="offset points",
                xytext=(0, 5),
                ha="center",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9),
            )

    # 1. Update Breakdown
    d_u, u_oo, u_pp, u_on, u_tot = get_metrics_list("update", "update")
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
        plt.title("Update Performance Breakdown")
        plt.grid(True)
        plt.legend()
        plt.savefig(os.path.join(plots_dir, "update_breakdown.png"))
        plt.close()

    # 2. Query Breakdown
    d_q, q_oo, q_pp, q_on, q_tot = get_metrics_list("query", "query")
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
        plt.title("Query Performance Breakdown")
        plt.grid(True)
        plt.legend()
        plt.savefig(os.path.join(plots_dir, "query_breakdown.png"))
        plt.close()

    # 3. Update Comparison (Side-by-Side)
    if d_u:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Single Phase
        ax1.plot(d_u, u_oo, marker="o", color="blue")
        annotate_points(d_u, u_oo, ax=ax1)
        ax1.set_xlabel("Depth")
        ax1.set_ylabel("Time (ms)")
        ax1.set_title("Update: Single Phase (OnlineOnly)")
        ax1.grid(True)

        # Two Phase
        ax2.plot(d_u, u_tot, marker="D", color="red")
        annotate_points(d_u, u_tot, ax=ax2)
        ax2.set_xlabel("Depth")
        ax2.set_ylabel("Time (ms)")
        ax2.set_title("Update: Two Phase (Preproc + Online)")
        ax2.grid(True)

        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "update_comparison.png"))
        plt.close()

    # 4. Query Comparison (Side-by-Side)
    if d_q:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Single Phase
        ax1.plot(d_q, q_oo, marker="o", color="blue")
        annotate_points(d_q, q_oo, ax=ax1)
        ax1.set_xlabel("Depth")
        ax1.set_ylabel("Time (ms)")
        ax1.set_title("Query: Single Phase (OnlineOnly)")
        ax1.grid(True)

        # Two Phase
        ax2.plot(d_q, q_tot, marker="D", color="red")
        annotate_points(d_q, q_tot, ax=ax2)
        ax2.set_xlabel("Depth")
        ax2.set_ylabel("Time (ms)")
        ax2.set_title("Query: Two Phase (Preproc + Online)")
        ax2.grid(True)

        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "query_comparison.png"))
        plt.close()

    console.print(
        "[green]Plots saved: update_breakdown.png, query_breakdown.png, update_comparison.png, query_comparison.png[/green]"
    )


if __name__ == "__main__":
    analyze_results()
