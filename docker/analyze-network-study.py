#!/usr/bin/env python3
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import re
import sys
from rich.console import Console
from rich.table import Table

console = Console()

def parse_network_params(filename):
    """
    Extracts latency and bandwidth from filename.
    Expected format: ..._lat{X}ms_bw{Y}mbit...csv
    """
    lat_match = re.search(r'lat(\d+)ms', filename)
    bw_match = re.search(r'bw(\d+)mbit', filename)
    
    lat = int(lat_match.group(1)) if lat_match else 0
    bw = int(bw_match.group(1)) if bw_match else 0
    return lat, bw

def load_data(log_dir):
    all_files = glob.glob(os.path.join(log_dir, "*.csv"))
    data = []
    
    for f in all_files:
        try:
            # Check if this file has network params
            if "lat" not in f or "bw" not in f:
                continue
                
            lat, bw = parse_network_params(f)
            df = pd.read_csv(f)
            
            # Extract key metrics
            # 1. Total Update Time
            update_time = df[(df['phase'] == 'update') & 
                           (df['metric'] == 'total_update_time')]['value'].mean()
            
            # 2. Total Query Time
            query_time = df[(df['phase'] == 'query') & 
                          (df['metric'] == 'total_query_time')]['value'].mean()
            
            # 3. Total Wall Clock (execution time) 
            # In segmentTree9.cpp, it logs 'total_all_updates_time' and 'total_all_queries_time'
            
            # Helper to get value safely
            def get_metric(phase, metric, default=0):
                rows = df[(df['phase'] == phase) & (df['metric'] == metric)]
                if not rows.empty:
                    return rows['value'].iloc[0]
                return default

            wall_update = get_metric('update_summary', 'total_all_updates_time')
            wall_query = get_metric('query_summary', 'total_all_queries_time')
            
            # Fallback: if summaries missing, sum up individual ops
            if wall_update == 0:
                 wall_update = df[(df['phase'] == 'update') & (df['metric'] == 'total_update_time')]['value'].sum()
            if wall_query == 0:
                 wall_query = df[(df['phase'] == 'query') & (df['metric'] == 'total_query_time')]['value'].sum()

            total_time = wall_update + wall_query
            
            data.append({
                'latency_ms': lat,
                'bandwidth_mbit': bw,
                'avg_update_time_ms': update_time,
                'avg_query_time_ms': query_time,
                'total_execution_time_ms': total_time,
                'filename': f
            })
        except Exception as e:
            console.print(f"[yellow]Skipping {f}: {e}[/yellow]")
            
    return pd.DataFrame(data)

def plot_metric_vs_param(df, x_param, y_param, fixed_param, fixed_value, title, ylabel, filename):
    """Helper to plot a single metric"""
    plt.figure(figsize=(10, 6))
    subset = df[df[fixed_param] == fixed_value].sort_values(x_param)
    
    if len(subset) < 2:
        return

    sns.lineplot(data=subset, x=x_param, y=y_param, marker='o', linewidth=2)
    
    # Add labels
    for _, row in subset.iterrows():
        plt.annotate(f"{row[y_param]:.1f}", 
                     (row[x_param], row[y_param]),
                     textcoords="offset points", xytext=(0,10), ha='center')

    plt.title(title)
    plt.xlabel('Latency (ms)' if x_param == 'latency_ms' else 'Bandwidth (Mbit)')
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.3)
    
    plt.savefig(filename, dpi=300)
    console.print(f"[green]Saved {filename}[/green]")
    plt.close()

def plot_latency_study(df, output_dir):
    """Plots metrics vs Latency for each fixed bandwidth"""
    bws = df['bandwidth_mbit'].unique()
    
    for bw in bws:
        # Plot 1: Avg Update Time
        plot_metric_vs_param(
            df, 'latency_ms', 'avg_update_time_ms', 
            'bandwidth_mbit', bw,
            f'Average Update Time vs Latency (Bandwidth: {bw} Mbit)',
            'Average Update Time (ms)',
            f"{output_dir}/study_latency_bw{bw}_update.png"
        )
        
        # Plot 2: Avg Query Time
        plot_metric_vs_param(
            df, 'latency_ms', 'avg_query_time_ms', 
            'bandwidth_mbit', bw,
            f'Average Query Time vs Latency (Bandwidth: {bw} Mbit)',
            'Average Query Time (ms)',
            f"{output_dir}/study_latency_bw{bw}_query.png"
        )

def plot_bandwidth_study(df, output_dir):
    """Plots metrics vs Bandwidth for each fixed latency"""
    lats = df['latency_ms'].unique()
    
    for lat in lats:
        # Plot 1: Avg Update Time
        plot_metric_vs_param(
            df, 'bandwidth_mbit', 'avg_update_time_ms', 
            'latency_ms', lat,
            f'Average Update Time vs Bandwidth (Latency: {lat} ms)',
            'Average Update Time (ms)',
            f"{output_dir}/study_bandwidth_lat{lat}_update.png"
        )
        
        # Plot 2: Avg Query Time
        plot_metric_vs_param(
            df, 'bandwidth_mbit', 'avg_query_time_ms', 
            'latency_ms', lat,
            f'Average Query Time vs Bandwidth (Latency: {lat} ms)',
            'Average Query Time (ms)',
            f"{output_dir}/study_bandwidth_lat{lat}_query.png"
        )

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 analyze-network-study.py <logs_dir>")
        sys.exit(1)
        
    log_dir = sys.argv[1]
    output_dir = os.path.join(log_dir, "plots")
    os.makedirs(output_dir, exist_ok=True)
    
    console.print(f"[bold cyan]Analyzing Network Study Logs in {log_dir}...[/bold cyan]")
    
    df = load_data(log_dir)
    
    if df.empty:
        console.print("[red]No tagged log files found (looking for *_latXms_bwYmbit*.csv)[/red]")
        sys.exit(1)
        
    console.print("\n[bold]Summary Data:[/bold]")
    console.print(df[['latency_ms', 'bandwidth_mbit', 'avg_update_time_ms', 'avg_query_time_ms']].sort_values(['latency_ms', 'bandwidth_mbit']).to_string())
    
    plot_latency_study(df, output_dir)
    plot_bandwidth_study(df, output_dir)

if __name__ == "__main__":
    main()
