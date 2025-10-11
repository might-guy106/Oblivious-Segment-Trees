#!/usr/bin/env python3
"""
Performance Analysis Script for Segment Tree Experiments

Usage:
    python3 analyze_performance.py logs/performance_exp001.csv
    python3 analyze_performance.py logs/performance_*.csv  # Compare multiple experiments
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import glob
import os
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
from rich.progress import track
from rich import box

# Initialize rich console
console = Console()

def load_performance_data(csv_path):
    """Load and return performance data from CSV"""
    df = pd.read_csv(csv_path)
    return df

def analyze_single_experiment(df):
    """Analyze a single experiment's performance"""
    exp_id = df['experiment_id'].iloc[0]
    depth = df['depth'].iloc[0]
    updates = df['num_updates'].iloc[0]
    queries = df['num_queries'].iloc[0]
    
    # Create header panel
    header = Panel(
        f"[bold blue]Experiment:[/bold blue] {exp_id}\n"
        f"[bold green]Configuration:[/bold green] depth={depth}, updates={updates}, queries={queries}",
        title="[bold cyan]Experiment Analysis[/bold cyan]",
        border_style="blue",
        box=box.ROUNDED
    )
    console.print(header)
    
    # Initialization stats
    init_data = df[df['phase'] == 'init']
    if not init_data.empty:
        init_time = init_data[init_data['metric'] == 'initialization_time']['value'].iloc[0]
        console.print(f"\n[bold yellow]Initialization Time:[/bold yellow] [green]{init_time:.2f} ms[/green]")
    
    # Update phase analysis
    update_data = df[df['phase'] == 'update']
    if not update_data.empty:
        console.print("\n[bold magenta]--- Update Operations ---[/bold magenta]")
        
        update_table = Table(title="Update Performance Metrics", box=box.ROUNDED)
        update_table.add_column("Metric", style="cyan", no_wrap=True)
        update_table.add_column("Average", style="green")
        update_table.add_column("Std Dev", style="yellow")
        update_table.add_column("Min", style="blue")
        update_table.add_column("Max", style="red")
        
        # Path computation times
        path_times = update_data[update_data['metric'] == 'path_computation_time']['value']
        if not path_times.empty:
            update_table.add_row(
                "Path Computation (ms)",
                f"{path_times.mean():.2f}",
                f"{path_times.std():.2f}",
                f"{path_times.min():.2f}",
                f"{path_times.max():.2f}"
            )
        
        # Parallel update times
        parallel_times = update_data[update_data['metric'] == 'parallel_updates_time']['value']
        if not parallel_times.empty:
            update_table.add_row(
                "Parallel Updates (ms)",
                f"{parallel_times.mean():.2f}",
                f"{parallel_times.std():.2f}",
                f"{parallel_times.min():.2f}",
                f"{parallel_times.max():.2f}"
            )
        
        console.print(update_table)
        
        # Total update time
        total_update = update_data[update_data['metric'] == 'total_update_time']['value']
        if not total_update.empty:
            console.print(f"[bold yellow]Total Update Time:[/bold yellow] [green]{total_update.sum():.2f}ms[/green] "
                         f"[dim](avg per update: {total_update.mean():.2f}ms)[/dim]")
    
    # Query phase analysis
    query_data = df[df['phase'] == 'query']
    if not query_data.empty:
        console.print("\n[bold magenta]--- Range Sum Query Operations ---[/bold magenta]")
        
        query_table = Table(title="Query Performance Metrics", box=box.ROUNDED)
        query_table.add_column("Metric", style="cyan", no_wrap=True)
        query_table.add_column("Average", style="green")
        query_table.add_column("Std Dev", style="yellow")
        query_table.add_column("Min", style="blue")
        query_table.add_column("Max", style="red")
        
        # Path computation times
        path_times = query_data[query_data['metric'] == 'path_computation_time']['value']
        if not path_times.empty:
            query_table.add_row(
                "Path Computation (ms)",
                f"{path_times.mean():.2f}",
                f"{path_times.std():.2f}",
                f"{path_times.min():.2f}",
                f"{path_times.max():.2f}"
            )
        
        # Direct sum computation times
        sum_times = query_data[query_data['metric'] == 'direct_sum_time']['value']
        if not sum_times.empty:
            query_table.add_row(
                "Direct Sum Computation (ms)",
                f"{sum_times.mean():.2f}",
                f"{sum_times.std():.2f}",
                f"{sum_times.min():.2f}",
                f"{sum_times.max():.2f}"
            )
        
        console.print(query_table)
        
        # Total query time
        total_query = query_data[query_data['metric'] == 'total_query_time']['value']
        if not total_query.empty:
            console.print(f"[bold yellow]Total Query Time:[/bold yellow] [green]{total_query.sum():.2f}ms[/green] "
                         f"[dim](avg per query: {total_query.mean():.2f}ms)[/dim]")
    
    # Network/computation stats
    for phase in ['update_summary', 'query_summary']:
        phase_data = df[df['phase'] == phase]
        if not phase_data.empty:
            console.print(f"\n[bold magenta]--- {phase.replace('_', ' ').title()} ---[/bold magenta]")
            
            stats_table = Table(box=box.SIMPLE)
            stats_table.add_column("Metric", style="cyan")
            stats_table.add_column("Value", style="green")
            
            for metric in ['messages_sent', 'message_bytes', 'aes_operations', 'memory_usage']:
                metric_data = phase_data[phase_data['metric'] == metric]
                if not metric_data.empty:
                    value = metric_data['value'].iloc[0]
                    unit = metric_data['unit'].iloc[0]
                    stats_table.add_row(
                        metric.replace('_', ' ').title(),
                        f"{value:,.0f} {unit}"
                    )
            
            if stats_table.rows:
                console.print(stats_table)

def compare_experiments(csv_files):
    """Compare multiple experiments"""
    header = Panel(
        "[bold cyan]Comparing Multiple Experiments[/bold cyan]",
        border_style="cyan",
        box=box.ROUNDED
    )
    console.print(header)
    
    results = []
    for csv_file in track(csv_files, description="Processing experiments..."):
        df = load_performance_data(csv_file)
        exp_id = df['experiment_id'].iloc[0]
        depth = df['depth'].iloc[0]
        updates = df['num_updates'].iloc[0]
        queries = df['num_queries'].iloc[0]
        
        # Get average times
        update_data = df[(df['phase'] == 'update') & (df['metric'] == 'total_update_time')]
        query_data = df[(df['phase'] == 'query') & (df['metric'] == 'total_query_time')]
        
        avg_update = update_data['value'].mean() if not update_data.empty else 0
        avg_query = query_data['value'].mean() if not query_data.empty else 0
        
        results.append({
            'experiment_id': exp_id,
            'depth': depth,
            'updates': updates,
            'queries': queries,
            'avg_update_time': avg_update,
            'avg_query_time': avg_query,
            'file': csv_file
        })
    
    comparison_df = pd.DataFrame(results)
    
    # Create rich comparison table
    comparison_table = Table(title="Experiment Comparison", box=box.ROUNDED)
    comparison_table.add_column("Experiment ID", style="cyan", no_wrap=True)
    comparison_table.add_column("Depth", style="magenta")
    comparison_table.add_column("Updates", style="blue")
    comparison_table.add_column("Queries", style="blue")
    comparison_table.add_column("Avg Update Time (ms)", style="green")
    comparison_table.add_column("Avg Query Time (ms)", style="yellow")
    
    for _, row in comparison_df.iterrows():
        comparison_table.add_row(
            row['experiment_id'],
            str(row['depth']),
            str(row['updates']),
            str(row['queries']),
            f"{row['avg_update_time']:.2f}",
            f"{row['avg_query_time']:.2f}"
        )
    
    console.print(comparison_table)
    
    # Generate comparison plots
    plot_comparison(comparison_df)
    
    return comparison_df

def plot_comparison(comparison_df, output_dir='logs/plots'):
    """Generate comparison plots for multiple experiments"""
    os.makedirs(output_dir, exist_ok=True)
    
    if len(comparison_df) < 2:
        console.print("[yellow]Warning: Need at least 2 experiments for comparison plots[/yellow]")
        return
    
    # Extract version and sort by depth
    def extract_version_and_depth(exp_id):
        """Extract version (st8, st7) and depth from experiment ID"""
        parts = exp_id.split('_')
        if len(parts) >= 2 and parts[0].startswith('st') and parts[1].startswith('d'):
            version = parts[0]  # e.g., "st8"
            depth = int(parts[1][1:])  # e.g., 20 from "d20"
            return version, depth
        return "unknown", 0
    
    # Add version and depth columns for sorting and grouping
    comparison_df['version'] = comparison_df['experiment_id'].apply(lambda x: extract_version_and_depth(x)[0])
    comparison_df['depth_num'] = comparison_df['experiment_id'].apply(lambda x: extract_version_and_depth(x)[1])
    
    # Sort by version first, then by depth
    comparison_df = comparison_df.sort_values(['version', 'depth_num'])
    
    # Group by version for line plots
    versions = comparison_df['version'].unique()
    colors = ['steelblue', 'coral', 'green', 'purple', 'orange', 'brown', 'pink', 'gray']
    color_map = {version: colors[i % len(colors)] for i, version in enumerate(versions)}
    
    sns.set_style("whitegrid")
    
    # Plot 1: Average Update Times Comparison (Line Plot)
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot lines for each version
    for version in versions:
        version_data = comparison_df[comparison_df['version'] == version].sort_values('depth_num')
        if not version_data.empty:
            ax.plot(version_data['depth_num'], version_data['avg_update_time'], 
                   marker='o', label=version, linewidth=2, markersize=8,
                   color=color_map[version])
            
            # Add value labels on points
            for _, row in version_data.iterrows():
                ax.annotate(f'{row["avg_update_time"]:.1f}ms', 
                           (row['depth_num'], row['avg_update_time']),
                           textcoords="offset points", xytext=(0,10), ha='center',
                           fontsize=9, fontweight='bold')
    
    ax.set_xlabel('Tree Depth', fontsize=12)
    ax.set_ylabel('Average Update Time (ms)', fontsize=12)
    ax.set_title('Average Update Time Comparison Across Experiments', fontsize=14, fontweight='bold')
    ax.legend(title='Segment Tree Version', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/comparison_update_times.png', dpi=300, bbox_inches='tight')
    console.print(f"[green]Saved: {output_dir}/comparison_update_times.png[/green]")
    plt.close()
    
    # Plot 2: Average Query Times Comparison (Line Plot)
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot lines for each version
    for version in versions:
        version_data = comparison_df[comparison_df['version'] == version].sort_values('depth_num')
        if not version_data.empty:
            ax.plot(version_data['depth_num'], version_data['avg_query_time'], 
                   marker='s', label=version, linewidth=2, markersize=8,
                   color=color_map[version])
            
            # Add value labels on points
            for _, row in version_data.iterrows():
                ax.annotate(f'{row["avg_query_time"]:.1f}ms', 
                           (row['depth_num'], row['avg_query_time']),
                           textcoords="offset points", xytext=(0,10), ha='center',
                           fontsize=9, fontweight='bold')
    
    ax.set_xlabel('Tree Depth', fontsize=12)
    ax.set_ylabel('Average Query Time (ms)', fontsize=12)
    ax.set_title('Average Range Query Time Comparison Across Experiments', fontsize=14, fontweight='bold')
    ax.legend(title='Segment Tree Version', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/comparison_query_times.png', dpi=300, bbox_inches='tight')
    console.print(f"[green]Saved: {output_dir}/comparison_query_times.png[/green]")
    plt.close()
    
    # Plot 3: Side-by-side comparison (Line Plots)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Update times subplot
    for version in versions:
        version_data = comparison_df[comparison_df['version'] == version].sort_values('depth_num')
        if not version_data.empty:
            ax1.plot(version_data['depth_num'], version_data['avg_update_time'], 
                    marker='o', label=version, linewidth=2, markersize=6,
                    color=color_map[version])
    
    ax1.set_xlabel('Tree Depth', fontsize=11)
    ax1.set_ylabel('Average Update Time (ms)', fontsize=11)
    ax1.set_title('Average Update Times', fontsize=12, fontweight='bold')
    ax1.legend(title='Version', fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    # Query times subplot
    for version in versions:
        version_data = comparison_df[comparison_df['version'] == version].sort_values('depth_num')
        if not version_data.empty:
            ax2.plot(version_data['depth_num'], version_data['avg_query_time'], 
                    marker='s', label=version, linewidth=2, markersize=6,
                    color=color_map[version])
    
    ax2.set_xlabel('Tree Depth', fontsize=11)
    ax2.set_ylabel('Average Query Time (ms)', fontsize=11)
    ax2.set_title('Average Range Query Times', fontsize=12, fontweight='bold')
    ax2.legend(title='Version', fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/comparison_side_by_side.png', dpi=300, bbox_inches='tight')
    console.print(f"[green]Saved: {output_dir}/comparison_side_by_side.png[/green]")
    plt.close()

def plot_performance(df, output_dir='logs/plots'):
    """Generate performance visualization plots"""
    os.makedirs(output_dir, exist_ok=True)
    exp_id = df['experiment_id'].iloc[0]
    
    sns.set_style("whitegrid")
    
    # Plot 1: Update operation breakdown
    update_data = df[df['phase'] == 'update']
    if not update_data.empty:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        path_times = update_data[update_data['metric'] == 'path_computation_time']
        parallel_times = update_data[update_data['metric'] == 'parallel_updates_time']
        
        if not path_times.empty and not parallel_times.empty:
            x = path_times['operation_id']
            ax.plot(x, path_times['value'], marker='o', label='Path Computation', linewidth=2)
            ax.plot(x, parallel_times['value'], marker='s', label='Parallel Updates', linewidth=2)
            
            ax.set_xlabel('Update Operation ID', fontsize=12)
            ax.set_ylabel('Time (ms)', fontsize=12)
            ax.set_title(f'Update Operation Breakdown - {exp_id}', fontsize=14)
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(f'{output_dir}/update_breakdown_{exp_id}.png', dpi=300)
            console.print(f"[green]Saved: {output_dir}/update_breakdown_{exp_id}.png[/green]")
            plt.close()
    
    # Plot 2: Query operation breakdown
    query_data = df[df['phase'] == 'query']
    if not query_data.empty:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        path_times = query_data[query_data['metric'] == 'path_computation_time']
        sum_times = query_data[query_data['metric'] == 'direct_sum_time']
        
        if not path_times.empty and not sum_times.empty:
            x = path_times['operation_id']
            ax.plot(x, path_times['value'], marker='o', label='Path Computation', linewidth=2)
            ax.plot(x, sum_times['value'], marker='s', label='Direct Sum Computation', linewidth=2)
            
            ax.set_xlabel('Query Operation ID', fontsize=12)
            ax.set_ylabel('Time (ms)', fontsize=12)
            ax.set_title(f'Query Operation Breakdown - {exp_id}', fontsize=14)
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(f'{output_dir}/query_breakdown_{exp_id}.png', dpi=300)
            console.print(f"[green]Saved: {output_dir}/query_breakdown_{exp_id}.png[/green]")
            plt.close()
    
    # Plot 3: Total time comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    if not update_data.empty:
        total_update = update_data[update_data['metric'] == 'total_update_time']
        if not total_update.empty:
            ax1.bar(total_update['operation_id'], total_update['value'], color='steelblue', alpha=0.7)
            ax1.set_xlabel('Update Operation ID', fontsize=11)
            ax1.set_ylabel('Time (ms)', fontsize=11)
            ax1.set_title('Total Update Time per Operation', fontsize=12)
            ax1.grid(True, alpha=0.3, axis='y')
    
    if not query_data.empty:
        total_query = query_data[query_data['metric'] == 'total_query_time']
        if not total_query.empty:
            ax2.bar(total_query['operation_id'], total_query['value'], color='coral', alpha=0.7)
            ax2.set_xlabel('Query Operation ID', fontsize=11)
            ax2.set_ylabel('Time (ms)', fontsize=11)
            ax2.set_title('Total Query Time per Operation', fontsize=12)
            ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/total_times_{exp_id}.png', dpi=300)
    console.print(f"[green]Saved: {output_dir}/total_times_{exp_id}.png[/green]")
    plt.close()

def main():
    if len(sys.argv) < 2:
        console.print("[red]Usage:[/red] python3 analyze_performance.py <csv_file(s)>")
        console.print("[yellow]Example:[/yellow] python3 analyze_performance.py logs/performance_exp001.csv")
        console.print("[yellow]Example:[/yellow] python3 analyze_performance.py logs/performance_*.csv")
        sys.exit(1)
    
    # Expand glob patterns
    csv_files = []
    for pattern in sys.argv[1:]:
        csv_files.extend(glob.glob(pattern))
    
    if not csv_files:
        console.print("[red]No CSV files found![/red]")
        sys.exit(1)
    
    # Welcome message
    welcome = Panel(
        f"[bold green]Performance Analysis Tool[/bold green]\n"
        f"Found {len(csv_files)} experiment file(s) to analyze",
        title="[bold cyan]Segment Tree Performance Analyzer[/bold cyan]",
        border_style="green",
        box=box.ROUNDED
    )
    console.print(welcome)
    
    # Analyze each experiment
    for csv_file in csv_files:
        console.print(f"\n[bold blue]Processing:[/bold blue] [cyan]{csv_file}[/cyan]")
        df = load_performance_data(csv_file)
        analyze_single_experiment(df)
        
        # Generate plots
        try:
            plot_performance(df)
        except Exception as e:
            console.print(f"[yellow]Warning: Could not generate plots: {e}[/yellow]")
    
    # Compare if multiple experiments
    if len(csv_files) > 1:
        console.print("\n" + "="*60)
        compare_experiments(csv_files)

if __name__ == "__main__":
    main()
