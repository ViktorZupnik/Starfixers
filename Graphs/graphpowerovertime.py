import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
df = pd.read_excel('power_graph.xlsx', index_col=0)

# Extract and clean the time row
time_str = df.loc['Time']
time = pd.to_numeric(time_str.astype(str).str.replace(',', ''), errors='coerce')

# Drop the 'Time' row
df = df.drop('Time')

repeated_phases = ['Outer Orbit', 'Observe', 'Prepare', 'De-Orbit']
all_phases = df.columns.tolist()
singular_phases = [phase for phase in all_phases if phase not in repeated_phases]

def plot_singular_with_gap(singular_phases, gap_between=('Startup', 'EOL'), gap_size=300):
    time_sing = time[singular_phases]
    gap_idx = singular_phases.index(gap_between[0])
    
    # Calculate bar positions with gap
    bar_positions = []
    current_pos = 0
    for i, phase in enumerate(singular_phases):
        bar_positions.append(current_pos)
        current_pos += time_sing[phase]
        if i == gap_idx:
            current_pos += gap_size
    
    bar_positions = np.array(bar_positions)
    df_sing = df[singular_phases].T

    plt.figure(figsize=(12, 7))
    bottom = np.zeros(len(df_sing))

    for subsystem in df_sing.columns:
        plt.bar(bar_positions, df_sing[subsystem], bottom=bottom,
                width=time_sing, align='edge', label=subsystem)
        bottom += df_sing[subsystem]

    plt.xlabel('Time (s)')
    plt.ylabel('Power (W)')
    plt.title('Non-repeated Mission Phases')
    plt.grid(axis='y')

    # Add correct phase labels only above actual bars
    for pos, phase, duration in zip(bar_positions, singular_phases, time_sing):
        plt.text(pos + duration / 2, bottom.max() * 0.9, phase,
                 ha='center', va='bottom', fontsize=10, rotation=45)

    # Add vertical dotted lines for the gap
    gap_start = bar_positions[gap_idx] + time_sing[gap_between[0]]
    gap_end = gap_start + gap_size
    plt.axvline(x=gap_start, color='black', linestyle=':', linewidth=1.5)
    plt.axvline(x=gap_end, color='black', linestyle=':', linewidth=1.5)

    # Label between dotted lines
    plt.text((gap_start + gap_end) / 2 + 900, bottom.max() * 0.9, 'Repeated Phases',
             ha='center', va='bottom', fontsize=12)

    plt.legend(title='Subsystem', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()


def plot_repeated(repeated_phases):
    time_rep = time[repeated_phases]
    cumulative_time = time_rep.cumsum()
    bar_positions = cumulative_time - time_rep
    
    df_rep = df[repeated_phases].T
    
    plt.figure(figsize=(12,7))
    bottom = np.zeros(len(df_rep))
    
    for subsystem in df_rep.columns:
        plt.bar(bar_positions, df_rep[subsystem], bottom=bottom,
                width=time_rep, align='edge', label=subsystem)
        bottom += df_rep[subsystem]
    
    plt.xlabel('Cumulative Time')
    plt.ylabel('Power')
    plt.title('Repeated Mission Phases')
    plt.grid(axis='y')
    
    # Add phase labels
    for pos, phase, duration in zip(bar_positions, df_rep.index, time_rep):
        plt.text(pos + duration/2, bottom.max() * 0.9, phase,
                 ha='center', va='bottom', fontsize=10, rotation=45)
    
    plt.legend(title='Subsystem', bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()
    plt.show()

# Plot singular phases with gap and repeated phases normally
plot_singular_with_gap(singular_phases)
plot_repeated(repeated_phases)
