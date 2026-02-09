from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import matplotlib as mpl
import mpltern  # noqa: F401
import numpy as np
from pathlib import Path
from mpi4py import MPI
from spparks import spparks

def plot_ternary(ax: Axes, ts, coverage, *, ts_min=-1, ts_max=-1):
    if ts_min == -1:
        vmin = ts.min()
    else:
        vmin = ts_min
    if ts_max == -1:
        vmax = ts.max()
    else:
        vmax = ts_max

    norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    cmap = plt.cm.viridis
    sc_handles = []

    sc = ax.scatter(
        coverage[:, 0],
        coverage[:, 1],
        coverage[:, 2],
        c=ts,
        cmap=cmap,
        norm=norm,
        s=4,          # point size
        edgecolors="none",
    )
    sc_handles.append(sc)

    #ax.quiver(space_converage[:,0],space_converage[:,1],space_converage[:,2],coverage_tang[:,0],coverage_tang[:,1],coverage_tang[:,2], width=0.005,zorder=10)
    ax.set_tlabel("$\\theta_*$")
    ax.set_llabel("$\\theta_A$")
    ax.set_rlabel("$\\theta_B$")

    # Gridlines
    ax.grid()

    # Color ticks, grids, tick-labels
    ax.taxis.set_tick_params(tick2On=True, colors='C0', grid_color='C0')
    ax.laxis.set_tick_params(tick2On=True, colors='C1', grid_color='C1')
    ax.raxis.set_tick_params(tick2On=True, colors='C2', grid_color='C2')

    # Color labels
    ax.taxis.label.set_color('C0')
    ax.laxis.label.set_color('C1')
    ax.raxis.label.set_color('C2')

    # Color spines
    ax.spines['tside'].set_color('C0')
    ax.spines['lside'].set_color('C1')
    ax.spines['rside'].set_color('C2')

    return sc_handles

def _coverage_from_site(site, nglobal, comm):
    local_counts = np.bincount(site, minlength=3).astype(np.int64)
    global_counts = np.zeros_like(local_counts)
    comm.Allreduce(local_counts, global_counts, op=MPI.SUM)
    return global_counts / float(nglobal)


def main(total_time=5, dt_min=1e-5):
    if hasattr(MPI,"Is_initialized") and not MPI.Is_initialized():
        MPI.Init()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    spk = spparks(comm=comm)

    input_file = Path(__file__).with_name("in.oscillations")
    spk.file(str(input_file))

    site = spk.extract_numpy("site", 1)
    nglobal = spk.extract("nglobal", 0)

    times = []
    coverages = []

    t = spk.extract("time", 3)
    while t < total_time - 1e-12:
        dt = min(dt_min, total_time - t)
        spk.command(f"run {dt}")
        t = spk.extract("time", 3)

        coverage = _coverage_from_site(site, nglobal, comm)
        times.append(t)
        coverages.append(coverage)

    times = np.asarray(times)
    coverages = np.vstack(coverages) if coverages else np.zeros((0, 3))

    if rank == 0:
        fig_ternary = plt.figure(constrained_layout=True, dpi=300)
        ax = fig_ternary.add_subplot(111, projection="ternary")
        sc_handles = plot_ternary(ax, times, coverages)
        cbar = fig_ternary.colorbar(sc_handles[0], ax=ax, shrink=0.50, pad=0.1)
        cbar.set_label("Time")
        plt.show()


if __name__ == "__main__":
    main()
