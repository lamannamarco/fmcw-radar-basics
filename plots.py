"""Figures of the README (see main.py for the simulation and the processing).

plot_scene, plot_if_signal, plot_range_doppler, plot_virtual_array and plot_range_angle draw
the figures of Section 2 of the README.
"""
import numpy as np
import matplotlib.pyplot as plt


def to_db(x, dyn=60):
    x = 10 * np.log10(x + 1e-12)
    return np.clip(x, x.max() - dyn, None) - x.max()


def draw_array(ax, p):
    """TX, RX and the resulting virtual array, on the axes `ax`."""
    N_tx, N_rx, lam, d, x_tx, x_rx = p.N_tx, p.N_rx, p.lam, p.d, p.x_tx, p.x_rx
    tx_colors = plt.cm.tab10(np.arange(N_tx))
    ax.scatter(x_tx / lam, np.full(N_tx, 2), marker="s", s=90, c=tx_colors)
    ax.scatter(x_rx / lam, np.full(N_rx, 1), marker="o", s=90, c="gray")
    virt_pos = (x_tx[:, None] + x_rx[None, :]) / lam
    ax.scatter(virt_pos.ravel(), np.zeros(N_tx * N_rx), marker="D", s=90,
               c=np.repeat(tx_colors, N_rx, axis=0))
    for k in range(N_tx):
        ax.annotate(f"TX{k}", (x_tx[k] / lam, 2.18), ha="center", color=tx_colors[k])
    for r in range(N_rx):
        ax.annotate(f"RX{r}", (x_rx[r] / lam, 1.18), ha="center", color="gray")
    ax.annotate("", (x_tx[1] / lam, 2.5), (x_tx[0] / lam, 2.5),
                arrowprops=dict(arrowstyle="<->"))
    ax.text(x_tx[1] / lam / 2, 2.55, f"{N_rx}·λ/2 = {N_rx * d / lam:.1f} λ", ha="center")
    ax.annotate("", (x_rx[1] / lam, 1.5), (x_rx[0] / lam, 1.5), arrowprops=dict(arrowstyle="<->"))
    ax.text(x_rx[1] / lam / 2, 1.55, "λ/2", ha="center")
    ax.text(virt_pos.max() / lam * lam + 0.25, 0, "virtual array:\nposition = TX + RX,\n"
            f"{N_tx * N_rx} elements, λ/2 apart\n(color = which TX)", va="center", fontsize=8)
    ax.set_yticks([0, 1, 2], ["virtual", "RX", "TX"])
    ax.set_xlabel("position along array (λ)")
    ax.set_ylim(-0.9, 3)
    ax.set_xlim(-0.5, virt_pos.max() + 2.2)
    ax.grid(alpha=0.3)


def draw_scene(ax, p):
    """The radar at the origin and the targets, seen from above, on the axes `ax`."""
    targets, R_max = p.targets, p.R_max
    tgt_colors = ["tab:red", "tab:blue"]

    a_ = np.radians(np.linspace(-90, 90, 181))
    for rr_ in (100, 200, 300):
        ax.plot(rr_ * np.sin(a_), rr_ * np.cos(a_), color="gray", lw=0.6)
        ax.text(0, rr_ + 4, f"{rr_} m", ha="center", fontsize=7, color="gray")
    ax.plot(0, 0, "k^", ms=12)
    ax.annotate("Radar", (0, 0), xytext=(0, -16), textcoords="offset points", ha="center")
    for i, ((R0, ang, *_), col) in enumerate(zip(targets, tgt_colors), start=1):
        x, y = R0 * np.sin(np.radians(ang)), R0 * np.cos(np.radians(ang))
        ax.plot(x, y, "o", color=col, ms=10)
        ax.annotate(f"Target {i}", (x, y), xytext=(0, 10), textcoords="offset points",
                    ha="center", color=col, fontweight="bold")
    ax.set_aspect("equal")
    ax.set_xlim(-R_max, R_max)
    ax.set_ylim(-0.08 * R_max, R_max * 1.05)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("The scene (top view)")


def plot_scene(p):
    """Step 1: the radar at the origin and the targets, seen from above."""
    fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
    draw_scene(ax, p)
    return fig


def plot_if_signal(p):
    """Step 2: I/Q samples of one chirp (left) and the magnitude of their FFT (right)."""
    targets, R_max, c, alpha, fs, N_adc = p.targets, p.R_max, p.c, p.alpha, p.fs, p.N_adc
    tgt_colors = ["tab:red", "tab:blue"]

    x0 = p.cube[0, 0]  # first virtual element, first chirp
    t_us = (p.T_adc_start + p.t_fast) * 1e6
    spec = np.abs(np.fft.fft(x0 * p.win_r))[:N_adc // 2]
    range_ = np.arange(N_adc // 2) * fs / N_adc * c / (2 * alpha)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.8), layout="constrained")
    ax1.plot(t_us, x0.real, lw=0.8, label="I (real)")
    ax1.plot(t_us, x0.imag, lw=0.8, label="Q (imag)")
    ax1.set_xlabel("time within chirp (µs)")
    ax1.set_ylabel("amplitude")
    ax1.set_title("Time domain")
    ax1.set_xlim(0, p.Tc * 1e6)
    ax1.legend(loc="upper right")
    ax1.grid(alpha=0.3)

    spec_db = 20 * np.log10(spec / spec.max())
    ax2.plot(range_, spec_db, lw=0.9, color="k")
    for i, ((R0, *_), col) in enumerate(zip(targets, tgt_colors), start=1):
        ax2.axvline(R0, color=col, ls="--", lw=1)
        ax2.annotate(f"Target {i}\n{R0:.0f} m", (R0, 0), xytext=(6, -4),
                     textcoords="offset points", va="top", color=col, fontweight="bold")
    ax2.set_xlim(0, R_max)
    ax2.set_ylim(-70, 3)
    ax2.set_xlabel("range (m)")
    ax2.set_ylabel("|FFT| (dB, normalized)")
    ax2.set_title("Frequency / range domain")
    ax2.grid(alpha=0.3)
    return fig


def plot_range_doppler(p):
    """Section 2.2: raw (I) data, FFT along fast time and FFT along slow time (moving targets)."""
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), layout="constrained")
    Q_demo, v_demo = p.Q_demo, p.v_demo
    t_us = (p.T_adc_start + p.t_fast) * 1e6

    axs[0].imshow(p.x_demo.real, aspect="auto", cmap="RdBu_r", origin="lower",
                  extent=[t_us[0], t_us[-1], -0.5, Q_demo - 0.5])
    axs[0].set_xlim(0, p.Tc * 1e6)
    axs[0].set_xlabel("fast time (µs)")
    axs[0].set_ylabel("slow time (chirp index)")
    axs[0].set_title("Raw data (I)")

    mag = 20 * np.log10(np.abs(p.rfft_demo) / np.abs(p.rfft_demo).max() + 1e-9)
    axs[1].imshow(np.clip(mag, -50, 0), aspect="auto", cmap="inferno", origin="lower",
                  extent=[p.range_axis[0], p.range_axis[-1], -0.5, Q_demo - 0.5])
    axs[1].set_xlabel("range (m)")
    axs[1].set_ylabel("slow time (chirp index)")
    axs[1].set_title("FFT along fast time")

    im = axs[2].pcolormesh(p.range_axis, p.vel_axis, to_db(np.abs(p.rdm) ** 2, dyn=50),
                           shading="auto", cmap="inferno")
    for (R0, *_), v in zip(p.targets, v_demo):
        axs[2].plot(R0, v, "o", mfc="none", mec="w", ms=14)
    axs[2].set_xlabel("range (m)")
    axs[2].set_ylabel("velocity (m/s)")
    axs[2].set_title("FFT along slow time: range-Doppler map")
    fig.colorbar(im, ax=axs[2], label="dB (normalized)")
    return fig


def plot_virtual_array(p):
    """Step 6: TX, RX and the resulting virtual array."""
    fig, ax = plt.subplots(figsize=(10, 3.6), layout="constrained")
    draw_array(ax, p)
    ax.set_title("Physical TX and RX antennas, and the resulting virtual array")
    return fig


def plot_range_angle(p):
    """Step 9: the scene (left) and the range-angle map reconstructed from the data (right)."""
    R_max = p.R_max
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6), layout="constrained")
    draw_scene(ax1, p)

    th, rr = np.meshgrid(np.radians(p.ang_axis), p.range_axis, indexing="ij")
    im = ax2.pcolormesh(rr * np.sin(th), rr * np.cos(th), to_db(p.ra_map, dyn=50),
                        shading="gouraud", cmap="inferno")
    ax2.set_aspect("equal")
    ax2.set_facecolor("k")
    ax2.set_xlim(-R_max, R_max)
    ax2.set_ylim(-0.08 * R_max, R_max * 1.05)
    ax2.set_xlabel("x (m)")
    ax2.set_ylabel("y (m)")
    ax2.set_title("Reconstructed image: range-angle map")
    fig.colorbar(im, ax=ax2, label="dB (normalized)", shrink=0.8)
    return fig
