"""FMCW radar basics: 77 GHz TDM-MIMO radar (4 TX x 3 RX) observing two point targets.

The code follows Section 2 of the README, and generates its figures (docs/):
  2.1  IF signal of one chirp and its FFT                    -> 01_scene, 03_if_signal
  2.2  Train of chirps: fast/slow time, range-Doppler map    -> 05_range_doppler
  2.3  MIMO virtual array, data cube, range-angle map        -> 07_virtual_array, 09_range_angle

Conventions: complex (I/Q) sampling, positive angle = to the right of boresight,
positive velocity = receding target. The notation of the README is used (Q chirps, N_adc samples,
l-th TX, m-th RX, p-th virtual element).
"""
from types import SimpleNamespace
import numpy as np
import matplotlib.pyplot as plt
from plots import (plot_if_signal, plot_range_angle, plot_range_doppler, plot_scene,
                   plot_virtual_array)
from scipy.constants import speed_of_light as c

rng = np.random.default_rng(0)

# %% 2.1 Radar parameters (Section 2.1.2)
fc = 77e9                # [Hz] Carrier frequency
lam = c / fc             # [m] Carrier wavelength
B = 120e6                # [Hz] Sweep bandwidth
Tc = 12e-6               # [s] Chirp duration
alpha = B / Tc           # [Hz/s] Chirp slope (10 MHz/us)
fs = 45e6                # [Hz] ADC sampling rate (I/Q, so the IF frequencies go up to fs)
Ts = 1 / fs              # [s] Sampling interval
N_adc = 512              # Samples per chirp
T_adc = N_adc * Ts       # [s] ADC window (11.4 us), sits inside the chirp
T_adc_start = 0.5e-6     # [s] Settling time before the ADC starts
R_max = 300.0            # [m] Range shown in the figures (f_IF,max = fs would be c*fs/(2*alpha) = 675 m)

Q = 128                  # Chirps per TX; they are sent back to back (one every Tc)
t_fast = np.arange(N_adc) * Ts  # [s] fast time (sample n)

# %% 2.3 Array: 4 TX and 3 RX on a line (x axis), see Section 2.3.3
N_tx, N_rx = 4, 3
N_ve = N_tx * N_rx                # Virtual elements
d = lam / 2                       # [m] Virtual element spacing
x_tx = np.arange(N_tx) * N_rx * d  # TX spaced N_rx * lam/2 apart
x_rx = np.arange(N_rx) * d         # RX spaced lam/2 apart
x_ve = (x_tx[:, None] + x_rx[None, :]).ravel()  # virtual element p = l*N_rx + m sits at x_tx[l] + x_rx[m]
T_p = N_tx * Tc                   # TDM period

B_adc = alpha * T_adc             # Bandwidth seen by the ADC
print(f"range resolution   : {c / (2 * B_adc):.2f} m (ADC window; c/2B = {c / (2 * B):.2f} m)")
print(f"IF at {R_max:.0f} m        : {2 * alpha * R_max / c / 1e6:.1f} MHz (fs = {fs / 1e6:.1f} MHz)")
print(f"angular resolution : ~{np.degrees(2 / N_ve):.1f} deg at broadside")

# Targets: (range m, angle deg, per-sample SNR dB). Static, unless stated otherwise.
# Add more lines if you want more targets!
targets = [
    (75.0, -30.0, 20.0),
    (150.0, 45.0, 10.0),
]


# %% 2.1-2.3 Simulation
def simulate(targets, velocities, n_chirps, x_ve):
    """Raw IF data as a cube, shape (len(x_ve), n_chirps, N_adc): (element, chirp, sample).

    The target at range R and velocity v is at R + v*q*Tc during chirp q. Virtual element p, at
    position x_ve[p], sees the extra phase 2*pi*x_ve[p]*sin(theta)/lam. Complex white Gaussian
    noise is added.
    """
    t = T_adc_start + t_fast  # time within the chirp
    q = np.arange(n_chirps)[:, None]  # chirp index (slow time)
    data = np.zeros((len(x_ve), n_chirps, N_adc), complex)
    for (R0, ang, snr), v in zip(targets, velocities):
        tau = 2 * (R0 + v * q * Tc) / c  # round-trip delay
        ph_time = 2 * np.pi * (fc * tau + alpha * tau * t)  # (chirp, sample)
        ph_space = 2 * np.pi * x_ve * np.sin(np.radians(ang)) / lam  # (element,)
        data += 10 ** (snr / 20) * np.exp(1j * (ph_time + ph_space[:, None, None]))
    noise = rng.standard_normal(data.shape) + 1j * rng.standard_normal(data.shape)
    return data + noise / np.sqrt(2)


win_r = np.hanning(N_adc)
dR_bin = c * fs / (2 * alpha * N_adc)  # [m] bin k <-> f_IF = k*fs/N_fft, N_fft = N_adc
n_bins = int(np.ceil(R_max / dR_bin))
range_axis = np.arange(n_bins) * dR_bin

# %% 2.2 Range-Doppler map (illustration only: here the targets move)
Q_demo = 32
v_demo = (20.0, -10.0)  # [m/s]
x_demo = simulate(targets, v_demo, Q_demo, np.zeros(1))[0]  # (chirp, sample)
rfft_demo = np.fft.fft(x_demo * win_r, axis=1)[:, :n_bins]  # FFT along fast time
rdm = np.fft.fftshift(np.fft.fft(rfft_demo * np.hanning(Q_demo)[:, None], axis=0), axes=0)  # slow time
vel_axis = (np.arange(Q_demo) - Q_demo // 2) / (Q_demo * Tc) * lam / 2

# %% 2.3 Data cube and range-angle map (static targets)
cube = simulate(targets, np.zeros(len(targets)), Q, x_ve)  # (N_ve, Q, N_adc) = (12, Q, N_adc)

# Range FFT along fast time, then average over the chirps (the zero-Doppler bin: static targets)
virt = np.fft.fft(cube * win_r, axis=2)[..., :n_bins].mean(axis=1)  # (N_ve, range)

# Angle FFT across the virtual elements: bin k_a <-> sin(theta) = 2*k_a/N_ang for d = lam/2
N_ang = 128
win_a = np.hanning(N_ve)
afft = np.fft.fftshift(np.fft.fft(virt * win_a[:, None], n=N_ang, axis=0), axes=0)
ra_map = np.abs(afft) ** 2  # (angle, range)
sin_axis = np.append(2 * (np.arange(N_ang) - N_ang // 2) / N_ang, 1.0)  # +1.0 closes the polar plot
ra_map = np.vstack([ra_map, ra_map[:1]])
ang_axis = np.degrees(np.arcsin(np.clip(sin_axis, -1, 1)))

# %% Plots
SAVE = True
SHOW = False

p = SimpleNamespace(**globals())
for name, fn in (("01_scene", plot_scene), ("03_if_signal", plot_if_signal),
                 ("05_range_doppler", plot_range_doppler),
                 ("07_virtual_array", plot_virtual_array), ("09_range_angle", plot_range_angle)):
    fig = fn(p)
    if SAVE:
        fig.savefig(f"docs/{name}.png", dpi=200)
    if not SHOW:
        plt.close(fig)
if SHOW:
    plt.show()
