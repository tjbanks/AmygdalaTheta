import numpy as np
import pandas as pd
import scipy.signal as ss
import matplotlib.pyplot as plt
import h5py
import sys
import os
from Analysis.build_input import get_populations, get_stim_cycle, t_stop

MODEL_PATH = os.path.join('..')
CONFIG = 'config.json'
pop_color = {'PyrA': 'blue', 'PyrC': 'green', 'PV': 'red', 'SOM': 'purple', 'CR': 'orange'}
pop_color = {p: 'tab:' + clr for p, clr in pop_color.items()}
pop_names = list(pop_color.keys())


def raster(pop_spike, pop_color, id_column='node_ids', s=0.1, ax=None):
    if ax is None:
        _, ax = plt.subplots(1, 1)
    ymin, ymax = [], []
    for p, spike_df in pop_spike.items():
        if len(spike_df) > 0:
            ids = spike_df[id_column].values
            ax.scatter(spike_df['timestamps'], ids, c=pop_color[p], s=s, label=p)
            ymin.append(ids.min())
            ymax.append(ids.max())
    ax.set_xlim(left=0.)
    ax.set_ylim([np.min(ymin) - 1, np.max(ymax) + 1])
    ax.set_title('Spike Raster Plot')
    ax.legend(loc='upper right', framealpha=0.9, markerfirst=False)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Cell ID')
    return ax


def firing_rate_histogram(pop_fr, pop_color, bins=30, min_fr=None,
                          logscale=False, stacked=True, ax=None):
    if min_fr is not None:
        pop_fr = {p: np.fmax(fr, min_fr) for p, fr in pop_fr.items()}
    fr = np.concatenate(list(pop_fr.values()))
    if logscale:
        fr = fr[np.nonzero(fr)[0]]
        bins = np.geomspace(fr.min(), fr.max(), bins + 1)
    else:
        bins = np.linspace(fr.min(), fr.max(), bins + 1)
    pop_names = list(pop_fr.keys())
    if ax is None:
        _, ax = plt.subplots(1, 1)
    if stacked:
        ax.hist(pop_fr.values(), bins=bins, label=pop_names,
                color=[pop_color[p] for p in pop_names], stacked=True)
    else:
        for p, fr in pop_fr.items():
            ax.hist(fr, bins=bins, label=p, color=pop_color[p], alpha=0.5)
    if logscale:
        ax.set_xscale('log')
        plt.draw()
        xt = ax.get_xticks()
        xtl = [x.get_text() for x in ax.get_xticklabels()]
        xt = np.append(xt, min_fr)
        xtl.append('0')
        ax.set_xticks(xt)
        ax.set_xticklabels(xtl)
    ax.set_xlim(bins[0], bins[-1])
    ax.legend(loc='upper right')
    ax.set_title('Firing Rate Histogram')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Count')
    return ax


def firing_rate(spikes_df, num_cells=None, time_windows=(0.,), frequency=True):
    """
    Count number of spikes for each cell.
    spikes_df: dataframe of node id and spike times (ms)
    num_cells: number of cells (that determines maximum node id)
    time_windows: list of time windows for counting spikes (second)
    frequency: whether return firing frequency in Hz or just number of spikes
    """
    if not spikes_df['timestamps'].is_monotonic:
        spikes_df = spikes_df.sort_values(by='timestamps')
    if num_cells is None:
        num_cells = spikes_df['node_ids'].max() + 1
    time_windows = 1000. * np.asarray(time_windows).ravel()
    if time_windows.size % 2:
        time_windows = np.append(time_windows, spikes_df['timestamps'].max())
    nspk = np.zeros(num_cells, dtype=int)
    n, N = 0, time_windows.size
    count = False
    for t, i in zip(spikes_df['timestamps'], spikes_df['node_ids']):
        while n < N and t > time_windows[n]:
            n += 1
            count = not count
        if count:
            nspk[i] = nspk[i] + 1
    if frequency:
        nspk = nspk / (total_duration(time_windows) / 1000)
    return nspk


def total_duration(time_windows):
    return np.diff(np.reshape(time_windows, (-1, 2)), axis=1).sum()


def pop_spike_rate(spike_times, time, frequeny=False):
    t = np.arange(*time)
    t = np.append(t, t[-1] + time[2])
    spike_rate, _ = np.histogram(np.asarray(spike_times), t)
    if frequeny:
        spike_rate = 1000 / time[2] * spike_rate
    return spike_rate


def xcorr_coeff(x, y, max_lag=None, dt=1., plot=True, ax=None):
    x = np.asarray(x)
    y = np.asarray(y)
    xcorr = ss.correlate(x, y) / x.size / x.std() / y.std()
    xcorr_lags = ss.correlation_lags(x.size, y.size)
    if max_lag is not None:
        lag_idx = np.nonzero(np.abs(xcorr_lags) <= max_lag / dt)[0]
        xcorr = xcorr[lag_idx]
        xcorr_lags = xcorr_lags[lag_idx]

    if plot:
        if ax is None:
            _, ax = plt.subplots(1, 1)
        ax.plot(dt * xcorr_lags, xcorr)
        ax.set_xlabel('Lag')
        ax.set_ylabel('Cross Correlation')
    return xcorr, xcorr_lags


def get_seg_on_stimulus(x, fs, on_time, off_time,
                        t_start, t=t_stop, tseg=None):
    x = np.asarray(x)
    in_dim = x.ndim
    if in_dim == 1:
        x = x.reshape(1, x.size)
    t = np.asarray(t)
    t_stop = t.size / fs if t.ndim else t
    if tseg is None:
        tseg = on_time # time segment length for PSD (second)
    t_cycle, n_cycle = get_stim_cycle(on_time, off_time, t_start, t_stop)

    nfft = int(tseg * fs) # steps per segment
    i_start = int(t_start * fs)
    i_on = int(on_time * fs)
    i_cycle = int(t_cycle * fs)
    nseg_cycle = int(np.ceil(i_on / nfft))
    x_on = np.zeros((x.shape[0], n_cycle * nseg_cycle * nfft))

    for i in range(n_cycle):
        m = i_start + i * i_cycle
        for j in range(nseg_cycle):
            xx = x[:, m + j * nfft:m + min((j + 1) * nfft, i_on)]
            n = (i * nseg_cycle + j) * nfft
            x_on[:, n:n + xx.shape[1]] = xx
    if in_dim == 1:
        x_on = x_on.ravel()

    stim_cycle = {'t_cycle': t_cycle, 'n_cycle': n_cycle,
                  't_start': t_start, 'on_time': on_time,
                  'i_start': i_start, 'i_cycle': i_cycle}
    return x_on, nfft, stim_cycle


def get_psd_on_stimulus(x, fs, on_time, off_time,
                        t_start, t=t_stop, tseg=None):
    x_on, nfft, stim_cycle = get_seg_on_stimulus(
        x, fs, on_time, off_time, t_start, t=t, tseg=tseg)
    f, pxx = ss.welch(x_on, fs=fs, window='boxcar', nperseg=nfft, noverlap=0)
    return f, pxx, stim_cycle


def get_coh_on_stimulus(x, y, fs, on_time, off_time,
                        t_start, t=t_stop, tseg=None):
    xy = np.array([x, y])
    xy_on, nfft, _ = get_seg_on_stimulus(
        xy, fs, on_time, off_time, t_start, t=t, tseg=tseg)
    f, cxy = ss.coherence(xy_on[0], xy_on[1], fs=fs,
        window='boxcar', nperseg=nfft, noverlap=0)
    return f, cxy


def plot_stimulus_cycles(t, x, stim_cycle, dv_n_sigma=5.,
                         var_label='LFP (mV)', ax=None):
    t_cycle, n_cycle = stim_cycle['t_cycle'], stim_cycle['n_cycle']
    t_start, on_time = stim_cycle['t_start'], stim_cycle['on_time']
    i_start, i_cycle = stim_cycle['i_start'], stim_cycle['i_cycle']
    dv = dv_n_sigma * np.std(x[i_start:])
    
    if ax is None:
        _, ax = plt.subplots(1, 1)
    ax.plot(t[:i_start], x[:i_start] + n_cycle * dv, 'k', label='pre stimulus')
    for i in range(n_cycle):
        m = i_start + i * i_cycle
        xx = x[m:m + i_cycle] + (n_cycle - i - 1) * dv
        ax.plot(t[:len(xx)], xx, label=f'stimulus {i + 1:d}')
    ax.axvline(on_time * 1000, color='gray', label='stimulus off')
    ax.set_xlim(0, 1000 * max(1.25 * t_cycle, t_start))
    ax.set_ylim(np.array((-2, n_cycle + 2)) * dv)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel(var_label)
    ax.legend(loc='lower right', frameon=False)


def fit_fooof(f, pxx, aperiodic_mode='fixed', dB_threshold=3., max_n_peaks=10,
              freq_range=None, peak_width_limits=None, report=False,
              plot=False, plt_log=False, plt_range=None, figsize=None):
    from fooof import FOOOF

    if aperiodic_mode != 'knee':
        aperiodic_mode = 'fixed'
    def set_range(x, upper=f[-1]):        
        x = np.array(upper) if x is None else np.array(x)
        return [f[2], x.item()] if x.size == 1 else x.tolist()
    freq_range = set_range(freq_range)
    peak_width_limits = set_range(peak_width_limits, np.inf)

    # Initialize a FOOOF object
    fm = FOOOF(peak_width_limits=peak_width_limits, min_peak_height=dB_threshold / 10,
               peak_threshold=0., max_n_peaks=max_n_peaks, aperiodic_mode=aperiodic_mode)
    # Fit the model
    fm.fit(f, pxx, freq_range)
    results = fm.get_results()

    if report:
        fm.print_results()
        if aperiodic_mode=='knee':
            ap_params = results.aperiodic_params
            if ap_params[1] <= 0:
                print('Negative value of knee parameter occurred. Suggestion: Fit without knee parameter.')
            knee_freq = np.abs(ap_params[1]) ** (1 / ap_params[2])
            print(f'Knee location: {knee_freq:.2f} Hz')
    if plot:
        plt_range = set_range(plt_range)
        fm.plot(plt_log=plt_log)
        plt.xlim(np.log10(plt_range) if plt_log else plt_range)
        if figsize:
            plt.gcf().set_size_inches(figsize)
        plt.show()
    return results, fm


def psd_residual(f, pxx, fooof_result, plot=False, plt_log=False, plt_range=None, ax=None):
    from fooof.sim.gen import gen_model

    full_fit, _, ap_fit = gen_model(f[1:], fooof_result.aperiodic_params,
                                    fooof_result.gaussian_params, return_components=True)
    full_fit, ap_fit = 10 ** full_fit, 10 ** ap_fit

    res_psd = np.insert(pxx[1:] - ap_fit, 0, 0.)
    res_fit = np.insert(full_fit - ap_fit, 0, 0.)

    if plot:
        if ax is None:
            _, ax = plt.subplots(1, 1)
        plt_range = np.array(f[-1]) if plt_range is None else np.array(plt_range)
        if plt_range.size == 1:
            plt_range = [f[1] if plt_log else 0., plt_range.item()]
        f_idx = (f >= plt_range[0]) & (f <= plt_range[1])
        ax.plot(f[f_idx], res_psd[f_idx], 'b', label='residual')
        ax.plot(f[f_idx], res_fit[f_idx], 'r', label='fit')
        if plt_log:
            ax.set_xscale('log')
        ax.set_xlim(plt_range)
        ax.legend(loc='upper right')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('PSD Residual')
    return res_psd, res_fit


def plot_spectrogram(x, fs, tseg, tres=np.inf, remove_aperiodic=None,
                     plt_log=False, plt_range=None, ax=None):
    nfft = int(tseg * fs) # steps per segment
    noverlap = int(max(tseg - tres, 0) * fs)
    f, t, sxx = ss.spectrogram(x, fs=fs, nperseg=nfft, noverlap=noverlap, window='boxcar')

    cbar_label = 'PSD' if remove_aperiodic is None else 'PSD Residual'
    if plt_log:
        with np.errstate(divide='ignore'):
            sxx = np.log10(sxx)
        cbar_label += ' log(power)'

    if remove_aperiodic is not None:
        from fooof.sim.gen import gen_aperiodic

        ap_fit = gen_aperiodic(f[1:], remove_aperiodic.aperiodic_params)
        sxx[1:, :] -= (ap_fit if plt_log else 10 ** ap_fit)[:, None]
        sxx[0, :] = 0.

    if ax is None:
        _, ax = plt.subplots(1, 1)
    plt_range = np.array(f[-1]) if plt_range is None else np.array(plt_range)
    if plt_range.size == 1:
        plt_range = [f[1] if plt_log else 0., plt_range.item()]
    f_idx = (f >= plt_range[0]) & (f <= plt_range[1])

    pcm = ax.pcolormesh(t * 1000., f, sxx, shading='gouraud')
    ax.set_ylim(plt_range)
    plt.colorbar(mappable=pcm, ax=ax, label=cbar_label)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Frequency (Hz)')

    return f, t, sxx


def plot(choose, spike_file=None, config=None, figsize=(6.4, 4.8)):

    global CONFIG
    if config is None:
        config = CONFIG
    else:
        CONFIG = config
    config = os.path.join(MODEL_PATH, config)
    if spike_file is None:
        spike_file = 'output/spikes.h5'

    if choose<=2:
        from bmtool.util import util
        
        nodes = util.load_nodes_from_config(config)
        network_name = 'BLA'
        cortex_df = nodes[network_name]

        with h5py.File(spike_file) as f:
            spikes_df = pd.DataFrame({
                'node_ids': f['spikes'][network_name]['node_ids'],
                'timestamps': f['spikes'][network_name]['timestamps']
            })
            spikes_df.sort_values(by='timestamps', inplace=True, ignore_index=True)

    if choose==1:
        spikes_df['pop_name'] = cortex_df.loc[spikes_df['node_ids'], 'pop_name'].values
        pop_spike = get_populations(spikes_df, pop_names)

        print("Plotting cortex spike raster")
        _, ax = plt.subplots(1, 1, figsize=figsize)
        raster(pop_spike, pop_color, ax=ax)
        plt.show()
        #_ = plot_raster(config_file=config, group_by='pop_name')
        return pop_spike

    elif choose==2:
        from bmtk.simulator import bionet

        conf = bionet.Config.from_json(config)
        t_stop = conf['run']['tstop'] / 1000

        frs = firing_rate(spikes_df, num_cells=len(cortex_df), time_windows=(0., t_stop))
        cortex_nodes = get_populations(cortex_df, pop_names, only_id=True)
        pop_fr = {p: frs[nid] for p, nid in cortex_nodes.items()}

        print('Firing rate: mean, std')
        for p, fr in pop_fr.items():
            print(f'{p}: {fr.mean():.3g}, {fr.std():.3g}')

        print("Plotting firing rates")
        min_fr = 0.5 / total_duration((0., t_stop)) # to replace 0 spikes
        _, ax = plt.subplots(1, 1, figsize=figsize)
        firing_rate_histogram(pop_fr, pop_color, bins=20, min_fr=min_fr,
                              logscale=True, stacked=False, ax=ax)
        plt.show()
        return frs

    elif choose==3:
        from bmtk.analyzer.compartment import plot_traces
        print("plotting voltage trace")
        _ = plot_traces(config_file=config, report_name='v_report',
                        group_by='pop_name', average=True)
    else:
        print('Choice ID not identified.')


if __name__ == '__main__':
    argv = sys.argv[1:]
    narg = len(argv)
    if narg > 0:
        argv[0] = int(argv[0])
    plot(*argv)