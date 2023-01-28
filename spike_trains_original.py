from bmtk.utils.reports.spike_trains import SpikeTrains
from bmtk.utils.sonata.utils import get_node_ids
import numpy as np
from six import string_types

class PoissonSpikeGenerator(SpikeTrains):
    """ A Class for generating spike-trains with a homogeneous and inhomogeneous Poission distribution.
    Uses the methods describe in Dayan and Abbott, 2001.
    """
    max_spikes_per_node = 10000000

    def __init__(self, population=None, seed=None, output_units='ms', **kwargs):
        if population is not None and 'default_population' not in kwargs:
            kwargs['default_population'] = population

        if seed:
            np.random.seed(seed)

        super(PoissonSpikeGenerator, self).__init__(units=output_units, **kwargs)
        # self.units = units
        if output_units.lower() in ['ms', 'millisecond', 'milliseconds']:
            self._units = 'ms'
            self.output_conversion = 1000.0
        elif output_units.lower() in ['s', 'second', 'seconds']:
            self._units = 's'
            self.output_conversion = 1.0
        else:
            raise AttributeError('Unknown output_units value {}'.format(output_units))

    def add(self, node_ids, firing_rate, population=None, times=(0.0, 1.0)):
        # TODO: Add refactory period
        if isinstance(node_ids, string_types):
            # if user passes in path to nodes.h5 file count number of nodes
            node_ids = get_node_ids(node_ids, population)
        if np.isscalar(node_ids):
            # In case user passes in single node_id
            node_ids = [node_ids]

        if np.isscalar(firing_rate):
            self._build_fixed_fr(node_ids, population, firing_rate, times)
        elif isinstance(firing_rate, (list, np.ndarray)):
            self._build_inhomegeous_fr(node_ids, population, firing_rate, times)

    def time_range(self, population=None):
        df = self.to_dataframe(populations=population, with_population_col=False)
        timestamps = df['timestamps']
        return np.min(timestamps), np.max(timestamps)

    def _build_fixed_fr(self, node_ids, population, fr, times):
        if np.isscalar(times) and times > 0.0:
            tstart = 0.0
            tstop = times
        else:
            tstart = times[0]
            tstop = times[-1]
            if tstart >= tstop:
                raise ValueError('Invalid start and stop times.')

        for node_id in node_ids:
            c_time = tstart
            while True:
                interval = -np.log(1.0 - np.random.uniform()) / fr
                c_time += interval
                if c_time > tstop:
                    break

                self.add_spike(node_id=node_id, population=population, timestamp=c_time*self.output_conversion)

    def _build_inhomegeous_fr(self, node_ids, population, fr, times):
        if np.min(fr) <= 0:
            raise Exception('Firing rates must not be negative')
        max_fr = np.max(fr)

        times = times
        tstart = times[0]
        tstop = times[-1]

        for node_id in node_ids:
            c_time = tstart
            time_indx = 0
            while True:
                # Using the prunning method, see Dayan and Abbott Ch 2
                interval = -np.log(1.0 - np.random.uniform()) / max_fr
                c_time += interval
                if c_time > tstop:
                    break

                # A spike occurs at t_i, find index j st times[j-1] < t_i < times[j], and interpolate the firing rates
                # using fr[j-1] and fr[j]
                while times[time_indx] <= c_time:
                    time_indx += 1

                fr_i = _interpolate_fr(c_time, times[time_indx-1], times[time_indx],
                                       fr[time_indx-1], fr[time_indx])

                if not fr_i/max_fr < np.random.uniform():
                    self.add_spike(node_id=node_id, population=population, timestamp=c_time*self.output_conversion)


def _interpolate_fr(t, t0, t1, fr0, fr1):
    # Used to interpolate the firing rate at time t from a discrete list of firing rates
    return fr0 + (fr1 - fr0)*(t - t0)/(t1 - t0)