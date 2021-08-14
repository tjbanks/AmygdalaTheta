from bmtk.analyzer.spike_trains import plot_raster
import sys

def plot_rast(config_file):
    _ = plot_raster(config_file=config_file)

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        plot_rast(sys.argv[-1])
    else:
        plot_rast('simulation_configECP.json')
