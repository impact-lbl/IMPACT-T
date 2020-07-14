import matplotlib
import matplotlib.pyplot
import numpy
import sys
import pathlib

def get_input_filename(bunch):
    """Return the filename of the input file for a particular bunch."""
    if bunch == 1:
        filename = 'ImpactT.in'
    else:
        filename = f'ImpactT{bunch}.in'
    if pathlib.Path(filename+'.rendered').is_file():
        filename += '.rendered'
    return filename

def get_bunch_count():
    """Get the number of bunches from the first input file."""
    input = read_input_file(get_input_filename(1))
    return int(input[1].split()[2])

def get_lattice():
    """Get the lattice from the first input file as a list of lists."""
    input = read_input_file(get_input_filename(1))
    return [line.split() for line in input[9:]]

def get_bpms(lattice, z_offset=None):
    """Get the location and file number of all BPMs as a list of tuples."""
    return [(get_position_as_text(float(elem[4]), z_offset), int(elem[2]))
            for elem in lattice if elem[3]=='-2']

def get_position_as_text(position, z_offset=None):
    """Return a text description of a numerical position."""
    if z_offset:
        position = position - z_offset
    if position > 1:
        unit = 'm'
    else:
        position = position*1000
        unit = 'mm'
    if position >= 1000:
        position = round(position)
    else:
        position = float(f'{position:.4g}')
    return f'{position} {unit}'

def get_bunch_counts(bunch_list):
    """Get the particle counts for each bunch as an array."""
    input_filenames = [get_input_filename(bunch) for bunch in bunch_list]
    counts = [get_particle_count(filename) for filename in input_filenames]
    return numpy.array(counts)

def get_particle_count(filename):
    """Get the macroparticle count from a given input file."""
    input = read_input_file(filename)
    return int(input[2].split()[1])

def get_mass(filename):
    """Get the particle mass from a given input file."""
    input = read_input_file(filename)
    return float(input[8].split()[2])

def is_mass_matched(bunch_list):
    """Check whether the particle mass values are the same for given bunches."""
    mass = [get_mass(get_input_filename(bunch)) for bunch in bunch_list]
    return len(set(mass)) == 1

def bunch_text(bunch_list):
    """Create a text string to describe the bunch list for titles."""
    if get_bunch_count() == 1:
        return 'single bunch'
    elif len(bunch_list) == 1:
        return f'bunch {bunch_list[0]}'
    elif len(bunch_list) == get_bunch_count():
        return 'all bunches'
    elif bunch_list[-1] - bunch_list[0] + 1 == len(bunch_list):
        return f'bunches {bunch_list[0]} to {bunch_list[-1]}'
    else:
        return ('bunches '
                + ', '.join([str(bunch) for bunch in bunch_list[:-1]])
                + ' and '
                + str(bunch_list[-1]))

def check_bunch_list(bunch_list):
    """Remove any bunches from the list that don't exist."""
    max_bunch = get_bunch_count()
    valid_bunches = [bunch for bunch in bunch_list if bunch <= max_bunch]
    invalid_bunches = [bunch for bunch in bunch_list if bunch > max_bunch]
    return valid_bunches, invalid_bunches

def read_input_file(filename):
    """Read input file and return list of input lines, excluding comments."""
    with open(filename, 'r') as f:
        return [line.strip() for line in f.readlines() if line[0] != '!']

def load_data_from_file(filename, header=0):
    """Load data from a text file, skipping header rows, and return an array."""
    with open(filename, 'r') as f:
        data = [[float(value) for value in line.split()]
                for line in f.readlines()[header:]]
    return numpy.array(data)

def load_experimental_results(filename='experimental_data.txt'):
    """Load z, rms beam size and error values (m) from file as an array."""
    return load_data_from_file(filename, header=1)

def load_statistics_data(bunch_list, z_offset=None):
    """Load x and y data per bunch from statistics files as arrays."""
    xdata = []
    ydata = []
    for bunch in bunch_list:
        x = load_data_from_file(f'fort.{bunch}024')
        y = load_data_from_file(f'fort.{bunch}024')
        if z_offset:
            x = shift_z(x, z_offset, 1)
            y = shift_z(y, z_offset, 1)
        xdata.append(x)
        ydata.append(y)
    return numpy.array(xdata), numpy.array(ydata)

def load_bunch_count_data(z_offset=None):
    """Load bunch counts vs time from `fort.11` as an array."""
    data = load_data_from_file('fort.11')
    data = numpy.delete(data, 0, 1) # remove first column (time-step #)
    if z_offset:
        data = shift_z(data, z_offset, 1)
    return data

def load_phase_space_data(filenumber, bunch_list, z_offset=None):
    """Load phase space data per bunch as a list of arrays."""
    data = []
    for bunch in bunch_list:
        this_data = load_data_from_file(f'fort.{filenumber+bunch-1}')
        if z_offset:
            this_data = shift_z(this_data, z_offset, 4)
        data.append(calculate_energies(this_data, bunch))
    return data

def shift_z(data, z_offset, z_col=1):
    """Shift all z values by the given offset."""
    data[:,z_col] = data[:,z_col] - z_offset
    return data

def calculate_energies(data, bunch):
    """Calculate the energies for per-bunch data and save as an extra column."""
    gamma = numpy.sqrt(1 + numpy.square(data.T[1])
                         + numpy.square(data.T[3])
                         + numpy.square(data.T[5]))
    mass = get_mass(get_input_filename(bunch))
    W = (gamma - 1)*mass
    new_data = numpy.zeros((data.shape[0], data.shape[1]+1))
    new_data[:,:-1] = data
    new_data[:,-1] = W
    return new_data

def combine_bunch_values(data, bunch_list):
    """Combine values of separate bunches into a single summary dataset."""
    if len(data) == 1:
        return data[0]
    npt = numpy.array(get_bunch_counts(bunch_list))
    # Decompose data
    t_data = data[:,:,0]
    z0_data = data[:,:,1]
    x0_data = data[:,:,2]
    xrms_data = data[:,:,3]
    px0_data = data[:,:,4]
    pxrms_data = data[:,:,5]
    xpx_data = -data[:,:,6]
    # All t data should be the same
    t = t_data[0]
    if not all([row == t.tolist() for row in t_data.tolist()]):
        raise ValueError('Time step values not in sync across bunches.')
    # Centroids can be combined by a weighted mean
    z0 = numpy.sum(z0_data.T.dot(numpy.diag(npt)).T, 0)/npt.sum()
    x0 = numpy.sum(x0_data.T.dot(numpy.diag(npt)).T, 0)/npt.sum()
    px0 = numpy.sum(px0_data.T.dot(numpy.diag(npt)).T, 0)/npt.sum()
    # To combine rms values we need the weighted sum of squares
    sumsqx = numpy.sum(numpy.square(xrms_data).T.dot(numpy.diag(npt)).T
                       + numpy.square(x0_data).T.dot(numpy.diag(npt)).T, 0)
    sumsqpx = numpy.sum(numpy.square(pxrms_data).T.dot(numpy.diag(npt)).T
                        + numpy.square(px0_data).T.dot(numpy.diag(npt)).T, 0)
    sumxpx = numpy.sum((xpx_data + x0_data*px0_data).T.dot(numpy.diag(npt)).T, 0)
    xrms = numpy.sqrt(sumsqx/npt.sum() - numpy.square(x0))
    pxrms = numpy.sqrt(sumsqpx/npt.sum() - numpy.square(px0))
    xpx = sumxpx/npt.sum() - (x0 * px0)
    epx = numpy.sqrt(xrms*xrms * pxrms*pxrms - xpx*xpx)
    # Return combined data as an array
    return numpy.array([t, z0, x0, xrms, px0, pxrms, -xpx, epx]).T

def combine_phase_space_data(data):
    """Combine per-bunch phase space data into single array."""
    return numpy.concatenate(data)

def get_xdata(data, xaxis):
    """Return data for the given x-axis: time t (ns) or location z (mm)."""
    if xaxis == 't':
        return data.T[0]*1.0e9
    elif xaxis == 'z':
        return data.T[1]*1.0e3
    else:
        raise ValueError(f'Invalid xaxis specifier: {xaxis}')

def get_xlabel(xaxis):
    """Return label for the given x-axis: time t (ns) or location z (mm)."""
    if xaxis == 't':
        return 'Time (ns)'
    elif xaxis == 'z':
        return 'Average z-position (mm)'
    else:
        raise ValueError(f'Invalid xaxis specifier: {xaxis}')

def plot_beam_size(axes, data, bunch_list, xaxis='z', title=None,
                   experiment_data=None, combined_data=None):
    """Create a plot of beam size per bunch."""
    if xaxis == 't' and experiment_data is not None:
        print('Cannot plot experimental results against time t.')
        print('Continuing without experimental results.')
        experiment_data = None
    if not title:
        title = 'Beam sizes for ' + bunch_text(bunch_list)
    axes.figure.suptitle(title)
    axes.set_xlabel(get_xlabel(xaxis))
    axes.set_ylabel('Beam size (mm)')
    if combined_data is None:
        combined_data = combine_bunch_values(data, bunch_list)
    if experiment_data is not None:
        plot_beam_size_experimental(axes, experiment_data)
    for i in range(len(bunch_list)):
        plot_beam_size_single(axes, data[i], xaxis, '--',
                              label=f'Bunch {bunch_list[i]} rms')
    if len(bunch_list) > 1:
        plot_beam_size_single(axes, combined_data, xaxis, 'r-',
                              label=f'Combined rms')
    axes.legend(fontsize='x-small')

def plot_beam_size_experimental(axes, data):
    """Plot experimental data points with error bars"""
    z = data.T[0]*1.0e3
    rms = data.T[1]*1.0e3
    error = data.T[2]*1.0e3
    axes.errorbar(z, rms, yerr=error,
                  fmt='ko', markersize=2.0, elinewidth=0.5, capsize=1.0,
                  label='Experimental rms')

def plot_beam_size_single(axes, data, xaxis, fmt, label):
    """Plot the rms beam size for the given data onto the given axes."""
    x = get_xdata(data, xaxis)
    rms = data.T[3]*1.0e3
    axes.plot(x, rms, fmt, linewidth=1, label=label)

def plot_emittance(axes, xdata, ydata, bunch_list, xaxis='t', title=None,
                   combined_xdata=None, combined_ydata=None):
    """Create a plot of average x and y emittance per bunch."""
    if not title:
        title = 'Emittance for ' + bunch_text(bunch_list)
    axes.figure.suptitle(title)
    axes.set_xlabel(get_xlabel(xaxis))
    axes.set_ylabel('Normalised rms emittance (π mm mrad)')
    for i in range(len(bunch_list)):
        plot_emittance_single(axes, xdata[i], ydata[i], xaxis,
                              '--', label=f'Bunch {bunch_list[i]}')
    if len(bunch_list) > 1:
        if combined_xdata is None:
            combined_xdata = combine_bunch_values(xdata, bunch_list)
        if combined_ydata is None:
            combined_ydata = combine_bunch_values(ydata, bunch_list)
        plot_emittance_single(axes, combined_xdata, combined_ydata, xaxis,
                              'r-', label='Combined')
    axes.legend(fontsize='x-small')

def plot_emittance_single(axes, xdata, ydata, xaxis, fmt, label):
    """Plot the average x and y emittance onto the given axes."""
    x = get_xdata(xdata, xaxis)
    emittance = (xdata.T[7] + ydata.T[7])/2 * 1e6
    axes.plot(x, emittance, fmt, linewidth=1, label=label)

def plot_emittance_growth(axes, xdata, ydata, bunch_list, xaxis='t', title=None,
                          combined_xdata=None, combined_ydata=None):
    """Create a plot of average x and y emittance growth per bunch."""
    if not title:
        title = 'Emittance growth for ' + bunch_text(bunch_list)
    axes.figure.suptitle(title)
    axes.set_xlabel(get_xlabel(xaxis))
    axes.set_ylabel('Average emittance growth in x and y (relative)')
    for i in range(len(bunch_list)):
        plot_emittance_growth_single(axes, xdata[i], ydata[i], xaxis,
                                     '--', label=f'Bunch {bunch_list[i]}')
    if len(bunch_list) > 1:
        if combined_xdata is None:
            combined_xdata = combine_bunch_values(xdata, bunch_list)
        if combined_ydata is None:
            combined_ydata = combine_bunch_values(ydata, bunch_list)
        plot_emittance_growth_single(axes, combined_xdata, combined_ydata, xaxis,
                                     'r-', label=f'Combined')
    axes.legend(fontsize='x-small')

def plot_emittance_growth_single(axes, xdata, ydata, xaxis, fmt, label):
    """Plot the average x and y emittance growth onto the given axes."""
    x = get_xdata(xdata, xaxis)
    emittance = (xdata.T[7] + ydata.T[7])/2 * 1e6
    initial_emittance = max(emittance[0], 1e-10)
    growth = emittance/initial_emittance - 1
    axes.plot(x, growth, fmt, linewidth=1, label=label)

def plot_bunch_count(axes, data, bunch_list, xaxis='t', title=None):
    """Plot the number of particles per bunch against 't' or 'z'."""
    if not title:
        title = 'Bunch counts for ' + bunch_text(bunch_list)
    axes.figure.suptitle(title)
    x = get_xdata(data, xaxis)
    xlabel = get_xlabel(xaxis)
    counts = data.T[3:][[bunch - 1 for bunch in bunch_list]]
    axes.stackplot(x, counts, labels=[f'Bunch {bunch}' for bunch in bunch_list])
    axes.set_xlim(left=0.0)
    axes.set_xlabel(xlabel)
    axes.set_ylabel('Number of macroparticles')
    axes.legend()

def plot_phase_spaces(axes, data, bunch_list, title=None, grid_size=100):
    """Plot four phase spaces onto the given array of axes."""
    if not title:
        title = 'Phase space for ' + bunch_text(bunch_list)
    axes[0,0].figure.suptitle(title)
    x = data.T[0]
    px = data.T[1]
    y = data.T[2]
    py = data.T[3]
    z = data.T[4]
    pz = data.T[5]
    W = data.T[6]
    xp = numpy.arctan(px/pz)
    yp = numpy.arctan(py/pz)
    plot_phase_space(axes[0,0], x*1e3, xp*1e3, 'x (mm)', 'x` (mrad)', grid_size)
    plot_phase_space(axes[0,1], y*1e3, yp*1e3, 'y (mm)', 'y` (mrad)', grid_size)
    plot_phase_space(axes[1,1], x*1e3, y*1e3, 'x (mm)', 'y (mm)', grid_size)
    plot_phase_space(axes[1,0], z*1e3, W/1e6, 'z (mm)', 'Energy (MeV)', grid_size)

def plot_phase_space(axes, xdata, ydata, xlabel, ylabel, grid_size=100):
    """Plot a single phase space onto the given axes."""
    if grid_size < 10:
        grid_size = 10
    axes.set_xlabel(xlabel, fontsize='x-small')
    axes.set_ylabel(ylabel, fontsize='x-small')
    axes.tick_params(labelsize='xx-small')
    hist2d = plot_phase_space_hist2d(axes, xdata, ydata, grid_size)
    add_plot_margins(axes, 0.1)
    plot_phase_space_hist1d(axes, hist2d, grid_size)

def plot_phase_space_hist2d(axes, x, y, grid_size=100):
    """Plot the 2d histogram part of the phase space plot."""
    colour_map = matplotlib.cm.get_cmap('jet')
    colour_map.set_under('white', 0.)
    return axes.hist2d(x, y, bins=grid_size, cmap=colour_map, cmin=1)

def plot_phase_space_hist1d(axes, hist2d, grid_size=100):
    """Plot 1d histograms on the axes of the phase space plot."""
    hist, xedges, yedges, _ = hist2d
    xmin, xmax, ymin, ymax = xedges[0], xedges[-1], yedges[0], yedges[-1]
    x0, x1, y0, y1 = axes.axis()
    xscale = numpy.array(range(grid_size)) / grid_size * (xmax - xmin) + xmin
    yscale = numpy.array(range(grid_size)) / grid_size * (ymax - ymin) + ymin
    xhist = numpy.nansum(hist, 1)
    yhist = numpy.nansum(hist, 0)
    xhist_scaled = xhist / xhist.max() * (y1 - y0) * 0.2 + y0
    yhist_scaled = yhist / yhist.max() * (x1 - x0) * 0.2 + x0
    axes.plot(xscale, xhist_scaled, color='green', linewidth=0.75)
    axes.plot(yhist_scaled, yscale, color='green', linewidth=0.75)

def plot_bunch_energies(axes, data, bunch_list, title=None, bins=100):
    """Plot per-bunch energy spectra histograms."""
    if not title:
        title = 'Energy spectra for ' + bunch_text(bunch_list)
    axes.figure.suptitle(title)
    W = [bunch.T[6]/1e6 for bunch in data]
    if len(bunch_list) > 1:
        axes.hist(numpy.concatenate(W), bins=bins, label='Total',
                  histtype='stepfilled', linewidth=1.0,
                  color='red', facecolor=(1,0,0,0.1), edgecolor=(1,0,0,1.0))
    axes.hist(W, bins=bins, histtype='stepfilled', alpha=0.5,
              label=[f'Bunch {bunch}' for bunch in bunch_list])
    axes.set_xlabel('Energy (MeV)')
    axes.set_ylabel('Number of macroparticles')
    handles, labels = axes.get_legend_handles_labels()
    handles.reverse()
    labels.reverse()
    axes.legend(handles, labels)

def plot_total_energy(axes, combined_data, bunch_list, title=None, bins=100):
    """Plot total energy spectrum histogram on log scale."""
    if not title:
        title = 'Total energy spectrum for ' + bunch_text(bunch_list)
    axes.figure.suptitle(title)
    W = combined_data.T[6]/1e6
    axes.hist(W, bins=bins, histtype='stepfilled', linewidth=1.0,
              color='red', facecolor=(1,0,0,0.1), edgecolor=(1,0,0,1.0))
    axes.set_xlabel('Energy (MeV)')
    axes.set_ylabel('Number of macroparticles')
    axes.set_yscale('log')

def add_plot_margins(axes, margin):
    """Adjust the axis limits to include a margin around the data."""
    xmin, xmax = axes.get_xlim()
    ymin, ymax = axes.get_ylim()
    xmargin = margin*(xmax - xmin)
    ymargin = margin*(ymax - ymin)
    axes.set_xlim(xmin - xmargin, xmax + xmargin)
    axes.set_ylim(ymin - ymargin, ymax + ymargin)

def plot_all(bunch_list, z_offset=None):
    """Run and save all plots consecutively."""
    full_bunch_list = list(range(1, int(get_bunch_count()) + 1))
    bunch_list, invalid_bunches = check_bunch_list(bunch_list)
    if invalid_bunches:
        print(f'Skipping invalid bunches: {invalid_bunches}')
    print('Loading experimental data...')
    try:
        experimental_results = load_experimental_results()
    except FileNotFoundError as err:
        print(f'Experimental results data file not found: {err}')
        print('Continuing without experimental data.')
        experimental_results = None
    print('Loading statistical data...')
    try:
        xdata, ydata = load_statistics_data(bunch_list, z_offset)
    except FileNotFoundError as err:
        print(f'Statistical data file not found: {err}')
        print('Skipping statistical plots.')
    else:
        combined_xdata = combine_bunch_values(xdata, bunch_list)
        combined_ydata = combine_bunch_values(ydata, bunch_list)
        print('Plotting beam size...')
        figure, axes = matplotlib.pyplot.subplots(dpi=300)
        plot_beam_size(axes, xdata, bunch_list, combined_data=combined_xdata)
        figure.savefig('beam-size')
        matplotlib.pyplot.close(figure)
        if experimental_results is not None:
            figure, axes = matplotlib.pyplot.subplots(dpi=300)
            plot_beam_size(axes, xdata, bunch_list,
                           combined_data=combined_xdata,
                           experiment_data=experimental_results)
            figure.savefig('beam-size-vs-experiment')
            matplotlib.pyplot.close(figure)
        print('Plotting emittance...')
        figure, axes = matplotlib.pyplot.subplots(dpi=300)
        plot_emittance(axes, xdata, ydata, bunch_list,
                       combined_xdata=combined_xdata,
                       combined_ydata=combined_ydata)
        figure.savefig('emittance')
        matplotlib.pyplot.close(figure)
        figure, axes = matplotlib.pyplot.subplots(dpi=300)
        plot_emittance_growth(axes, xdata, ydata, bunch_list,
                              combined_xdata=combined_xdata,
                              combined_ydata=combined_ydata)
        figure.savefig('emittance-growth')
        matplotlib.pyplot.close(figure)
    print('Loading bunch count data...')
    try:
        data = load_bunch_count_data(z_offset)
    except FileNotFoundError as err:
        print(f'Bunch count data file not found: {err}')
        print('Skipping bunch count plot.')
    else:
        print('Plotting bunch counts...')
        figure, axes = matplotlib.pyplot.subplots(dpi=300)
        plot_bunch_count(axes, data, bunch_list)
        figure.savefig('bunch-count')
        matplotlib.pyplot.close(figure)
    print('Loading initial phase space data...')
    try:
        full_data = load_phase_space_data(40, full_bunch_list, z_offset)
    except FileNotFoundError as err:
        print(f'Phase space data file not found: {err}')
        print('Skipping initial phase space step.')
    else:
        data = [full_data[bunch-1] for bunch in bunch_list]
        combined_data = combine_phase_space_data(data)
        print('Plotting initial phase space data...')
        figure, axes = matplotlib.pyplot.subplots(nrows=2, ncols=2, dpi=300)
        plot_phase_spaces(axes, combined_data, bunch_list, grid_size=300,
            title=f'Initial phase space for {bunch_text(bunch_list)}')
        figure.savefig('phase-space-initial')
        matplotlib.pyplot.close(figure)
        if len(full_bunch_list) > 1:
            for bunch in full_bunch_list:
                figure, axes = matplotlib.pyplot.subplots(2, 2, dpi=300)
                plot_phase_spaces(axes, full_data[bunch-1], [bunch],
                    title=f'Initial phase space for {bunch_text([bunch])}',
                    grid_size=300)
                figure.savefig(f'phase-space-initial-bunch{bunch}')
                matplotlib.pyplot.close(figure)
        print('Plotting initial energy spectra...')
        figure, axes = matplotlib.pyplot.subplots(dpi=300)
        plot_bunch_energies(axes, data, bunch_list, bins=300,
            title=f'Initial energy spectra for {bunch_text(bunch_list)}')
        figure.savefig('energies-initial')
        matplotlib.pyplot.close(figure)
        figure, axes = matplotlib.pyplot.subplots(dpi=300)
        plot_total_energy(axes, combined_data, bunch_list, bins=300,
            title=f'Initial total energy spectrum for {bunch_text(bunch_list)}')
        figure.savefig('energy-initial')
        matplotlib.pyplot.close(figure)
    print('Loading final phase space data...')
    try:
        full_data = load_phase_space_data(50, full_bunch_list, z_offset)
    except FileNotFoundError as err:
        print(f'Phase space data file not found: {err}')
        print('Skipping final phase space step.')
    else:
        data = [full_data[bunch-1] for bunch in bunch_list]
        combined_data = combine_phase_space_data(data)
        print('Plotting final phase space data...')
        figure, axes = matplotlib.pyplot.subplots(nrows=2, ncols=2, dpi=300)
        plot_phase_spaces(axes, combined_data, bunch_list, grid_size=300,
            title=f'Final phase space for {bunch_text(bunch_list)}')
        figure.savefig('phase-space-final')
        matplotlib.pyplot.close(figure)
        if len(full_bunch_list) > 1:
            for bunch in full_bunch_list:
                figure, axes = matplotlib.pyplot.subplots(2, 2, dpi=300)
                plot_phase_spaces(axes, full_data[bunch-1], [bunch],
                    title=f'Final phase space for {bunch_text([bunch])}',
                    grid_size=300)
                figure.savefig(f'phase-space-final-bunch{bunch}')
                matplotlib.pyplot.close(figure)
        print('Plotting final energy spectra...')
        figure, axes = matplotlib.pyplot.subplots(dpi=300)
        plot_bunch_energies(axes, data,  bunch_list, bins=300,
            title=f'Final energy spectra for {bunch_text(bunch_list)}')
        figure.savefig('energies-final')
        matplotlib.pyplot.close(figure)
        figure, axes = matplotlib.pyplot.subplots(dpi=300)
        plot_total_energy(axes, combined_data, bunch_list, bins=300,
            title=f'Final total energy spectrum for {bunch_text(bunch_list)}')
        figure.savefig('energy-final')
        matplotlib.pyplot.close(figure)
    print('Getting list of BPMs...')
    try:
        lattice = get_lattice()
        bpm_list = get_bpms(lattice, z_offset)
    except FileNotFoundError as err:
        print(f'Input file not found: {err}')
        print('Skipping BPM plot steps.')
    else:
        for location, filenumber in bpm_list:
            print(f'Loading BPM {filenumber} phase space data...')
            try:
                full_data = load_phase_space_data(filenumber, full_bunch_list,
                                                  z_offset)
            except FileNotFoundError as err:
                print(f'BPM data file not found: {err}')
                print('Skipping this BPM plot step.')
            else:
                data = [full_data[bunch-1] for bunch in bunch_list]
                combined_data = combine_phase_space_data(data)
                print(f'Plotting BPM {filenumber} phase space data...')
                figure, axes = matplotlib.pyplot.subplots(2, 2, dpi=300)
                plot_phase_spaces(axes, combined_data, bunch_list,
                                  title=(f'Phase space at z = {location} '
                                         f'for {bunch_text(bunch_list)}'),
                                  grid_size=300)
                figure.savefig(f'phase-space-{filenumber}')
                matplotlib.pyplot.close(figure)
                if len(full_bunch_list) > 1:
                    for bunch in full_bunch_list:
                        figure, axes = matplotlib.pyplot.subplots(2, 2, dpi=300)
                        plot_phase_spaces(axes, full_data[bunch-1], [bunch],
                                        title=(f'Phase space at z = {location} '
                                                f'for {bunch_text([bunch])}'),
                                        grid_size=300)
                        figure.savefig(f'phase-space-{filenumber}-bunch{bunch}')
                        matplotlib.pyplot.close(figure)
                print(f'Plotting BPM {filenumber} energy spectra...')
                figure, axes = matplotlib.pyplot.subplots(dpi=300)
                plot_bunch_energies(axes, data, bunch_list, bins=300,
                                    title=(f'Energy spectra at z = {location} '
                                           f'for {bunch_text(bunch_list)}'))
                figure.savefig(f'energies-{filenumber}')
                matplotlib.pyplot.close(figure)
                figure, axes = matplotlib.pyplot.subplots(dpi=300)
                plot_total_energy(axes, combined_data, bunch_list, bins=300,
                                  title=(f'Energy spectrum at z = {location} '
                                         f'for {bunch_text(bunch_list)}'))
                figure.savefig(f'energy-{filenumber}')
                matplotlib.pyplot.close(figure)

if __name__ == '__main__':
    matplotlib.use('agg') # Use the AGG renderer to produce PNG output
    if len(sys.argv) > 1:
        bunch_parameter = sys.argv[1]
        if bunch_parameter.isdigit():
            bunch_list = list(range(1, int(bunch_parameter) + 1))
        else:
            try:
                bunch_list = [int(value) for value in sys.argv[1].split(',')]
            except:
                bunch_list = list(range(1, int(get_bunch_count()) + 1))
    else:
        bunch_list = list(range(1, int(get_bunch_count()) + 1))
    if len(sys.argv) > 2:
        z_offset = sys.argv[2]
        try:
            z_offset = float(z_offset)
        except:
            z_offset = None
    else:
        z_offset = None
    plot_all(bunch_list, z_offset)
