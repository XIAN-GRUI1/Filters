[el, info1] = Butterworth(15, 2*10e8, 'lowpass', 'Steps', 2, 'Metric', 'max', 'Verbose', true);
[~, info2] = Chebyshev(33, 1e8, 'lowpass', 'Ripple', 0.5, 'Verbose', true);