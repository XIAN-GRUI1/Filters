[el, info1] = Butterworth(5, 8e8, 'lowpass', 'Steps', 2, 'Metric', 'max', 'Verbose', true);
[~, info2] = Chebyshev(5, 8e8, 'lowpass', 'Ripple', 0.5, 'Verbose', true);