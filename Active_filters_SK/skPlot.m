function skPlot(F)
%SKPLOT  Plot the designed filter response (actual vs ideal).
%   SKPLOT(F)  draws |H| in dB and the phase of the realized filter
%   together with the ideal prototype response.
%
%   See also SKDESIGN.

f = F.fc;
w = 2*pi*logspace(log10(f*1e-4), log10(f*1e4), 2001);
[H, Hi] = skResponse(F, w);
HdB = 20*log10(max(abs(H), 1e-12));
HidB = 20*log10(max(abs(Hi), 1e-12));
fv = w / (2*pi);

figure('Name', sprintf('SK filter: %s-%d %s', upper(F.type), F.order, upper(F.passtype)));
subplot(2, 1, 1);
semilogx(fv, HdB, 'LineWidth', 1.5); hold on;
semilogx(fv, HidB, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
grid on;
xlabel('Frequency (Hz)'); ylabel('|H| (dB)');
legend('realized (finite GBW)', 'ideal prototype', 'Location', 'southwest');
title(sprintf('%s order %d %s, fc = %s Hz%s', upper(F.type), F.order, ...
    upper(F.passtype), skFmtSI(f), ...
    sprintf(', GBW = %s Hz', skFmtSI(F.gbw/(2*pi)))));
ylim([min(-80, min(HdB)-10), max(5, max(HdB)+5)]);
xline(f, ':', 'fc');

subplot(2, 1, 2);
semilogx(fv, unwrap(angle(H))*180/pi, 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
xline(f, ':', 'fc');
end
