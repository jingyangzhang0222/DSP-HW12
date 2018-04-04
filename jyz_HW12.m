close all
clear

% load image and basic process ¡ý
img=imread('daytona_512.png');
[height, width] = size(img);
gimg = img; 
dimg = double(gimg);

% set basic parameters¡ý
wp= 0.15 * pi; fp = wp / pi;
ws= 0.2 * pi;  fs = ws / pi;
wc = (wp + ws) / 2; fc = wc / pi;
delp = 0.05; Rp = -20 * log10(1 - delp);% delp : pass-band ripple
dels = 0.05; Rs = -20 * log10(dels);    % dels : stop-band ripple

[NIIR, nwp] = ellipord(fp, fs, Rp, Rs);
[bIIR, aIIR] = ellip(NIIR, Rp, Rs, nwp);
[HIIR, ~] = freqz(bIIR, aIIR);

[NFIR, f0, m0, w] = firpmord([wp ws] / pi,[1 0], [delp, dels]);
bFIR = firpm(NFIR, f0, m0, w);
[HFIR, ~]=freqz(bFIR, 1);

del = 0.0001; % a small number
rts = roots(bFIR);
inrts = rts(abs(rts)<=(1+del));%INSIDE or ON the U.C.
outrts = rts(abs(rts)>(1+del));%OUTSIDE the U.C.
reflectedrts = 1 ./ conj(outrts);
newrts = [inrts; reflectedrts];
newrts = leja(newrts); %Reorder roots using Leja to make poly work more stable
[minbFIR] = poly(newrts); %convert roots to polynomical coefficients
minbFIR = minbFIR * sum(bFIR) / sum(minbFIR);
[HminFIR, om] = freqz(minbFIR, 1);

disp(['The order of IIR filter: ', num2str(NIIR)])
disp(['The order of FIR filter: ', num2str(NFIR)])

% apply filters to each row of the image
IIRimg = zeros(height, width);      %initialization
IIRffimg = zeros(height, width);    %initialization
FIRimg = zeros(height, width);      %initialization
FIRffimg = zeros(height, width);    %initialization
minFIRimg = zeros(height, width);   %initialization
minFIRffimg = zeros(height, width); %initialization
for i = 1 : height
    IIRimg(i, :) = filter(bIIR, aIIR, dimg(i, :));
    IIRffimg(i, :) = filtfilt(bIIR, aIIR, dimg(i, :));
    FIRimg(i, :) = filter(bFIR, 1, dimg(i, :));
    FIRffimg(i, :) = filtfilt(bFIR, 1, dimg(i, :));
    minFIRimg(i, :) = filter(minbFIR,1, dimg(i, :));
    minFIRffimg(i, :) = filtfilt(minbFIR,1, dimg(i, :));
end

figure (1)
subplot(2, 2, 1);
plot(om/pi, abs(HIIR), 'r', 'Linewidth', 1)
xlabel('\omega/\pi');ylabel('|H(\omega)|')
ylim([0 1.1]);grid
title('Magnitude Response of IIR Ellipitic Filter')
subplot(2, 2, 2);
plot(om/pi, abs(HFIR), 'g', 'Linewidth', 1)
xlabel('\omega/\pi');ylabel('|H(\omega)|')
ylim([0 1.1]);grid
title('Magnitude Response of FIR Filter')
subplot(2, 2, 3);
plot(om/pi,abs(HminFIR), 'b', 'Linewidth', 1)
xlabel('\omega/\pi');ylabel('|H(\omega)|')
ylim([0 1.1]);grid
title('Magnitude Response of Minimum Phase FIR Filter')
subplot(2, 2, 4);
plot(om/pi, abs(HIIR),':r', om/pi, abs(HFIR),'-g', om/pi, abs(HminFIR), '--b', 'Linewidth', 1)
xlabel('\omega/\pi')
ylim([0 1.1]);grid
legend('IIR Ellipitic', 'FIR Filter', 'Minimum Phase FIR')
title('Magnitude Responses')

figure (2)
subplot(2, 2, 1);
plot(om / pi, unwrap(angle(HIIR) * 2) / 2, 'r', 'Linewidth', 1)
xlabel('\omega/\pi')
ylabel('\angleH(\omega)');grid
title('Phase Response of IIR Ellipitic Filter ')
subplot(2, 2, 2);
plot(om / pi, unwrap(angle(HFIR) * 2) / 2, 'g', 'Linewidth', 1)
xlabel('\omega/\pi')
ylabel('\angleH(\omega)');grid
title('Phase Response of FIR Filter ')
subplot(2, 2, 3);
plot(om / pi, unwrap(angle(HminFIR) * 2) / 2, 'b', 'Linewidth', 1)
xlabel('\omega/\pi')
ylabel('\angleH(\omega)');grid
title('Phase Response of Minimum Phase FIR Filter ')
subplot(2, 2, 4);
plot(om/pi,unwrap(angle(HIIR)*2)/2,':r',om/pi,unwrap(angle(HFIR)*2)/2,'-g',om/pi,unwrap(angle(HminFIR)*2)/2,'--b', 'Linewidth', 1);
xlabel('\omega/\pi');grid
legend('IIR Ellipitic', 'FIR Filter', 'Minimum Phase FIR')
title('Phase Responses')

figure (3)
subplot(2, 2, 1)
zplane(bIIR, aIIR);
title('Zero-Pole Plot of IIR Elliptic Filter')
subplot(2, 2, 2)
zplane(bFIR, 1);
title('Zero-Pole Plot of FIR Filter')
subplot(2, 2, 3)
zplane(minbFIR, 1);
title('Zero-Pole Plot of Minimum Phase FIR Filter')
subplot(2, 2, 4)
zplane(minbFIR, 1);
title('Minimum Phase FIR Filter (Zoomed in)')
axis([-.1 1.1 -.5 .5])

figure (4)
imshow(img, []);
title('Original Image');

figure (5)
subplot(1, 2, 1)
imshow(IIRimg, []);
title('IIR Filtered Image');
subplot(1, 2, 2)
imshow(IIRffimg, []);
title('IIR F-B Filtered Image');

figure (6)
subplot(1, 2, 1)
imshow(FIRimg, []);
title('FIR Filtered Image');
subplot(1, 2, 2)
imshow(FIRffimg, []);
title('FIR  F-B Filtered Image');

figure (7)
subplot(1, 2, 1)
imshow(minFIRimg, []);
title('Minimum Phase FIR Filtered Image')
subplot(1, 2, 2)
imshow(minFIRffimg, []);
title('Minimum Phase FIR F-B Filtered Image')