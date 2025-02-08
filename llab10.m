%% LAB:10 
% AIM: Consider an image and perform the following operations: 
% 1. Discrete Fourier Transform (DFT) using Equation 
% 2. Verify the implemented function of DFT with inbuilt function (Hint: np.fft.fft2, np.fft.fftshift or cv2.dft(), cv2.idft()) 
% 3. Plot magnitude and Phase spectrum (Hint: magnitude : 20*np.log(np.abs()), Phase: np.angle()) 4. Apply shifting operation and observe the magnitude and phase spectrum 
% 5. Apply rotation and observe the magnitude and phase spectrum (Comment on problem 4 and 5) (Hint: np.rot90()) 
% 6. Consider two images. Find the magnitude and phase spectrum of both images. Reconstruct the image using magnitude of first image and phase of second image and vice versa.


clc;
clear all;
close all;
datetime 


function F = dft_manual(image)
    [M, N] = size(image);
    F = zeros(M, N);
    for u = 1:M
        for v = 1:N
            sum_val = 0;
            for x = 1:M
                for y = 1:N
                    angle = -2 * pi * ((u * x / M) + (v * y / N));
                    sum_val = sum_val + image(x, y) * exp(1j * angle);
                end
            end
            F(u, v) = sum_val;
        end
    end
end



% Read and convert image to grayscale if not already
image = checkerboard(8) > 0.5; % 8x8 checker squares, creating a small binary checkerboard

if size(image, 3) == 3
    image = rgb2gray(image);
end
image = double(image) * 255; % Scale to 255 for grayscale compatibility

% Manual DFT
dft_manual_result = dft_manual(image);

% Inbuilt DFT
dft_builtin = fft2(image);

% Verify the results
comparison = isequal(round(dft_manual_result, 4), round(dft_builtin, 4));
disp("Manual DFT and Inbuilt DFT are the same: " + comparison);



% Compute magnitude and phase spectra
magnitude_spectrum = 20 * log10(abs(dft_builtin) + 1); % Adding 1 to avoid log(0)
phase_spectrum = angle(dft_builtin);

% Plot magnitude and phase spectrum
figure;
subplot(1, 2, 1);
imshow(magnitude_spectrum, []);
title('Magnitude Spectrum');
subplot(1, 2, 2);
imshow(phase_spectrum, []);
title('Phase Spectrum');


% Shift zero frequency component to center
shifted_dft = fftshift(dft_builtin);
shifted_magnitude = 20 * log10(abs(shifted_dft) + 1);
shifted_phase = angle(shifted_dft);

% Plot shifted magnitude and phase spectra
figure;
subplot(1, 2, 1);
imshow(shifted_magnitude, []);
title('Shifted Magnitude');
subplot(1, 2, 2);
imshow(shifted_phase, []);
title('Shifted Phase');




% Rotate image by 90 degrees
rotated_image = rot90(image);
rotated_dft = fft2(rotated_image);
rotated_magnitude = 20 * log10(abs(fftshift(rotated_dft)) + 1);
rotated_phase = angle(fftshift(rotated_dft));

% Plot rotated magnitude and phase spectra
figure;
subplot(1, 2, 1);
imshow(rotated_magnitude, []);
title('Rotated Magnitud');
subplot(1, 2, 2);
imshow(rotated_phase, []);
title('Rotated Phase');



% Load second image
% Use a different generated pattern for image2
image2 = zeros(64, 64);           % Create a 64x64 blank image
radius = 10;
center = [32, 32];
[x, y] = meshgrid(1:64, 1:64);
mask = sqrt((x - center(1)).^2 + (y - center(2)).^2) <= radius;
image2(mask) = 255;               % Set pixels inside the circle to white


if size(image2, 3) == 3
    image2 = rgb2gray(image2);
end
image2 = double(image2);

% Compute DFTs
dft_image1 = fft2(image);
dft_image2 = fft2(image2);

% Magnitude and phase of both images
magnitude1 = abs(dft_image1);
phase1 = angle(dft_image1);
magnitude2 = abs(dft_image2);
phase2 = angle(dft_image2);

% Reconstruct images by swapping magnitude and phase
reconstructed1 = real(ifft2(magnitude1 .* exp(1j * phase2)));
reconstructed2 = real(ifft2(magnitude2 .* exp(1j * phase1)));

% Plot reconstructed images
figure;
subplot(2, 1, 1);
imshow(reconstructed1, []);
title('Magnitude from Image 1, Phase from Image 2');
subplot(2, 1, 2);
imshow(reconstructed2, []);
title('Magnitude from Image 2, Phase from Image 1');
