# Discrete Fourier Transform (DFT) Analysis

## AIM
Perform various operations on an image using Discrete Fourier Transform (DFT):
1. Implement DFT manually using the mathematical equation.
2. Verify the manual DFT implementation with the inbuilt function (e.g., `np.fft.fft2`, `np.fft.fftshift`, `cv2.dft()`, `cv2.idft()`).
3. Compute and plot the magnitude and phase spectrum.
4. Apply frequency shifting and observe the changes in magnitude and phase spectrum.
5. Apply image rotation and observe the corresponding changes in magnitude and phase spectrum.
6. Use two images to reconstruct new images by swapping magnitude and phase information.

## Requirements
- MATLAB
- Basic understanding of image processing and Fourier Transform

## Implementation Details

### 1. **Manual DFT Implementation**
- Implemented a function `dft_manual(image)` that computes the DFT using the mathematical equation.
- Applied this function to an 8x8 checkerboard image.

### 2. **Verification of Manual DFT with Inbuilt Function**
- Used `fft2(image)` to compute DFT and compared the results with the manual implementation.
- Used `isequal` to check accuracy.

### 3. **Magnitude and Phase Spectrum Calculation**
- Magnitude spectrum calculated as `20 * log10(abs(DFT) + 1)`.
- Phase spectrum obtained using `angle(DFT)`.
- Displayed using `imshow()`.

### 4. **Frequency Shifting**
- Applied `fftshift()` to bring the zero frequency component to the center.
- Observed the impact on magnitude and phase.

### 5. **Effect of Rotation on Frequency Spectrum**
- Rotated the image using `rot90()` and computed the corresponding DFT.
- Observed how the magnitude and phase spectrum change upon rotation.

### 6. **Image Reconstruction by Swapping Magnitude and Phase**
- Used a second image (a circular pattern) and computed its DFT.
- Swapped magnitude and phase components between the two images.
- Reconstructed images using `ifft2()` and displayed results.

## Results
- Verified that manual DFT matches the inbuilt function.
- Observed the effect of shifting and rotation on magnitude and phase spectra.
- Reconstructed meaningful images by swapping magnitude and phase between two images.

## Usage
1. Run the MATLAB script in a MATLAB environment.
2. Observe the outputs in the generated figures.

## References
- MATLAB Documentation on `fft2()`, `ifft2()`, `fftshift()`.
- Signal and Image Processing Textbooks for Fourier Transform concepts.

