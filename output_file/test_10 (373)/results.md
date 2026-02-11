Found Sideband Peak at: x=1857 y=1102

=== step 3: ASM Propagation ===
Applying Gaussian smoothing (σ=2.5)...

=== Step 4 & 5: Unwrap & Scaling ===
Unwrapping Phase (Simple 2D)...
Applying Plane Subtraction (Tilt Removal)...
Max Height/Depth: 2716.0897621646272 nm
Min Height/Depth: -2607.4811354226485 nm

Done! Press Enter to exit...

=== PROCESS COMPLETE ===
PS D:\HoloScratch-Julia> julia --project=. holoyscratch-deterction.jl
=== System Configuration ===
Wavelength: 632.0 nm
Distance:   11.0 cm
Pixel Size: 3.45 μm

Image found: data-test/373.png

=== step 1: Setup & Load ===
Resolution: 3840 x 2160 pixels
Bit Depth:  Converted to 64-bit Float for processing

=== step 2: FFT & Spectral Analysis ===
Found Sideband Peak at: x=1857 y=1102

=== step 3: ASM Propagation ===

--- Auto-Focus Sweep ---
Coarse best z = 5.0 cm
Fine   best z = 4.0 cm  (score=0.349643)
Using z = 4.0 cm for reconstruction
Applying Gaussian smoothing (σ=2.5)...

=== Step 4 & 5: Unwrap & Scaling ===
Unwrapping Phase (Simple 2D)...
Applying Plane Subtraction (Tilt Removal)...

Focused at z = 4.0 cm
Max Height/Depth: 2090.8670024423263 nm
Min Height/Depth: -2454.152861283233 nm

Done! Press Enter to exit...

=== PROCESS COMPLETE ===