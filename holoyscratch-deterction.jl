using Images, FileIO, Plots, FFTW, Statistics
gr()

const wavelength = 632e-9
const dx_pixel_size = 3.45e-6
const distance = 0.1
const dy_pixel_size = dx_pixel_size
const λ = wavelength
const z = distance

println("=== System Configuration ===")
println("Wavelength: $(wavelength * 1e9) nm")
println("Distance:   $(distance * 100) cm")
println("Pixel Size: $(dx_pixel_size * 1e6) μm")


# Load Data & Process
img_path = "data-test/373.png"

if isfile(img_path)
    println("\nImage found: $img_path")
    println("\n=== step 1: Setup & Load ===")
    img = load(img_path)
    img_gray = Gray.(img)
    u_hologram = Float64.(img_gray)

    (M, N) = size(u_hologram)
    println("Resolution: $N x $M pixels")
    println("Bit Depth:  Converted to 64-bit Float for processing")

    # Plot Raw Hologram
    p1 = heatmap(u_hologram, 
        title = "Raw Hologram",
        color = :gray,
        aspect_ratio = :equal,
        axis = nothing
    )
    
    println("\n=== step 2: FFT & Spectral Analysis ===")
    # --- FFT & Spectral Analysis ---
    F = fft(u_hologram)
    F_shifted = fftshift(F)
    spectrum_mag = log.(abs.(F_shifted) .+ 1)

    # plot Frequency Spectrum
    p2 = heatmap(spectrum_mag,
        c = :viridis,
        title = "Frequency Spectrum (Log)",
        aspect_ratio = :equal,
        axis = nothing
    )

    # Mark DC Term
    (h, w) = size(spectrum_mag)
    cx, cy = w/2 +1, h/2 +1
    plot!(p2, [cx], [cy], seriestype=:scatter, color=:red, label="DC Center", markersize=5)

    # === Spaatial Filtering ===
    cxx , cyy = Int(w/2), Int(h/2)
    search_spectrum = copy(spectrum_mag)

    # close ตรงกลางภาพ (DC term) เพื่อไม่ให้คอม เข้าใจผิดว่าจุดกลางคือ Sideband 
    mask_radius_dc = 50
    for i in 1:h, j in 1:w
        if (i - cyy)^2 + (j - cxx)^2 <mask_radius_dc^2
            search_spectrum[i, j] = 0
        end
    end

    # หา pos ที่สว่างที่สุดที่เหลืออยู่ (Sideband)
    max_val, max_idx = findmax(search_spectrum)
    peak_y,  peak_x = max_idx[1], max_idx[2]
    println("Found Sideband Peak at: x=$peak_x y=$peak_y")

    # creat image Circular filter Mask
    filter_redius = 150 # <--- ขนาดรูรับแสง (ถ้าภาพเบลอให้ปรับค่านี้)
    mask = zeros(ComplexF64, h, w)    
    
    for i in 1:h, j in 1:w
        dist = sqrt((i - peak_y)^2 + (j - peak_x)^2)

        # ให้เก็บค่าความถี่ไว้ ถ้าอยู่ในวง
        if dist < filter_redius
            mask[i, j] = 1.0
        end
    end

    # cute and moive to Center
    F_filtered = F_shifted .* mask
    shift_y = cyy - peak_y
    shift_x = cxx - peak_x
    F_centered = circshift(F_filtered, (shift_y, shift_x))

    # แปลงกลับเป็นคลื่น (Inverse FFT)
    U_filtered = ifft(ifftshift(F_centered))

    p3 = heatmap(
        log.(abs.(F_filtered) .+ 1),
        c=:viridis,
        title="Filtered sidedand",
        axis=nothing
    )

    p4 = heatmap(
        log.(abs.(F_centered) .+ 1),
        c=:viridis,
        title="Filtered sidedand",
        axis=nothing
    )
    
    println("\n=== step 3: ASM Propagation ===")
    # Creat Frequency Coordinates and Creat asxis Frequency fx, fy ให้ตรงกับภาพ
    (M, N) = size(U_filtered)
    fx = range(-1/(2*dx_pixel_size), 1/(2*dx_pixel_size), length=N)
    fy = range(-1/(2*dy_pixel_size), 1/(2*dy_pixel_size), length=M)

    # Creat Grid เพื่อให้คำนวณพิกเซลพร้อมกัน
    FX = repeat(reshape(fx, 1, :), M, 1)
    FY = repeat(reshape(fy, :, 1), 1, N)

    # Creat Transfer Function (H)
    k = 2 * π / λ
    root_term = Complex.(1 .- (λ .* FX).^2 .- (λ .* FY).^2)

    # สร้างแผ่นกรองแสง H
    H = exp.(1im .* k .* z .* sqrt.(root_term))
    # ตัดคลื่น Evanescent (คลื่นที่เดินทางไม่ได้) ทิ้งไป เพื่อลด Noise
    H[real.(root_term) .< 0] .= 0

    # ย้ายแสงไปที่วัตถุ (Propagation)
    A_0 = fftshift(fft(U_filtered))
    A_z = A_0 .* H

    # แปลงกลับเป็นภาพ (Inverse FFT)
    U_recon = ifft(ifftshift(A_z))

    # แยกข้อมูล Amplitude และ Phase
    amp_img = abs.(U_recon) # ความสว่าง (ไว้ดูความชัด)
    phase_img = angle.(U_recon) # เฟส (ไว้คำนวณความลึก)

    p5 = heatmap(amp_img,
        c=:gray,
        title="Reconstructed Amplitude (z=$(z*100) cm)",
        aspect_ratio=:equal,
        axis=nothing
    )

    p6 = heatmap(phase_img,
        c=:viridis,
        title="Reconstructed Phase (Wrapped)", 
        aspect_ratio=:equal,
        axis=nothing
    )


    # # show All Plots
    plt = plot(p5, p6, layout=(1, 2), size=(1000, 500))
    output_path = "ASM_Propagation.png"
    savefig(plt, output_path)

    display(plt)
    println("\nDone! Press Enter to exit...")
    readline()

else
    println("ERROR: File not found at path -> $img_path")    
end