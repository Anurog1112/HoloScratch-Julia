using Images, FileIO, Plots, FFTW, Statistics, ImageFiltering
gr()

const wavelength = 633e-9
const dx_pixel_size = 5.807e-6
const distance = 0.1
const dy_pixel_size = dx_pixel_size
const λ = wavelength
const z = distance

# ============================================
# Helper Function: Spherical & Tilt Compensation (2nd-Order Fit)
# ============================================
function remove_spherical_aberration(img)
    rows, cols = size(img)
    X = repeat(1:cols, 1, rows)'
    Y = repeat(1:rows, 1, cols)
    
    # 2nd-order Polynomial: z = c1*x^2 + c2*y^2 + c3*xy + c4*x + c5*y + c6
    X1 = reshape(X, :)
    Y1 = reshape(Y, :)
    A = [X1.^2  Y1.^2  X1.*Y1  X1  Y1  ones(length(X1))]
    b = reshape(img, :)
    
    coeffs = A \ b
    
    background = coeffs[1] .* X.^2 .+ coeffs[2] .* Y.^2 .+ coeffs[3] .* (X .* Y) .+ 
        coeffs[4] .* X .+ coeffs[5] .* Y .+ coeffs[6]
    return img .- background
end

# ============================================
# Helper Function: Auto-Focus (Sharpness Metric)
# ============================================
function auto_focus(U_field, H_template, fx_grid, fy_grid, z_range;
    dx=dx_pixel_size, wl=wavelength)
    k = 2π / wl
    A_0 = fftshift(fft(U_field))

    best_z = z_range[1]
    best_score = -Inf

    for z_try in z_range
        root = Complex.(1 .- (wl .* fx_grid).^2 .- (wl .* fy_grid).^2)
        H_try = exp.(1im .* k .* z_try .* sqrt.(root))
        H_try[real.(root) .< 0] .= 0

        U_try = ifft(ifftshift(A_0 .* H_try))
        amp = abs.(U_try)

        # Focus Metric: Normalised Variance (สูง = คมชัด)
        μ = mean(amp)
        score = var(amp) / (μ^2 + 1e-30)

        if score > best_score
            best_score = score
            best_z = z_try
        end
    end
    return best_z, best_score
end

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
        color = :viridis,
        aspect_ratio = :equal,
        axis = nothing
    )

    plt_1 = plot(p1, layout=(1, 1), size=(1000, 600))
    savefig(plt_1, "Raw Hologram.png")

    # display(plt_1)
    # println("\nDone! Press Enter to exit...")
    # readline()
    
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

    plt_2 = plot(p2, layout=(1, 1), size=(1000, 600))
    savefig(plt_2, "Frequency Spectrum (Log).png")

    # display(plt_2)
    # println("\nDone! Press Enter to exit...")
    # readline()

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

    # creat Gaussian Soft-Edge Mask (ขอบนุ่ม ลด Ringing Artifact)
    filter_redius = 20 # <--- ขนาดรูรับแสง (ถ้าภาพเบลอให้ปรับค่านี้)
    mask = zeros(ComplexF64, h, w)    
    
    for i in 1:h, j in 1:w
        dist = sqrt((i - peak_y)^2 + (j - peak_x)^2)

        # Gaussian window: ค่อยๆ จางลงที่ขอบ แทนการตัดแบบคมชัด
        mask[i, j] = exp(-dist^2 / (2 * filter_redius^2))
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
        title="Filtered Sideband",
        axis=nothing
    )

    p4 = heatmap(
        log.(abs.(F_centered) .+ 1),
        c=:viridis,
        title="Centered Sideband",
        axis=nothing
    )

    plt_3 = plot(p3, p4, layout=(1, 2), size=(1000, 600))
    savefig(plt_3, "Filtered sidedand.png")

    # display(plt_3)
    # println("\nDone! Press Enter to exit...")
    # readline() 
    
    println("\n=== step 3: ASM Propagation ===")
    # Creat Frequency Coordinates and Creat asxis Frequency fx, fy ให้ตรงกับภาพ
    (M, N) = size(U_filtered)
    fx = range(-1/(2*dx_pixel_size), 1/(2*dx_pixel_size), length=N)
    fy = range(-1/(2*dy_pixel_size), 1/(2*dy_pixel_size), length=M)

    # Creat Grid เพื่อให้คำนวณพิกเซลพร้อมกัน
    FX = repeat(reshape(fx, 1, :), M, 1)
    FY = repeat(reshape(fy, :, 1), 1, N)

    # ======== Image Plane Holography Mode ========
    println("\n--- Image Plane Reconstruction ---")
    # เลนส์ 35mm ทำหน้าที่สร้างภาพแล้ว ระยะแพร่คลื่น (z) จึงต้องใกล้ 0
    # ใช้ค่า z เล็กน้อยมากเพื่อชดเชย Defocus (0.0 cm ถึง 0.5 cm)
    z_best = 0.00  # <--- ตั้งเป็นศูนย์ หรือปรับแค่ 0.001 - 0.005 ถ้ากล้องโฟกัสคลาดเคลื่อน
    println("Using z = $(round(z_best*100, digits=2)) cm (Image Plane mode, 35mm lens)")

    # Creat Transfer Function (H) ด้วย z_best
    k = 2 * π / λ
    root_term = Complex.(1 .- (λ .* FX).^2 .- (λ .* FY).^2)

    # สร้างแผ่นกรองแสง H
    H = exp.(1im .* k .* z_best .* sqrt.(root_term))
    # ตัดคลื่น Evanescent (คลื่นที่เดินทางไม่ได้) ทิ้งไป เพื่อลด Noise
    H[real.(root_term) .< 0] .= 0

    # ย้ายแสงไปที่วัตถุ (Propagation)
    A_0 = fftshift(fft(U_filtered))
    A_z = A_0 .* H

    # แปลงกลับเป็นภาพ (Inverse FFT)
    U_recon = ifft(ifftshift(A_z))

    # Gaussian Smoothing เพื่อลด Phase Noise (σ=5.0)
    println("Applying Gaussian smoothing (σ=5.0)...")
    U_recon_smooth = imfilter(U_recon, Kernel.gaussian(5.0))

    # Amplitude จาก Raw (ไว้ดูโฟกัส) และ Phase จาก Smoothed (ไว้คำนวณความลึก)
    amp_img = abs.(U_recon)            # ความสว่างจาก Raw field
    phase_img = angle.(U_recon_smooth) # เฟสจาก Smoothed field

    # Log-scale + viridis เพื่อเพิ่ม Contrast (แก้ภาพสีเทาว่างเปล่า)
    amp_log = log.(amp_img .+ 1e-10)
    p5 = heatmap(amp_log,
        c=:viridis,
        title="Reconstructed Amplitude (z=$(round(z_best*100, digits=3)) cm)",
        aspect_ratio=:equal,
        axis=nothing
    )

    p6 = heatmap(phase_img,
        c=:viridis,
        title="Reconstructed Phase (Wrapped)", 
        aspect_ratio=:equal,
        axis=nothing
    )

    plt_4 = plot(p5, p6, layout=(1, 2), size=(1000, 600))
    savefig(plt_4, "Reconstructed Amplitude and  Phase.png")

    # display(plt_4)
    # println("\nDone! Press Enter to exit...")
    # readline() 

    println("\n=== Step 4 & 5: Unwrap & Scaling ===")
    function unwrap_1d(phase_vec)
        unwrapped = copy(phase_vec)
        diffs = diff(phase_vec)

        jumps = round.(diffs ./ (2*π))
        correction = cumsum([0; -jumps .* 2π])

        return unwrapped .+ correction
    end

    function unwrap_2d_simple(img_phase)
        (h, w) = size(img_phase)
        unwrapped_img = zeros(h, w)

        # Unwrap แถวกลาง (Middle Row) เพื่อเป็นแกนอ้างอิง
        mid_y = Int(h/2)
        unwrapped_img[mid_y, :] = unwrap_1d(img_phase[mid_y, :])

        # Unwrap แนวตั้ง (Columns) โดยอิงจากแถวกลางออกไปบน-ล่าง
        for x in 1:w
            col_up = reverse(img_phase[1:mid_y, x])
            unwrapped_col_up = unwrap_1d(col_up)

            offset_up = unwrapped_img[mid_y, x] - unwrapped_col_up[1]
            unwrapped_img[1:mid_y, x] = reverse(unwrapped_col_up .+ offset_up)

            col_down = img_phase[mid_y:end, x]
            unwrapped_col_down = unwrap_1d(col_down)
            
            offset_down = unwrapped_img[mid_y, x] - unwrapped_col_down[1]
            unwrapped_img[mid_y:end, x] = unwrapped_col_down .+ offset_down
        end

        return unwrapped_img
    end

    # start do Unwrap
    println("Unwrapping Phase (Simple 2D)...")
    phi_unwrapped = unwrap_2d_simple(phase_img)

    # Spherical Aberration Compensation: ลบความโค้งของเลนส์ (2nd-Order Fit)
    println("Applying Spherical Aberration Compensation...")
    phi_flat = remove_spherical_aberration(phi_unwrapped)

    # Scaling: แปลง Phase (rad) -> Height (nm)
    height_map_nm = (phi_flat .* wavelength) ./ (4 * π) .* 1e9

    # Median Filter : กำจัดจุดกระโดด (Spikes) บน Height Map
    println("Applying Median Filter (11×11) on height map...")
    height_map_nm = mapwindow(median, height_map_nm, (11, 11))

    println("\nFocused at z = $(round(z_best*100, digits=3)) cm")
    
    # ปรับระดับพื้นผิวให้เริ่มที่ 0 (ใช้ median เพื่อความเสถียร)
    height_map_nm .-= median(height_map_nm)

    # กลับทิศ (Optional): ถ้ารอยขีดข่วนมันนูนขึ้นแทนที่จะบุ๋มลง ให้แก้เครื่องหมาย
    # height_map_nm = -height_map_nm

    println("Max Height/Depth: $(maximum(height_map_nm)) nm")
    println("Min Height/Depth: $(minimum(height_map_nm)) nm")

    # Show 3D Surface (Output final)
    scale_factor = 4
    h_small = height_map_nm[1:scale_factor:end, 1:scale_factor:end]
    
    p_final = surface(h_small,
        title = "Reconstructed 3D Surface Defect",
        camera = (30, 45),
        shading = :flat,
        c = :viridis,
        zlabel = "Height (nm)",
        xlabel = "X", ylabel = "Y",
    )

    # # show All Plots
    plt = plot(p_final, layout=(1, 1), size=(800, 600))
    output_path = "Surface_final.png"
    savefig(plt, output_path)

    display(plt)
    println("\nDone! Press Enter to exit...")
    readline()
    println("=== PROCESS COMPLETE ===")

else
    println("ERROR: File not found at path -> $img_path")    
end