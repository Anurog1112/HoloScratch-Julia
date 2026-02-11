# ═══════════════════════════════════════════════════════════════════
#  HoloScratch GUI — Digital Holography Scratch Detection
#  Native Desktop Application using GLMakie (No Web)
#  Features: File Picker, Run Button, Dashboard, Interactive 3D,
#            Save All Images
# ═══════════════════════════════════════════════════════════════════

using GLMakie
using FileIO                          # load()
using Images: Gray                    # Gray color type only (avoid Axis conflict)
using FFTW, Statistics, ImageFiltering

GLMakie.activate!(title = "HoloScratch — Scratch Detection System")

# ───────────────────────────────────────────────
# Physical Constants
# ───────────────────────────────────────────────
const WL = 632e-9       # Wavelength (m) — He-Ne laser
const DX = 3.45e-6      # Pixel pitch (m)
const DY = DX

# ───────────────────────────────────────────────
# Helper Functions
# ───────────────────────────────────────────────
function remove_tilt(img)
    rows, cols = size(img)
    X = repeat(1:cols, 1, rows)'
    Y = repeat(1:rows, 1, cols)
    A_m = [reshape(X, :) reshape(Y, :) ones(rows * cols)]
    coeffs = A_m \ reshape(img, :)
    return img .- (coeffs[1] .* X .+ coeffs[2] .* Y .+ coeffs[3])
end

function auto_focus(U, FX, FY, zs)
    k  = 2π / WL
    A0 = fftshift(fft(U))
    bz, bs = zs[1], -Inf
    for zt in zs
        rt = Complex.(1 .- (WL .* FX).^2 .- (WL .* FY).^2)
        H  = exp.(1im * k * zt .* sqrt.(rt))
        H[real.(rt) .< 0] .= 0
        a  = abs.(ifft(ifftshift(A0 .* H)))
        μ  = mean(a)
        s  = var(a) / (μ^2 + 1e-30)
        s > bs && (bs = s; bz = zt)
    end
    return bz, bs
end

unwrap_1d(v) = v .+ cumsum([0; -round.(diff(v) ./ 2π) .* 2π])

function unwrap_2d(ph)
    h, w = size(ph)
    out  = zeros(h, w)
    mid  = h ÷ 2
    out[mid, :] = unwrap_1d(ph[mid, :])
    for x in 1:w
        u = unwrap_1d(reverse(ph[1:mid, x]))
        out[1:mid, x] = reverse(u .+ (out[mid, x] - u[1]))
        d = unwrap_1d(ph[mid:end, x])
        out[mid:end, x] = d .+ (out[mid, x] - d[1])
    end
    return out
end

# ───────────────────────────────────────────────
# Windows Native File Dialogs (PowerShell)
# ───────────────────────────────────────────────
function open_file_dialog()
    Sys.iswindows() || return ""
    scr = raw"""
Add-Type -AssemblyName System.Windows.Forms
$d = New-Object System.Windows.Forms.OpenFileDialog
$d.Filter = "Images (*.png;*.jpg;*.bmp;*.tif)|*.png;*.jpg;*.bmp;*.tif|All (*.*)|*.*"
$d.Title = "Select Hologram Image"
if ($d.ShowDialog() -eq 'OK') { Write-Output $d.FileName }
"""
    f = tempname() * ".ps1"; write(f, scr)
    r = try strip(read(`powershell -NoProfile -ExecutionPolicy Bypass -File $f`, String))
    catch; "" end
    rm(f, force=true); return r
end

function pick_save_folder()
    Sys.iswindows() || return ""
    scr = raw"""
Add-Type -AssemblyName System.Windows.Forms
$d = New-Object System.Windows.Forms.FolderBrowserDialog
$d.Description = "Save all result images to..."
if ($d.ShowDialog() -eq 'OK') { Write-Output $d.SelectedPath }
"""
    f = tempname() * ".ps1"; write(f, scr)
    r = try strip(read(`powershell -NoProfile -ExecutionPolicy Bypass -File $f`, String))
    catch; "" end
    rm(f, force=true); return r
end

# ───────────────────────────────────────────────
# Main Processing Pipeline
# ───────────────────────────────────────────────
function run_pipeline!(path, axes, lscene, obs, data)
    ax1, ax2, ax3, ax4, ax5, ax6 = axes
    st, iz, ih = obs

    # ── Step 1: Load ──
    st[] = "Loading image..."; yield()
    img = Float64.(Gray.(load(path)))
    M, N = size(img)

    empty!(ax1)
    heatmap!(ax1, img, colormap=:viridis)
    ax1.title = "1. Raw Hologram ($N×$M)"
    data["1_Raw_Hologram"] = img

    # ── Step 2: FFT ──
    st[] = "Step 2: FFT & Spectral Analysis..."; yield()
    Fs   = fftshift(fft(img))
    spec = log.(abs.(Fs) .+ 1)
    h, w = size(spec)

    empty!(ax2)
    heatmap!(ax2, spec, colormap=:viridis)
    scatter!(ax2, [w ÷ 2 + 1], [h ÷ 2 + 1], color=:red, markersize=8)
    ax2.title = "2. Frequency Spectrum"
    data["2_Frequency_Spectrum"] = spec

    # ── Spatial Filtering ──
    st[] = "Step 2: Spatial Filtering..."; yield()
    cx, cy = w ÷ 2, h ÷ 2
    ss = copy(spec)
    for i in 1:h, j in 1:w
        (i - cy)^2 + (j - cx)^2 < 2500 && (ss[i, j] = 0)   # mask_radius_dc = 50
    end
    _, mi = findmax(ss)
    py, px = mi[1], mi[2]

    # Gaussian soft-edge mask (r=35)
    r = 35
    mask = [exp(-((i - py)^2 + (j - px)^2) / (2 * r^2)) for i in 1:h, j in 1:w]
    Ff = Fs .* mask
    Fc = circshift(Ff, (cy - py, cx - px))
    Uf = ifft(ifftshift(Fc))

    fv = log.(abs.(Ff) .+ 1)
    cv = log.(abs.(Fc) .+ 1)
    empty!(ax3); heatmap!(ax3, fv, colormap=:viridis)
    ax3.title = "3. Filtered Sideband"
    empty!(ax4); heatmap!(ax4, cv, colormap=:viridis)
    ax4.title = "4. Centered"
    data["3_Filtered_Sideband"] = fv
    data["4_Centered_Sideband"] = cv

    # ── Auto-Focus ──
    st[] = "Step 3: Auto-Focus (coarse)..."; yield()
    M2, N2 = size(Uf)
    fx = range(-1 / (2DX), 1 / (2DX), length=N2)
    fy = range(-1 / (2DY), 1 / (2DY), length=M2)
    FX = repeat(reshape(fx, 1, :), M2, 1)
    FY = repeat(reshape(fy, :, 1), 1, N2)

    z1, _ = auto_focus(Uf, FX, FY, range(0.05, 0.15, step=0.005))
    st[] = "Step 3: Auto-Focus (fine)..."; yield()
    zb, _ = auto_focus(Uf, FX, FY, range(max(0.01, z1 - 0.01), z1 + 0.01, step=0.0005))

    # ── ASM Propagation ──
    st[] = "Step 3: ASM Propagation (z=$(round(zb*100, digits=2)) cm)..."; yield()
    k  = 2π / WL
    rt = Complex.(1 .- (WL .* FX).^2 .- (WL .* FY).^2)
    H  = exp.(1im * k * zb .* sqrt.(rt))
    H[real.(rt) .< 0] .= 0
    Ur = ifft(ifftshift(fftshift(fft(Uf)) .* H))

    # Gaussian smoothing σ=5.0 on complex field
    Us  = imfilter(Ur, Kernel.gaussian(5.0))
    amp = abs.(Ur)           # Amplitude from raw (for focus check)
    ph  = angle.(Us)         # Phase from smoothed (for depth)

    al = log.(amp .+ 1e-10)
    empty!(ax5); heatmap!(ax5, al, colormap=:viridis)
    ax5.title = "5. Amplitude (z=$(round(zb*100, digits=1)) cm)"
    empty!(ax6); heatmap!(ax6, ph, colormap=:viridis)
    ax6.title = "6. Phase (Wrapped)"
    data["5_Amplitude"] = al
    data["6_Phase_Wrapped"] = ph

    # ── Phase Unwrapping + Height ──
    st[] = "Step 4: Phase Unwrapping..."; yield()
    phu = unwrap_2d(ph)
    phf = remove_tilt(phu)

    st[] = "Step 5: Height Map + Smoothing..."; yield()
    hm  = (phf .* WL) ./ (4π) .* 1e9
    hm  = mapwindow(median, hm, (11, 11))
    hm  = imfilter(hm, Kernel.gaussian(3.0))
    hm .-= median(hm)

    mx = round(maximum(hm), digits=1)
    mn = round(minimum(hm), digits=1)
    iz[] = "Focus: z = $(round(zb * 100, digits=2)) cm"
    ih[] = "Height: $mn ~ $mx nm"

    # ── 3D Surface ──
    st[] = "Rendering 3D Surface..."; yield()
    sf = 4
    hs = Float64.(hm[1:sf:end, 1:sf:end])

    # Clear and re-draw 3D
    while length(lscene.scene.plots) > 0
        delete!(lscene.scene, lscene.scene.plots[end])
    end
    surface!(lscene, 1:size(hs, 1), 1:size(hs, 2), hs, colormap=:viridis)

    data["7_3D_Surface"] = hm
    st[] = "[DONE] z=$(round(zb*100, digits=2)) cm | Height: $mn ~ $mx nm | Drag 3D to rotate, scroll to zoom"
    return nothing
end

# ───────────────────────────────────────────────
# Build GUI Layout
# ───────────────────────────────────────────────
function create_gui()
    fig = Figure(size=(1500, 1000), fontsize=13)

    # ── Title ──
    Label(fig[0, :], "HoloScratch - Digital Holography Scratch Detection",
          fontsize=22, font=:bold, color=:gray20)

    # ── Control Bar ──
    cg = fig[1, 1:4] = GridLayout(tellwidth=false)
    path_obs = Observable("No image selected")
    Label(cg[1, 1], "Image:", fontsize=14, font=:bold)
    Label(cg[1, 2], path_obs, fontsize=12, color=:gray40,
          tellwidth=false, halign=:left)
    btn_browse = Button(cg[1, 3], label=" Browse ",  fontsize=13)
    btn_run    = Button(cg[1, 4], label=" > Run ",   fontsize=13,
                        buttoncolor=RGBf(0.65, 0.92, 0.65))
    btn_save   = Button(cg[1, 5], label=" Save All ", fontsize=13)

    # ── 2D Plots (Row 2: 4 plots, Row 3: 2 plots + info) ──
    ax1 = Axis(fig[2, 1], title="1. Raw Hologram",      aspect=DataAspect())
    ax2 = Axis(fig[2, 2], title="2. Frequency Spectrum", aspect=DataAspect())
    ax3 = Axis(fig[2, 3], title="3. Filtered Sideband",  aspect=DataAspect())
    ax4 = Axis(fig[2, 4], title="4. Centered Sideband",  aspect=DataAspect())
    ax5 = Axis(fig[3, 1], title="5. Amplitude",          aspect=DataAspect())
    ax6 = Axis(fig[3, 2], title="6. Phase (Wrapped)",    aspect=DataAspect())
    for a in [ax1, ax2, ax3, ax4, ax5, ax6]
        hidedecorations!(a)
    end

    # ── Info Panel ──
    ig = fig[3, 3:4] = GridLayout()
    Label(ig[1, 1], "[Parameters]", fontsize=16, font=:bold, halign=:center)
    z_obs = Observable("Focus: z = — cm")
    h_obs = Observable("Height: — nm")
    Label(ig[2, 1], z_obs, fontsize=14, halign=:center, color=:gray20)
    Label(ig[3, 1], h_obs, fontsize=14, halign=:center, color=:gray20)
    Label(ig[4, 1], "WL = $(WL*1e9) nm   dx = $(DX*1e6) um",
          fontsize=11, color=:gray50, halign=:center)
    Label(ig[5, 1], "Mask: Gauss(r=35)  Smooth: s=5 > Med(11) > s=3",
          fontsize=10, color=:gray55, halign=:center)

    # ── 3D Surface (Interactive: Drag=Rotate, Scroll=Zoom) ──
    lscene = LScene(fig[4, 1:4], show_axis=true)

    # ── Status Bar ──
    st_obs = Observable("Ready -- Click Browse to select a hologram image")
    Label(fig[5, 1:4], st_obs, fontsize=13, color=:dodgerblue,
          halign=:left, tellwidth=false)

    # ── Internal State ──
    axes      = (ax1, ax2, ax3, ax4, ax5, ax6)
    obs       = (st_obs, z_obs, h_obs)
    sel_path  = Ref("")
    plot_data = Dict{String, Any}()
    is_busy   = Ref(false)

    # ────────── Button: Browse ──────────
    on(btn_browse.clicks) do _
        p = open_file_dialog()
        if !isempty(p) && isfile(p)
            sel_path[] = p
            path_obs[] = basename(p)
            st_obs[] = "[OK] $(basename(p)) selected -- Click Run to process"
        end
    end

    # ────────── Button: Run ──────────
    on(btn_run.clicks) do _
        is_busy[] && (st_obs[] = "[!] Already running!"; return)
        isempty(sel_path[]) && (st_obs[] = "[!] Select an image first!"; return)
        !isfile(sel_path[]) && (st_obs[] = "[!] File not found!"; return)
        is_busy[] = true
        @async begin
            try
                empty!(plot_data)
                run_pipeline!(sel_path[], axes, lscene, obs, plot_data)
            catch e
                st_obs[] = "[ERR] Error: $(sprint(showerror, e))"
                @error "Pipeline error" exception=(e, catch_backtrace())
            finally
                is_busy[] = false
            end
        end
    end

    # ────────── Button: Save All ──────────
    on(btn_save.clicks) do _
        isempty(plot_data) && (st_obs[] = "[!] Run processing first!"; return)
        folder = pick_save_folder()
        isempty(folder) && return
        @async begin
            try
                st_obs[] = "Saving images..."; yield()

                # Save full dashboard
                Makie.save(joinpath(folder, "00_dashboard.png"), fig, px_per_unit=2)
                n = 1

                # Save individual plots
                for (name, d) in sort(collect(plot_data))
                    tf = Figure(size=(800, 600))
                    if name == "7_3D_Surface"
                        s = 4; ds = d[1:s:end, 1:s:end]
                        ls = LScene(tf[1, 1])
                        surface!(ls, 1:size(ds, 1), 1:size(ds, 2),
                                 Float64.(ds), colormap=:viridis)
                    else
                        ax = Axis(tf[1, 1],
                                  title=replace(name, "_" => " "),
                                  aspect=DataAspect())
                        hm = heatmap!(ax, d, colormap=:viridis)
                        Colorbar(tf[1, 2], hm)
                    end
                    Makie.save(joinpath(folder, "$name.png"), tf)
                    n += 1
                end
                st_obs[] = "[OK] Saved $n images to: $(basename(folder))/"
            catch e
                st_obs[] = "[ERR] Save error: $(sprint(showerror, e))"
            end
        end
    end

    display(fig)
    return fig
end

# ═══════════════════════════════════════════════
# Launch Application
# ═══════════════════════════════════════════════
println("╔═══════════════════════════════════════════════════════╗")
println("║  Loading HoloScratch GUI...                          ║")
println("║  (First launch compiles GLMakie — may take 1-3 min)  ║")
println("╚═══════════════════════════════════════════════════════╝")

fig = create_gui()

println()
println("  ✅ HoloScratch GUI is running!")
println()
println("  ┌──────────────────────────────────────┐")
println("  │  📁 Browse  → Select hologram image  │")
println("  │  ▶  Run     → Process hologram       │")
println("  │  🖱  3D View → Drag=Rotate Scroll=Zoom│")
println("  │  💾 Save    → Export all images       │")
println("  └──────────────────────────────────────┘")
println()
println("  Press Enter in this terminal to exit.")
readline()
