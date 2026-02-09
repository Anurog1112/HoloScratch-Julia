import Pkg
Pkg.add(["Images", "FFTW", "ImageFiltering", "ImageSegmentation", "Plots", "FileIO"])


# define seting env
const wavelength = 532e-9
const dx_pixel_size = 3.45e-6
const distance = 0.1
const dy_pixel_size = dx_pixel_size

# show for check Configuration
println("=== System Configuration ===")
println("wavelength: $(wavelength * 1e9) nm")
println("distance: $(distance * 100) cm")
println("Pixel Size $(dx * 1e6) um")

# loda data
img_path = "data-test/145.png"

if isfile(img_path)
    # check path
    # println("path imag is $imag_path")
    img = load(img_path)
    img_gray = Gray.(img)
    u_hologram = Float64.(img_gray)

    # ดึงขนาดภาพจริงจากไฟล์
    (M, N) = size(u_hologram)

    # show imgae (Raw) to be confident
    p1 = heatmap(u_hologram, 
        title ="Raw Hologram ($M x $N)",
        xlabel = "X (pixels)",
        yield = "Y (pixels)",
        color = :gray,
        aspect_ratio = :equal,
        size = (600 , 600)
    )
    display(p1)
else
    println("ERROR: File not found! : $imag_path")    
end