using VortexLattice
using Plots
using Statistics

using ForwardDiff 

#///////////////////////////////////////////////////////////////////////
# 1. FUNCTION TO CREATE WING AND RUN VLM (Same as before)
#///////////////////////////////////////////////////////////////////////

function Create_wing_and_analyze(c_vec)
    c_root = c_vec[1]
    c_mid = c_vec[2]
    c_tip = c_vec[3]

    # --- Geometry Setup ---
    
    b_span = 10.0

    yle = [0.0, b_span/4, b_span/2] 
    zle = [0.0, 0.0, 0.0]
    chords = [c_root, c_mid, c_tip]
    xle = [-c_root/2, -c_mid/2, -c_tip/2] 

    ns_total = 20 
    nc = 1       
    
    theta = [2.0*pi/180, 2.0*pi/180, 2.0*pi/180] .* 0 #.* 0 wasn't there before
    phi = [0.0, 0.0, 0.0]
    fc = fill((xc)->0, 3) 
    
    spacing_s = Cosine() 
    spacing_c = Uniform()

    grid, ratio = wing_to_grid(xle, yle, zle, chords, theta, phi, ns_total, nc; 
        fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)  
    
    system = System([grid]; ratios=[ratio])

    
    Sref = 50.0
    cref = 2.0
    bref = 10.0
    Vinf = 1.0
    alpha = 1.0*pi/180
    rref = [0.50, 0.0, 0.0]
    symmetric = false


    ref = Reference(Sref, cref, bref, rref, Vinf)
    fs = Freestream(Vinf, alpha, 0.0, [0.0, 0.0, 0.0])
    

    steady_analysis!(system, ref, fs; symmetric=symmetric)
    

    cf, cm = lifting_line_coefficients(system; frame=Body())
    CDiff = far_field_drag(system) # Induced Drag (CD_i)
    CL = cf[1][3]                 # Lift Coefficient (Cz)
    
    return CDiff, CL, system #system is there to put the wing on paraview. Otherwise take it out
end 












#/////////////////////////////////////////////////////////////////////// FROM CHAT 2. and 3. Need to go and fully understand it or change it
# 2. FUNCTION TO PERFORM 1D SWEEP AND PLOT
#///////////////////////////////////////////////////////////////////////

function Sweep_One_Variable_And_Plot(variable_index, range_values, fixed_c_vec)
    
    # Initialize array to store the Induced Drag results
    CDi_values = zeros(length(range_values))
    
    # Perform the sweep
    for (i, val) in enumerate(range_values)
        # Create the current chord vector by copying the fixed vector
        c_vec = deepcopy(fixed_c_vec) #was copy
        
        # Set the swept variable (c1, c2, or c3) to the current value
        c_vec[variable_index] = val
        
        # Calculate the objective function (CDi)
        CDi, _, system = Create_wing_and_analyze(c_vec) #System is there to create wing on paraview
        CDi_values[i] = CDi
    end

    # Plotting based on the swept variable
    if variable_index == 1
        label = "Root Chord, \$c_1\$ (m)"
        fixed_str = "(Fixed: c2=$(fixed_c_vec[2])m, c3=$(fixed_c_vec[3])m)"
    elseif variable_index == 2
        label = "Mid Chord, \$c_2\$ (m)"
        fixed_str = "(Fixed: c1=$(fixed_c_vec[1])m, c3=$(fixed_c_vec[3])m)"
    else # variable_index == 3
        label = "Tip Chord, \$c_3\$ (m)"
        fixed_str = "(Fixed: c1=$(fixed_c_vec[1])m, c2=$(fixed_c_vec[2])m)"
    end
    
    p = Plots.plot(range_values, CDi_values, 
        xlabel=label, 
        ylabel="Induced Drag Coefficient (\$C_{D_i}\$)", 
        title="1D Sweep of \$C_{D_i}\$ vs. Chord Position \n $fixed_str",
        linewidth=3,
        legend=false)

    display(p)
    return p
end

#///////////////////////////////////////////////////////////////////////
# 3. RUN ALL THREE SWEEPS
#///////////////////////////////////////////////////////////////////////

# Define the ranges for the design variables
c1_range = range(0.5, stop=6.0, length=20) #was 4.0, stop=6.0
c2_range = range(0.5, stop=6.0, length=20) #was 3.0, stop=5.0
c3_range = range(0.5, stop=6.0, length=20) #was 0.5, stop=3.0

# Define a reasonable baseline chord vector (initial guess)
# The optimal elliptic chord shape has a ratio of C_root/C_tip between 3 and 6,
# with C_mid being somewhere in between.
BASELINE_CHORDS = [5.0, 3.5, 1.0] #Was 5.0, 5.0, 5.0

println("Starting 1D sweeps to verify objective function behavior...")

# Sweep 1: Vary Root Chord (c1), fix c2 and c3
p1 = Sweep_One_Variable_And_Plot(1, c1_range, BASELINE_CHORDS)

# Sweep 2: Vary Mid Chord (c2), fix c1 and c3
p2 = Sweep_One_Variable_And_Plot(2, c2_range, BASELINE_CHORDS)

# Sweep 3: Vary Tip Chord (c3), fix c1 and c2
p3 = Sweep_One_Variable_And_Plot(3, c3_range, BASELINE_CHORDS)

println("All 1D plots displayed.")




#This is to look at the wing to see if it looks right
a, b , system = Create_wing_and_analyze(BASELINE_CHORDS)
write_vtk("HWfourWing", system)


#///////////////////////////////////////////////////////////////////////
# 4. Jacobian
#///////////////////////////////////////////////////////////////////////


function sweep_then_Jacobian(range_values, fixed_c_vec)
        # Initialize array to store the Induced Drag results
    CDi_values = zeros(length(range_values))
        # Perform the sweep
    
    for i in range_values
        # Create the current chord vector by copying the fixed vector
        c_vec = deepcopy(fixed_c_vec) #was copy
        
        # Calculate the objective function (CDi)
        CDi, _, system = Create_wing_and_analyze(c_vec)
        push!(CDi_values, CDi)
    end


    jacobian = ForwardDiff(CDi_values)
    return jacobian 
end 

jacobian = sweep_then_Jacobian(c3_range, BASELINE_CHORDS) #Varies the tip chords
print(jacobian)

