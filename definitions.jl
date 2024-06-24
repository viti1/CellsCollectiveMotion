## ~~~~~~~~~~~ Definitions file - Defining and Recieving all simulation parameters ~~~~~~~~~~~~~~~ ##

# In order to run the simulation : in command line go to the foler where all .jl files are stored
# >> <full_path_to_julia_compiler>\julia  main.jl <your_simulation_name> param1=value1 param2=value2 .....
# If using linux , it's of course useful to have an alias to <full_path_to_julia_compiler>\julia

# Possible Parameters list : 
# start_from_last 
# seed - random seed 
# n - number of iterations 
# print_step - number of iteration to pring into output files
# print_force - pring to output files also the forces on each cell
# two_sided - (true/false)
# strips - (true/false) 
# grad - gradual strips (true/false)
# periodic_x - default is false, true if parameter exist 
# box - default is false, true if parameter exist 
# u1 - interaction pulling force strength 
# a1 - interaction pulling force distance
# u3 - interaction pulling force strength
# a3 - interaction pulling force distance
# ufact - factor to multiply u1 and u3 from initial values 
# beta  - viscek magnitude coefficient
# rbeta - relative beta to beta 60 
# rel_noise - change the noise magnitude to achieve standart velocity distribution according to a given beta
# rnoise - relative value of sigma0 to 150 ( sigma1 will be twice sigma0  )
# sigma  - absolute value of sigma0. sigma1 will be twice sigma0  
# const_noise - constant noise (with no dependancy on density )
# hmax  - maximum curvature value in Fborder graph. 
# slope - slope/1000 of Fborder between curvature  0 and -hmax (for 25000 type 25)
# rslope  - relative slope value to 25,000
# fmin -  Fborder minimum level
# rfmin - Fborder relative minimum level (F_MIN = 75*rfmin)
# fmax -  Fborder Fmin
# kappa - Resoring force coefficient 
# cords - Actin cables magnitude
# max_to_sigma - convert 
# no_division - no prolifiration 
# min_div - T0 parameter in division rate formula
# thermalization_time - time before stencil removal (only right side for one-sided simulation, and both sides of two-sided simulation)
# dx - average distance between cells in x axis in initial cells position
# dy - average distance between cells in x axis in initial cells position
# x_len - length in microns of y-axis in initial cells position
# y_len - length in microns of y-axis in initial cells position
# near_nbr - the radius in wich neighbor will be searched
# gwidth - smoothing gaussian width 
include("beta2sigma.jl")

#-------- run directory  ---------
PRINT_TO_STDOUT = true
if length(ARGS)==2 && ARGS[2]=="59"
  ARGS = String[];
  running_from_editor = true;
else
  running_from_editor = false;
end
if isempty(ARGS) #default
    directory = "../..//sims/run1/"
else
    s = ARGS[1];
    if ismatch(r"=",s) ; error("Illegal simulation name '$s'"); end
    if ismatch(r"scratch",s)
        directory = "$s/";
        PRINT_TO_STDOUT = false;
    else
        directory = "../../sims/$s/"
    end
end

if  !isempty(ARGS) && ~ismatch(r"scratch",s) && ~isdir(directory) 
	mkdir(directory)
end



println(" <<<<<  dir = $directory >>>>>")

# ----------------------------
    if ismatch(r"scratch",directory)
        PRINT_EVERY_I = false;
    else
        PRINT_EVERY_I = true;
    end


    NUM_OF_CYCLES = 2
    if inargs("n") ; NUM_OF_CYCLES=extract_param("n") ; end

    start_from_last = false
    PRINT_STEP = 100
    if NUM_OF_CYCLES<=5 PRINT_STEP=1 end

    if inargs("start_from_last") ; start_from_last=true ; PRINT_FORCE=true; end
    if inargs("print_step") ; PRINT_STEP=extract_param("print_step") ; end


#-------- Strips parameters -----------------
STRIPS = true;
if inargs("strips"); STRIPS=extract_bool_param("strips"); end

if STRIPS
    STRIP_DELIMITER_WIDTH = 120;
    STRIP_WIDTH = [30];
    if inargs("strip_width"); STRIP_WIDTH=extract_param("strip_width"); end
    if inargs("delimiter"); STRIP_DELIMITER_WIDTH=extract_param("delimiter"); end

    STRIPS_BORDERS = zeros(length(STRIP_WIDTH),2)
    STRIPS_BORDERS[1,1] = 0.5*STRIP_DELIMITER_WIDTH;
    STRIPS_BORDERS[1,2] = STRIPS_BORDERS[1,1] + STRIP_WIDTH[1] ;
    for j = 2:length(STRIP_WIDTH)
      STRIPS_BORDERS[j,1] = STRIPS_BORDERS[j-1,2] + STRIP_DELIMITER_WIDTH;
      STRIPS_BORDERS[j,2] = STRIPS_BORDERS[j,1] + STRIP_WIDTH[j] ;
    end


	CHANNEL_LENGTH = Inf;
	if inargs("channel_len")
	  CHANNEL_LENGTH = extract_param("channel_len");
	end

	CORNER_RADIUS = 5.0;
	if inargs("corner_radius")
	  CORNER_RADIUS = extract_param("corner_radius");
	end

	if STRIP_DELIMITER_WIDTH / 2 < CORNER_RADIUS 
		error("STRIP_DELIMITER_WIDTH/2 must be larger than CORNER_RADIUS ! ")
	end
		
	WALL_FORCE = 100.0;
	if inargs("wall_force")
		WALL_FORCE = extract_param("wall_force");
	end

	WALL_FORCE_CUTOFF = 12.0; # micron
	if inargs("wall_force_cutoff")
		WALL_FORCE_CUTOFF = extract_param("wall_force_cutoff");
	end
end
# ------ gradual Strip Parameters -------------
GRADUAL_STRIP = false;
if inargs("grad")
  GRADUAL_STRIP = extract_bool_param("grad");
end

if GRADUAL_STRIP
  STRIPS = false;
end

if GRADUAL_STRIP
  #STRIP_WIDTH = 150; # base width
  #STRIP_SLOPE = -200/900; #dy/dx  +is converging
  #STRIP_DELIMITER_WIDTH = 700

  STRIP_WIDTH = 600; # base width
  STRIP_SLOPE = 200/900; #dy/dx  +is converging
  STRIP_DELIMITER_WIDTH =200
  if inargs("strip_base"); STRIP_WIDTH=extract_param("strip_base"); end
  if inargs("strip_slope_sign");  STRIP_SLOPE=extract_param("strip_slope_sign")*abs(STRIP_SLOPE); end
  if inargs("strip_slope");  STRIP_SLOPE=extract_param("strip_slope"); end
  if inargs("delimiter");  STRIP_DELIMITER_WIDTH=extract_param("delimiter"); end


  STRIPS_BORDERS = zeros(length(STRIP_WIDTH),2)
  STRIPS_BORDERS[1,1] = STRIP_DELIMITER_WIDTH/2;
  STRIPS_BORDERS[1,2] = STRIP_DELIMITER_WIDTH/2 + STRIP_WIDTH;

end

#--------   general parameters --------------
    TWO_SIDED = false
    NO_BORDERS = false
    if inargs("two_sided") TWO_SIDED=extract_bool_param("two_sided") ;end
    

    RANDOM_SEED =  198;
    if inargs("seed") ; RANDOM_SEED=int(extract_param("seed")) ; end

    DELETE_DEATACHED = false;
    MAX_VALUE_V = 10000
    SMOOTH_GAUSSIAN_WIDTH = 25
    if inargs("gwidth"); SMOOTH_GAUSSIAN_WIDTH = extract_param("gwidth"); end

#-------- interaction force -----------------
    BETA = 60
    if inargs("beta") ; BETA=extract_param("beta"); end
    if inargs("rbeta") ; BETA=10*extract_param("rbeta"); println("BETA=$BETA"); end

    TAU  = 1.39
    ALPHA  = 1.42  #h^-1

    U0  = 2650   # um^2/h
    A0  = 8       # um
    U2  = 30
    A2  = 2

    U1  = 2 #2.8   #1/h^2
    if inargs("u1"); U1=extract_param("u1"); end
    A1  = 25    # um
    if inargs("a1"); A1=extract_param("a1"); end
    U3  = 0.7
    if inargs("u3"); U3=extract_param("u3"); end
    A3  = A1 # original 26
    if inargs("a3") A3=extract_param("a3"); end

    if inargs("ufact")
       U_FACTOR = extract_param("ufact");
       U1 = U_FACTOR*U1;
       U3 = U_FACTOR*U3;
    end
#-------- Noise force -----------------------
    SIGMA0 = 150.0
    SIGMA1 = 300.0

    # *** Three Options: **
    # 1) rel_noise  -  (adjust noise to keep P(v) - relevant for beta!=60 )
    # 2) rnoise=[desired relative value of sigma]
    # 3) sigma =[desired absolute value of sigma]

    if inargs("rel_noise") # Automatic detection of coefficient of noise to beta
        beta_to_sigma_coeff = beta2sigma(BETA);
        SIGMA0 = SIGMA0 * beta_to_sigma_coeff;
        SIGMA1 = SIGMA1 * beta_to_sigma_coeff;
        if beta_to_sigma_coeff!=1.0 println("beta_to_sigma_coeff =",beta_to_sigma_coeff); end
        if inargs("max_to_sigma") error("cat have both 'relative_noise' and 'max_to_sigma'"); end
    end

    if inargs("rnoise") # User defined coefficient of noise to beta, must have an "=" argument
        rel_noise = extract_param("rnoise");

        if inargs("max_to_sigma") error("cat have both 'relative_noise' and 'max_to_sigma'"); end

        if inargs("dsigma")
            DELTA_SIGMA = extract_param("dsigma")
            SIGMA0_0 = 0;
            if inargs("sigma00")
                SIGMA0_0 = extract_param("sigma00")
            end
            SIGMA0=SIGMA0_0 + rel_noise*DELTA_SIGMA;
            SIGMA1 = SIGMA0*2;
        else
            SIGMA0 = SIGMA0 * rel_noise;
            SIGMA1 = SIGMA1 * rel_noise;
        end

        println("SIGMA0=$SIGMA0")
    end

    if inargs("sigma")
        SIGMA0 = extract_param("sigma");
        SIGMA1 = SIGMA0*2;
    end

    if inargs("const_noise"); SIGMA1=SIGMA0; println("Constant noise $SIGMA1"); end

#-------- Contour forces --------------------
    F_H_MAX = 0.05 # um-1
    F_ALPHA = 25 * 1000 ; #m^2/h^2
    if inargs("hmax"); F_H_MAX=extract_param("hmax"); end
    if inargs("slope")  ;  F_ALPHA = extract_param("slope") * 1000          ; end
    if inargs("rslope");  F_ALPHA = 10 * 1000 * extract_param("rslope")   ; end
    F_CELL_MIN = 0;
    # two optionx for F_CELL_MIN:
    # fmin = [value]
    # rfmin = [relative to] - not valid right now!
    # if inargs("rfmin") ; F_CELL_MIN=F_CELL_MAX*extract_param("rfmin") ; end
    if inargs("fmin") ; F_CELL_MIN = extract_param("fmin") ; end
    if inargs("rfmin")
            DELTA_FMIN = 75;
            if inargs("dfmin")
                  DELTA_FMIN = extract_param("dfmin")
            end

            F_CELL_MIN=DELTA_FMIN*extract_param("rfmin") ;
    end
    F_CELL_MAX = F_CELL_MIN + F_ALPHA*F_H_MAX   # um/h^2

    if inargs("fmax")
        F_CELL_MAX = extract_param("fmax");
        F_ALPHA = (F_CELL_MAX-F_CELL_MIN)/F_H_MAX;
    end

    F_BETTA = 0.1 ;
    F_LAMDA_MAX = 250
    F_Q_MAX = 2*pi/F_LAMDA_MAX
    F_KAPPA = 1.5e7 #max2kappa(F_ALPHA) #F_ALPHA*( 1 - F_BETTA ) / (2*F_Q_MAX^2) # -> kappa/eta #2.5e7
    if inargs("kappa") ; F_KAPPA=extract_param("kappa") ; end

    F_CORD_MAG = 350
    if inargs("cords") ; F_CORD_MAG=extract_param("cords") ; end

    if inargs("max_to_sigma")
#            if inargs("sigma")
#                    error("Cant have both sigma and max_to_sigma");
#            end

#            if  inargs("rel_noise") ||  inargs("rnoise")
#                    # automatic or by hand sigma proportionality to beta
#                    println("Warning! Sigma defined by both max and beta!")
#            end

            if inargs("base")
                BASE_POINT = extract_param("base")
            else
                BASE_POINT = 25;
            end
            max_to_sigma_coeff = F_ALPHA / (BASE_POINT * 1000 ) #max2sigma(F_ALPHA)
            SIGMA0 = SIGMA0 * max_to_sigma_coeff;
            SIGMA1 = SIGMA1 * max_to_sigma_coeff;
            if inargs("max"); println("max_to_sigma_coeff =",max_to_sigma_coeff); end
    end



#-------  Cell division -------------------
    DIVIDE  = true
    if inargs("no_division") DIVIDE=false; end

    RHO0 = 4e-3
    MIN_DIV =  35 #T0
    RHO1 = 1/170 # =0.0059
    if inargs("rho1"); RHO1= extract_param("rho1"); end
    MIN_CURV_FOR_DIVIDE = 0.03 #originally -0.01
    if inargs("min_div"); MIN_DIV = extract_param("min_div"); end
    NEW_CELL_DIST = 5

#------ Inside a Box ---------------
    BOX=false;
    if inargs("box")
            BOX=true;
            println("input parameter: box option")
    end


# --------- basic parameters -------
    DT = 0.001
    THERMALIZATION_TIME = 0.0
    if inargs("thermalization_time"); # in hours
        THERMALIZATION_TIME=extract_param("thermalization_time")/DT;
    end
    dX = 16.0
    dY = 16.0
    if inargs("dx"); dX = extract_param("dx"); end
    if inargs("dy"); dY = extract_param("dy"); end

    X_LEN = 30
    if STRIPS || GRADUAL_STRIP
      Y_LEN = int(floor((STRIP_DELIMITER_WIDTH*(length(STRIP_WIDTH)) + sum(STRIP_WIDTH) ) / dY ))
    else
      Y_LEN = 70
    end
    if inargs("x_len") ; X_LEN=int(round(extract_param("x_len"))/dX) ; end
    if inargs("y_len") ; Y_LEN=int(round(extract_param("y_len"))/dY) ; end

  #  if !STRIPS && inargs("x_width") ; X_LEN=int(extract_param("x_width")) ; end
  #  if !STRIPS && inargs("y_width") ; Y_LEN=int(extract_param("y_width")) ; end

    NEAR_NBR = 70.0
    if inargs("near_nbr"); NEAR_NBR=extract_param("near_nbr"); end
    CRITICAL_ANGLE = 118 # in degrees

	if STRIPS && CHANNEL_LENGTH < NEAR_NBR
		error("channel_len must be larger than near_neighbor distance (".NEAR_NBR , ")")
	end

#  --- initial state ------

  	CIRCULAR_INITIAL_LOCATION = false
    if inargs("circular"); CIRCULAR_INITIAL_LOCATION = true; println("circular initial state"); end
    CR = 3*dX;
    if inargs("cr");  CR = extract_param("cr"); end

    CX = X_LEN*dX + CR + 40 ;
    CY = 70;

    if inargs("cx"); CX = extract_param("cx"); end
    if inargs("cy"); CY = extract_param("cy"); end
#----- the rest ---------------------
WEIGHTED_SMOOTH = true
CORDS = true
if ~isdefined(:PRINT_FORCE); PRINT_FORCE = false; end
if inargs("print_force") PRINT_FORCE=extract_bool_param("print_force"); end
CORDS_RIP = false # taken from here, not the loaded parametersSTRI
RIP_P = [1439,1456];
RIP = [610,725,430,495]; #[XMIN XMAX YMIN YMAX]
CORDS_STYLE = 2;
PERIODIC_Y = true
PERIODIC_X = false
if inargs("periodic_x") PERIODIC_X=true; println("periodic_x");
end
if PERIODIC_X; BOX=true; end #just to prevent contour forces


 # ----------- print / load parameters -----------------
if start_from_last
#        if inargs("origin_dir")
#           origin_dir = extract_string_param("origin_dir");
#           ;run(`cp "$(origin_dir)/last_pos.txt.txt `$(directory))
#           ;run(`cp "$(origin_dir)/parameters.txt $(directory)`)
#        end
        ;run(` awk -v print_step="$PRINT_STEP" -v n="$NUM_OF_CYCLES" -f awk_parameters "$(directory)parameters.txt"`  |> "$(directory)temp_file" )

        ;run(` mv $(directory)temp_file $(directory)parameters.txt`)
        #;run(`cat $(directory)parameters.txt`)
        ;run(` awk '/THERMALIZATION_TIME/ {print "THERMALIZATION_TIME=0;"; next}; {print}' "$(directory)parameters.txt"` |> "$(directory)temp_file")
        ;run(` mv $(directory)temp_file $(directory)parameters.txt`)
        include("$(directory)parameters.txt")
else
        param_file = open("$(directory)parameters.txt", "w");

		println(param_file,"% code directory: ",pwd());
        println(param_file,"NUM_OF_CYCLES =     ", NUM_OF_CYCLES   ,";" );
        println(param_file,"TWO_SIDED =     ", TWO_SIDED   ,";" );
        println(param_file,"CORDS =         ", CORDS       ,";" );
        println(param_file,"RANDOM_SEED =   ", RANDOM_SEED ,";" );
        println(param_file,"PRINT_STEP =    ", PRINT_STEP  ,";" );
        println(param_file,"DELETE_DEATACHED = ", DELETE_DEATACHED,";" );

        println(param_file,"\nBETA = "  ,BETA  ,";" );
        println(param_file,"TAU  = "    ,TAU   ,";" );
        println(param_file,"ALPHA  = "  ,ALPHA ,";" );
        println(param_file,"U0 = "      ,U0    ,";" );
        println(param_file,"U1 = "      ,U1    ,";" );
        println(param_file,"U2 = "      ,U2    ,";" );
        println(param_file,"U3 = "      ,U3    ,";" );
        println(param_file,"A0 = "      ,A0    ,";" );
        println(param_file,"A1 = "      ,A1    ,";" );
        println(param_file,"A2 = "      ,A2    ,";" );
        println(param_file,"A3 = "      ,A3    ,";" );

        println(param_file,"\nSIGMA0 = "  ,SIGMA0,";" );
        println(param_file,"SIGMA1 = "  ,SIGMA1,";" );

        println(param_file,"\nDIVIDE  =       ", DIVIDE      ,";" );
        println(param_file,"MIN_DIV = ",MIN_DIV    ,";" );
        println(param_file,"RHO0    = ",RHO0  ,";" );
        println(param_file,"RHO1 = ",RHO1    ,";" );

        println(param_file,"\nDT       = ",DT      ,";" );find_sigma_for_beta_45
        println(param_file,"THERMALIZATION_TIME = ",THERMALIZATION_TIME    ,";" );
        println(param_file,"dX  = ",dX      ,";" );DT
        println(param_file,"dY  = ",dY      ,";" );
        println(param_file,"X_LEN    = ",X_LEN   ,";" );
        println(param_file,"Y_LEN    = ", Y_LEN  ,";" );
        println(param_file,"BOARD_BORDER_X = ",dX*X_LEN ,";" );
        println(param_file,"BOARD_BORDER_Y = ",dY*Y_LEN ,";" );
        println(param_file,"PERIODIC_X = ", PERIODIC_X ,";" );
        println(param_file,"PERIODIC_Y = ", PERIODIC_Y ,";" );
        println(param_file,"BOX = ", PERIODIC_X ,";" );
        println(param_file,"NEAR_NBR = ",NEAR_NBR,";" );
        println(param_file,"CRITICAL_ANGLE = ",CRITICAL_ANGLE  ,";" );

        println(param_file,"\nF_H_MAX =    ",F_H_MAX   ,";")
        println(param_file,"F_CELL_MAX = ",F_CELL_MAX,";")
        println(param_file,"F_CELL_MIN = ",F_CELL_MIN,";")
        println(param_file,"F_ALPHA =    ",F_ALPHA   ,";")
        println(param_file,"F_BETTA =    ",F_BETTA   ,";")
        println(param_file,"F_LAMDA_MAX = ",F_LAMDA_MAX,";")
        println(param_file,"F_Q_MAX =    ",F_Q_MAX   ,";")
        println(param_file,"F_KAPPA =    ",F_KAPPA   ,";")
        println(param_file,"RD_MAG  = ",F_CORD_MAG,";")
        println(param_file,"WEIGHTED_SMOOTH = "  ,WEIGHTED_SMOOTH,";" );
        println(param_file,"SMOOTH_GAUSSIAN_WIDTH = ", SMOOTH_GAUSSIAN_WIDTH ,";" );

        # old stuff
        #println(param_file,"CORDS_RIP =     ", CORDS_RIP   ,";" );
        #println(param_file,"RIP = [ $(RIP[1]) $(RIP[2]) $(RIP[3]) $(RIP[4]) ] ;" );
        #println(param_file,"MAX_VALUE_V = ", MAX_VALUE_V,";" );
        #println(param_file,"CORDS_STYLE = ", CORDS_STYLE,";" );
        #println(param_file,"START_TIME = 0;");

        # columns stuff
        println(param_file,"STRIPS = ",STRIPS,";")
        if STRIPS || GRADUAL_STRIP
          println(param_file,"STRIP_DELIMITER_WIDTH = ",STRIP_DELIMITER_WIDTH,";") ;
          println(param_file,"STRIP_WIDTH = ",STRIP_WIDTH,";") ;
          print(param_file,"STRIPS_BORDERS = [ ")
          for stripid = 1:size(STRIPS_BORDERS,1)-1
            print(param_file,STRIPS_BORDERS[stripid,:],";") ;
          end
          println(param_file,STRIPS_BORDERS[end,:],"];");
          println(param_file,"CHANNEL_LENGTH = ",CHANNEL_LENGTH,";")
		  
		  println(param_file,"CORNER_RADIUS = ",CORNER_RADIUS,";")
		  println(param_file,"WALL_FORCE = ",WALL_FORCE,";")
		  println(param_file,"WALL_FORCE_CUTOFF = ",WALL_FORCE_CUTOFF,";")
        end

        println(param_file,"GRADUAL_STRIP = ",GRADUAL_STRIP,";")
        if GRADUAL_STRIP
            println(param_file,"STRIP_SLOPE = ",STRIP_SLOPE,";") ;
        end
        println(param_file,";")
        close(param_file);
end

#------------ Other sim constants --4560----------------------------
#must be afater load parameters

srand(RANDOM_SEED);

half_pi = pi/2
BOARD_BORDER_X = dX*X_LEN;
BOARD_BORDER_Y = dY*Y_LEN;
HALF_BOARD_BORDER_X = BOARD_BORDER_X / 2 ;
HALF_BOARD_BORDER_Y = BOARD_BORDER_Y / 2 ;
SQR_NEAR_NBR = (NEAR_NBR*NEAR_NBR)
COS_CRITICAL_ANGLE = cosd(CRITICAL_ANGLE);
CriticalAngle_Rad = CRITICAL_ANGLE * pi / 250;
Noise_Coeff = 1/sqrt(DT) ;
dt_divtau = DT/TAU;

# ---  Old parameters ---
TWO_P = false;
P1 = 1; #P2 = 2; P3 = 3; P4 = 4;
print();

#------------------- definition of type -------------------------
println("0")

if ~isdefined(:DEF_INCLUDED)
   type particle_t
                x::Float64
                y::Float64
                v::Array{Float64,1}
                eta::Array{Float64,1}
                nbrs::Array{Int,1}
                angles::Array{Float64,1}
                rho::Number
   end
   type corner_t
	    x::Float64
        y::Float64
		xsign::Float64
		ysign::Float64	
		r::Float64
		d::Float64 # wall cutoff
	end  

	type wall_t		
		p1::Array{Float64,1} 
		p2::Array{Float64,1} 
		d::Float64;  # wall cutoff
	end

end
DEF_INCLUDED = true;

# -----    help constant - corners -

if STRIPS
	CORNERS = corner_t[];  # make improvement - rounded corners at the entrance to the strip

	for j=1:length(STRIP_WIDTH)*4
	temp = corner_t(0.0,0.0,0.0,0.0,CORNER_RADIUS,WALL_FORCE_CUTOFF);
		push!(CORNERS,temp);
	end

	for j = 1:length(STRIP_WIDTH)			
		CORNERS[(j-1) * 4 + 1].xsign = 1.0;
		CORNERS[(j-1) * 4 + 1].ysign = 1.0;
			
		CORNERS[(j-1) * 4 + 2].xsign = 1.0;
		CORNERS[(j-1) * 4 + 2].ysign = -1.0;
		
		CORNERS[(j-1) * 4 + 3].xsign = -1.0;
		CORNERS[(j-1) * 4 + 3].ysign = 1.0;
		
		CORNERS[(j-1) * 4 + 4].xsign = -1.0;
		CORNERS[(j-1) * 4 + 4].ysign = -1.0;
	end

	for j=1:length(CORNERS)
		#println("j=",j)
		if CORNERS[j].xsign == -1
		
			CORNERS[j].x = BOARD_BORDER_X 				   - CORNERS[j].xsign*CORNERS[j].r ;
		else
			CORNERS[j].x = BOARD_BORDER_X + CHANNEL_LENGTH - CORNERS[j].xsign*CORNERS[j].r ;
		end
		
		if CORNERS[j].ysign == -1
			CORNERS[j].y =  STRIPS_BORDERS[ceil(j/4),2] - CORNERS[j].ysign* CORNERS[j].r ;
		else
			CORNERS[j].y =  STRIPS_BORDERS[ceil(j/4),1] - CORNERS[j].ysign* CORNERS[j].r ;
		end
		
		#println("x = ",CORNERS[j].x ,"   ; y = ",CORNERS[j].y)

	end
		
	WALLS = wall_t[];
	for j=1:length(STRIP_WIDTH)
		# vertical lower wall
		temp = wall_t( 	[ BOARD_BORDER_X+CORNER_RADIUS 					, STRIPS_BORDERS[j,1] ] , 
						[ BOARD_BORDER_X+CHANNEL_LENGTH-CORNER_RADIUS 	, STRIPS_BORDERS[j,1] ] ,
						  WALL_FORCE_CUTOFF);
		push!(WALLS,temp);
		
		# vertical upper wall
		temp = wall_t( 	[ BOARD_BORDER_X+CHANNEL_LENGTH-CORNER_RADIUS	, STRIPS_BORDERS[j,2] ] , 
						[ BOARD_BORDER_X+CORNER_RADIUS 					, STRIPS_BORDERS[j,2] ] , 
						  WALL_FORCE_CUTOFF);
		push!(WALLS,temp);
		
		# horizontal lower wall
		temp = wall_t( [ BOARD_BORDER_X+CHANNEL_LENGTH , STRIPS_BORDERS[j,1] - CORNER_RADIUS 			] , 
					   [ BOARD_BORDER_X+CHANNEL_LENGTH , STRIPS_BORDERS[j,1] - STRIP_DELIMITER_WIDTH/2 	] ,
						WALL_FORCE_CUTOFF);
		push!(WALLS,temp);
		
		# horizontal upper wall
		temp = wall_t( 	[ BOARD_BORDER_X+CHANNEL_LENGTH , STRIPS_BORDERS[j,2] + STRIP_DELIMITER_WIDTH/2 ] ,
						[ BOARD_BORDER_X+CHANNEL_LENGTH , STRIPS_BORDERS[j,2] + CORNER_RADIUS 			] ,
						 WALL_FORCE_CUTOFF);
		push!(WALLS,temp);
	end
		
end



