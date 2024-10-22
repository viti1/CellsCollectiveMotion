include("update_map.jl")
include("Printings.jl")
include("smooth.jl")

# calculate force using a potential
function calc_force_between_paritcles!(r_cur::Array{particle_t,1}, j::Int, i::Int, angle::Float64 , F::Array{Float64,1})
        r = sqrt(sqr_dist(r_cur,j,i));
        Fr = -(U0/(A0*A0)) * ( 2*r*exp(-(r/A0)^2) + U2*exp(-r/A2) );
        if ( r > A1 ) ; Fr += 2*U1*(r-A1); end
        if ( r > A3 ) ; Fr += U3*(r-A3)^2; end
        F[1] = Fr*cos(angle);
        F[2] = Fr*sin(angle);
end

AMP  = 300
#function  stencil_circle_force!(p::particle_t, F::Array{Float64,1} )
#        CX  = BOARD_BORDER_X/2
#        CY  = BOARD_BORDER_Y/2
#        R0  = 100


#        dx = p.x - CX;
#        dy = p.y - CY;
#        r = sqrt( dx*dx + dy*dy );
#        Fr = 0;

#        #if (j==90) printf("Fr = %f \n",Fr);end
#        if ( r > R0 )1:size(x,1)
#                Fr = AMP*(r - R0);
#                F[1] = -Fr/r*dx;
#                F[2] = -Fr/r*dy;
#        else
#                F[1] = 0;
#                F[2] = 0;
#        end
#end

function stencil_add_x_right_border_force!( p::particle_t, F::Array{Float64,1} )
        if  p.x > BOARD_BORDER_X
                F[1] += -AMP * ( p.x - BOARD_BORDER_X ) ;
        end
end

function stencil_add_x_left_border_force!( p::particle_t , F::Array{Float64,1} )
        if  p.x < 0
                F[1] += ( - AMP * p.x ) ;
        end
end

function wall_attractive_force(r_cur::particle_t)
	if isnan(r_cur.x) || isnan(r_cur.x)
		println("wall_attractive_force: NaN coordinates ")
	end
	
	rvec = [r_cur.x ; r_cur.y];
	for j=1:length(WALLS)
		w = WALLS[j].p2 - WALLS[j].p1 ; # wall vector
		v = rvec-WALLS[j].p1;           # vector from first point of the wall to the cell needed in order to calculate distance to the wall)
		svdot = dot(v,w);
		if  svdot > 0 && svdot < dot(w,w) 
			l_vector = svdot/dot(w,w)*w - v ;
			l_dist = norm(l_vector); # distance from the wall 
			if l_dist < WALLS[j].d && l_dist > 0
				F = WALL_FORCE*l_vector/l_dist; 
				if any(isnan(F))
					pr("l_vector : " , l_vector')
					pr("w : " , w')
					pr("v : " , v')
					v';
				end
				return F;
			end
		end
	end

	
	for j=1:length(CORNERS)
		if is_in_corner_force_influence(r_cur,CORNERS[j])
			l_vector = [CORNERS[j].x ; CORNERS[j].y] - rvec;  # vector from the cell to the corner circle center
			l_dist = norm(l_vector); # distance from the wall 
			if l_dist < CORNERS[j].r + CORNERS[j].d && l_dist > CORNERS[j].r
				F = WALL_FORCE*l_vector/l_dist; 
				if any(isnan(F))
					pr("l_vector : " , l_vector')
					pr("w : " , w')
					pr("v : " , v')
					v';
				end

				return F;
			end
		end
	end
	
	F = [0.0 ; 0.0]
	return F;
end

function calc_curvature(particles, cl)
        #check
        if length(cl) < 3
                return (zeros(length(cl),2) , zeros(length(cl)) , zeros(length(cl)),
                        zeros(length(cl),2) , zeros(length(cl)) , zeros(length(cl)),
                        zeros(length(cl),2) , zeros(length(cl))  );
        end

        #initialize
        x = Array(Float64,length(particles));
        y = Array(Float64,length(particles));
        i=1
        for p in particles
            x[i] = p.x ;
            y[i] = p.y ;
            i+=1;
        end

        cl_x = x[cl];
        cl_y = y[cl];

        # calc tangent vector ( x' )
        t_vec = [ diff(cl_x)  diff(cl_y)  ] ;
        t_vec = [t_vec ; [ ( cl_x[1]-cl_x[end])  (cl_y[1]-cl_y[end]) ] ] ;
        if PERIODIC_Y
            idx = find( t_vec[:,2] .< - HALF_BOARD_BORDER_Y);
            t_vec[idx,2] = BOARD_BORDER_Y * ones(length(idx)) - abs( t_vec[idx,2] );
            idx = find(t_vec[:,2] .> BOARD_BORDER_Y/2);
            t_vec[idx,2] = t_vec[idx,2]-BOARD_BORDER_Y*ones(size(idx));
        end

        #normalize t_vec ( x'/|x'| )
        ds = sqrt(t_vec[:,1].^2+t_vec[:,2].^2);
        t_vec = t_vec./[ds ds];

        #calc curvature vector (dt/dt)
        ds_avg = 0.5*( [ds[end]; ds[1:end-1]] + ds );
        H = [ t_vec[1,:]-t_vec[end,:] ; diff(t_vec)];
        H = H./[ds_avg ds_avg];
        np = t_vec + [t_vec[end,:] ; t_vec[1:end-1,:]] ;
        H_sign = ( (H[:,1].*np[:,2]-H[:,2].*np[:,1]) .> 0 )*2 - 1 ;
        H_abs = sqrt(H[:,1].^2+H[:,2].^2) ;

        H_normal = H./[H_abs  H_abs ];


        if STRIPS
			half_channel = (BOARD_BORDER_X + CHANNEL_LENGTH/2);
		
           	outside_strip = collect(cl_y .< STRIPS_BORDERS[1,1])  & collect(cl_x .< half_channel) ;

			for j=1:size(STRIPS_BORDERS,1)-1
               outside_strip =  outside_strip | ( ( collect(cl_y .> STRIPS_BORDERS[j,2]) & collect(cl_y .< STRIPS_BORDERS[j+1,1]) ) & collect(cl_x .< half_channel) ); 
           	end
			outside_strip = outside_strip | ( collect(cl_y .> STRIPS_BORDERS[end,2])  & collect(cl_x .< half_channel ) );
		
           H_abs[outside_strip] = 0;
           H[outside_strip,:] = 0 ;
        end

        if GRADUAL_STRIP
           outside_strip =  ( collect(cl_y .<= STRIPS_BORDERS[1,1]) | collect(cl_y .>= STRIPS_BORDERS[1,2]) ) &  collect(cl_x .<= BOARD_BORDER_X ) ;
           H_abs[outside_strip] = 0;
           H[outside_strip,:] = 0 ;
        end

        if any(H_abs.==0.0)
            H_normal[H_abs.==0,1] = 0 ;
            H_normal[H_abs.==0,2] = 0 ;
        end

        # smoothing
        Hmod_orig = H_abs.*H_sign;
        Hmod_orig[Hmod_orig .> F_H_MAX] = F_H_MAX;
        Hmod_orig[Hmod_orig .< -F_H_MAX] = -F_H_MAX;
        Hmod = smooth(H_abs.*H_sign,ds);
        H = [Hmod Hmod] .* [H_sign H_sign] .* H_normal

        if STRIPS || GRADUAL_STRIP
             H[outside_strip,:] = 0 ;
             H_abs[outside_strip] = 0;
        end

        # second derivative of H and its smoothing
        dd_Hmod_orig  = ds_double_derivative(Hmod, ds);
        dd_Hmod = smooth(dd_Hmod_orig[:],ds);
        ddH = [H_sign H_sign].* [ dd_Hmod dd_Hmod ] .* H_normal ;

        if STRIPS || GRADUAL_STRIP
             ddH[outside_strip,:] = 0 ;
             dd_Hmod[outside_strip] = 0;
        end

        nanidx = findfirst(isnan(Hmod));
        if nanidx!=0 pr("H_mod[$nanidx] is NaN ! cl[$nanidx] = ",cl[nanidx]) end
        nanidx = findfirst(isnan(dd_Hmod));
        if nanidx!=0 pr("dd_Hmod[$nanidx] is NaN ! cl[$nanidx] = ",cl[nanidx]) end
        return (H, Hmod, Hmod_orig , ddH, dd_Hmod, dd_Hmod_orig , t_vec, ds)
end

function ds_double_derivative( arr::Array{Float64}, ds::Array{Float64,1} )

    first_der  = [ diff(arr) ;  arr[1,:]-arr[end,:] ]
    for s = 1:size(arr,2) ; first_der[:,s] = first_der[:,s]./ds ; end

    ds_avg = 0.5*( [ ds[end]; ds[1:end-1] ] + ds );
    second_der = [ first_der[1,:]-first_der[end,:] ; diff(first_der)];

    for s = 1:size(arr,2)
           second_der[:,s] = second_der[:,s]./ ds_avg
    end
    return second_der
end


function restoring_force(cl::Array{Int,1},H::Array{Float64,2}, ddH::Array{Float64,2} )
    H_abs2 = H[:,1].^2 + H[:,2].^2;
    ret = ( - 3/2*F_KAPPA*H_abs2.*H - F_KAPPA*ddH );

    for i=1:length(cl)
      if length(find(cl.==cl[i])) > 1
          ret[i,:] = 0;
      end
    end
    return ret'
	ret'
end

function Fcell_force(cl::Array{Int,1}, H::Array{Float64,2}, H_mod::Array{Float64,1})
    H_norm = H ./ [ H_mod H_mod ];
    H_norm[H_mod.==0,:] = 0;
    ret = Array(Float64,size(H));

    for i=1:length(H_mod)
        if H_mod[i] > 0
             ret[i,:] = F_CELL_MIN * H_norm[i,:]
        elseif H_mod[i] > -F_H_MAX
             ret[i,:] = (F_CELL_MIN - F_ALPHA*H_mod[i])* H_norm[i,:]
        else
             ret[i,:] = F_CELL_MAX * H_norm[i,:]
        end
        if isnan(H_norm[i,1]) || isnan(H_norm[i,2])
          error("(H_norm[$i]=",H_norm[i,:],"  H_mod = ",H_mod[i], "  H[1]=", H[i,:])
        end

        # in case of 1-cell-width finger
        if length(find(cl.==cl[i])) > 1
          ret[i,:] = 0;
        end
    end
    return ret';
	ret'
end


function Fcord_force(H_mod::Array{Float64,1},t::Array{Float64,2},cline::Array{Int,1},particles::Array{particle_t,1})
        f_cord =zeros(size(t))
        if ~CORDS
            return f_cord' ;				
        end

        len = length(H_mod);
        prev = [len ; 1:(len-1)]
        next = [2:len ; 1]
        if CORDS_STYLE == 1 #nir
            for i=1:len
                if H_mod[i] > 0
                      f_cord[i,:] = F_CORD_MAG*( t[i,:] - t[prev[i],:])
                else
                      if H_mod[prev[i]] > 0
                         f_cord[i,:] = -F_CORD_MAG*t[prev[i],:]
                      end
                      if H_mod[next[i]] > 0
                         f_cord[i,:] += F_CORD_MAG* t[i,:]
                      end
                end
            end
        else #vika
            for i=1:len
                if H_mod[i] > 0 || ( H_mod[prev[i]] > 0 && H_mod[next[i]] > 0 )
                      f_cord[i,:] = F_CORD_MAG*( t[i,:] - t[prev[i],:])
                end
            end
        end

        if CORDS_RIP
#            ( Xmin ,Xmax,Ymin,Ymax ) = RIP
#            for j=1:length(cline)
#                p = r_cur[cline[j]];
#                if   ( Xmin < p.x < Xmax ) && ( Ymin < p.y < Ymax )
#                    f_cord[j,:] = 0
#                end
#            end
             possible_rip_particles = RIP_P;

             for rp in RIP_P
                   if rp <= length(particles)
                        possible_rip_particles  = [possible_rip_particles;particles[rp].nbrs];
                   else
                        pr("Warning: RIP_P contains unexisting particles")
                   end
             end

             for j=1:length(cline)
                     if cline[j] in possible_rip_particles
                        f_cord[j,:] = 0;
                        if f_cord[prev[j],1]!=0 #change it to be projection on t[j]
                             f_cord[prev[j],:] -= F_CORD_MAG*t[prev[j],:]
                        end
                        if f_cord[next[j],1]!=0
                             f_cord[next[j],:] += F_CORD_MAG* t[j,:]
                        end
                     end
             end
        end

        return f_cord';  # array of 2D vectors
end

#  =========================== VICSEK MOVEMENT ==========================================
function cells_movement!(r_cur::Array{particle_t,1},
                          r_new::Array{particle_t,1},
                          cline_r::Array{Int,1},
                          cline_l::Array{Int,1},
                          f_files_bulk::Array{IOStream,1},
                          f_files_r::Array{IOStream,1},
                          f_files_l::Array{IOStream,1},
                          iteration:: Int )

        print_iter = (iteration-1)%PRINT_STEP==0
        #for i = 1:length(r_cur)
        #    checkinput(i,r_cur[i],"cur")
        #end

        # ----- allocate force vars -----------------------
        sum_dv      = Array(Float64,2);
        fij         = Array(Float64,2);
        f_interaction = Array(Float64,2);
        f_stencil   = Array(Float64,2);
        #f_friction  = Array(Float64,2);
        #f_noise     = Array(Float64,2);
        #f_vicsek    = zeros(2); 
		Hmodr=Float64[]; Hmodl=Float64[];
        # ---- calculate cline forces --------------------

        if  ~isempty(cline_r)
            # H[x/y,partice_idx]
            if  iteration < THERMALIZATION_TIME || BOX
                (H,  ddH , t_vec) = init_arrays_2D_to_zero(3,length(cline_r))
                ( Hmod,Hmod_orig,ddHmod,ddHmod_orig,ds,Hmodr ) = init_arrays_1D_to_zero(6,length(cline_r));
                ( Fres_r , Fcord_r , Fcell_r ) = init_arrays_2D_to_zero(3,length(cline_r))
            else

                ( H, Hmod,Hmod_orig, ddH ,ddHmod, ddHmod_orig, t_vec,ds) = calc_curvature(r_cur, cline_r)
                Fres_r = restoring_force(cline_r,H,ddH)
                Fcell_r = Fcell_force(cline_r,H,Hmod)
                Fcord_r = Fcord_force(Hmod,t_vec,cline_r,r_cur)
                Hmodr = Hmod;
            end

            if print_iter
                print_cl_forces(f_files_r,
                                Fres_r, Fcell_r, Fcord_r,
                                H, Hmod, ddH,
                                ds, ddHmod, Hmod_orig, ddHmod_orig )
            end
        end

        if  ~isempty(cline_l)
            pr(cline_l)
            if  iteration < THERMALIZATION_TIME || BOX
                (H,  ddH , t_vec) = init_arrays_2D_to_zero(3,length(cline_l))
                ( Hmod,Hmod_orig,ddHmod,ddHmod_orig,ds, Hmodl) = init_arrays_1D_to_zero(6,length(cline_l));
                ( Fres_l , Fcord_l , Fcell_l ) = init_arrays_2D_to_zero(3,length(cline_l))
            else
                ( H, Hmod,Hmod_orig, ddH ,ddHmod,ddHmod_orig, t_vec,ds) = calc_curvature(r_cur, cline_l)
                Fres_l = restoring_force(cline_l,H,ddH)
                Fcell_l = Fcell_force(cline_l,H,Hmod)
                Fcord_l  = Fcord_force(Hmod,t_vec,cline_l,r_cur)
                Hmodl = Hmod;
            end

            if print_iter
                print_cl_forces(f_files_l,
                                Fres_l, Fcell_l, Fcord_l,
                                H, Hmod, ddH)
                                #ds, ddHmod, Hmod_orig, ddHmod_orig)
            end
        end

        # ~~~~~~~~|  LOOP over All Particles |~~~~~~~~~
        #``````````````````````````````````````````````

        for ip = 1:length(r_cur)
                p_cur = r_cur[ip];
                p_new = r_new[ip];
				
				# -----  f1 - friction force --------------------------------------
                f_friction[:] = -ALPHA*p_cur.v;
				
                # -----  f2 - sum of (Vj-Vi) where j is the n.n. of p[i] ----------
                sum_dv[:] = 0;
				if iteration == 200 && ip > 20 && ip < 30
							println("i_particle = ",ip, " V= " , p_cur.v) 
						end
                for n_idx in p_cur.nbrs
						if iteration == 200 && ip > 20 && ip < 30
							println("dv = ",r_cur[n_idx].v - p_cur.v) 
						end
                        sum_dv += ( r_cur[n_idx].v - p_cur.v );
                end

                if ( !isempty(p_cur.nbrs))
                        f_vicsek = sum_dv * BETA / length(p_cur.nbrs);
                end

                # -------  f3 - sum of the forces applied by n.n.  ( depends on the distance between the particles) -------
                f_interaction[:] = 0;
                ind = 1;
                for n_idx in p_cur.nbrs
                        calc_force_between_paritcles!(r_cur, ip, n_idx, p_cur.angles[ind] , fij);
                        f_interaction += fij;
                        ind += 1;
                end

                # -----  f4 - random noise -------------------------
				# p_new.eta is actually the eta one iteration before p_cur.
				eta_old = p_new.eta;
                f_noise[:] = ( SIGMA0+(SIGMA1-SIGMA0)*(1-p_cur.rho/RHO0) ) * eta_old;
                p_new.eta[:] = p_cur.eta + dt_divtau*(-p_cur.eta + Noise_Coeff * randn(2));

                # ---  f5 - stencil force --------
                f_stencil[:] = 0;
                if TWO_P
                    stencil_circle_force(p_cur, f_stencil);
                else
                    # right stencil
                    if  ( iteration < THERMALIZATION_TIME || BOX ) && !PERIODIC_X && !NO_BORDERS
                            stencil_add_x_right_border_force!(p_cur, f_stencil);
                    end

                    # left stencil
                    if ( !TWO_SIDED || iteration < THERMALIZATION_TIME || BOX ) && !PERIODIC_X && !NO_BORDERS
                            stencil_add_x_left_border_force!(p_cur, f_stencil);
                    end
                end

				# ---  f6 - walls attractive force --------
				f_wall = wall_attractive_force(p_cur);
				
				# if f_wall[1] !=0 || f_wall[1] !=0
					# pr("i=",iteration,"  f_wall= ",f_wall')
					# iteration'
				# end
		
		
                # ---- calc new velocity ----
				
                p_new.v[:] = p_cur.v + DT*( f_friction + f_interaction + f_vicsek + f_noise + f_stencil + f_wall) ;
                icl = findfirst(ip.==cline_r)
                if icl > 0 # if i in cline right
                     p_new.v += ( Fres_r[:,icl] + Fcell_r[:,icl] + Fcord_r[:,icl] ) * DT
                end

                icl = findfirst(ip.==cline_l)
                if icl > 0 # if i in cline leftf
                      p_new.v += ( Fres_l[:,icl] + Fcell_l[:,icl] + Fcord_l[:,icl] ) * DT
                end

                if ( abs(p_new.v[1]) > MAX_VALUE_V )
                        p_new.v[1] = sign(p_new.v[1])*MAX_VALUE_V
                end

                if ( abs(p_new.v[2] ) > MAX_VALUE_V )
                        p_new.v[2] = sign(p_new.v[2])*MAX_VALUE_V
                end
                # --- calc new loctaion -------
                p_new.x = p_cur.x + (DT * p_new.v[1] );
                p_new.y = p_cur.y + (DT * p_new.v[2] );

                if  PERIODIC_Y
                    if ( p_new.y >= BOARD_BORDER_Y )
                            p_new.y -= BOARD_BORDER_Y
                    elseif (  p_new.y < 0 )
                            p_new.y += BOARD_BORDER_Y
                    end
                end

                if  PERIODIC_X
                    if  p_new.x >= BOARD_BORDER_X
                            p_new.x -= BOARD_BORDER_X
                    elseif   p_new.x < 0
                            p_new.x += BOARD_BORDER_X
                    end
                end

                if STRIPS || GRADUAL_STRIP; slide_to_borders!(p_new,p_cur,ip); end
				
                # --- print all forces -----------
                if print_iter
                   print_bulk_forces_of_particle(f_files_bulk,f_friction, f_vicsek, f_interaction,f_noise,f_wall  )
                end

                # --- check x and v isnan ------
                if  isnan(p_new.v[1]) || isnan(p_new.v[2])
                        println("V[$ip] = ",p_new.v)

                        print_bulk_forces_of_particle(f_friction, f_vicsek, f_interaction,f_noise,f_wall )
                        if  ~isempty(cline_r)
                              print_cl_forces_of_particle( ip,cline_r, Fres_r, Fcell_r, Fcord_r,"Right" )
                        end
                        if ~isempty(cline_l)
                              print_cl_forces_of_particle( ip,cline_l,Fres_l, Fcell_l, Fcord_l,"Left" )
                        end
                        close_all_files()
                        error("Nan Velocity!!!")
                end

                if  isnan(p_new.x) || isnan(p_new.y)
                        close_all_files()
                        error("Nan coordinate: r_new[$ip].x = ",p_new.x," r_new[$i].y = ",p_new.y )
                end
        end

        # ~~~~~~~~~~~~~~~~~~| END LOOP |~~~~~~~~~~~~~~~~~~
        #`````````````````````````````````````````````````
        if print_iter print_new_line_to_files(f_files_bulk); end
        #pr("length before update ", length(r_new) )

        update_map(r_new,r_cur,[cline_r;cline_l],[Hmodr[:];Hmodl[:]],iteration); # devide and update nearest neibours and their angles in r_cur
end


function checkinput(ip::Int, p::particle_t, s::ASCIIString)
        if isnan(p.v[1]) || isinf(p.v[1]) pr("checkinput: v_$s[$ip].x = ",p.v[1] ) end
        if isnan(p.v[2]) || isinf(p.v[2]) pr("checkinput: v_$s[$ip].y = ",p.v[2] ) end
        if isnan(p.x) || isinf(p.x) pr("checkinput: r_$s[$ip].x = ",p.x ) end
        if isnan(p.y) || isinf(p.y) pr("checkinput: r_$s[$ip].y = ",p.y ) end
end

function slide_to_borders!(particle::particle_t, old_particle::particle_t,ptcl_id::Int64)
	particle_bu = deepcopy(particle);
	was_slide = false;
	
	# update coordinates
	if STRIPS
		was_slide = slide_to_borders_coordinate!(particle, old_particle);
	end
	
	if GRADUAL_STRIP
		was_slide = slide_to_gradual_borders_coordinates!(particle, old_particle);
	end
	
	# update velocity
	was_slide2 = particle_bu.x != particle.x || particle_bu.y !=  particle.y;
	if was_slide != was_slide2
		println("pid = ",ptcl_id ,": was_slide != was_slide2 " , was_slide, " != ", was_slide2)
	
		if particle_bu.x != particle.x
			println(" old.x = " ,particle_bu.x, "new_x = " ,particle.x)
		end
		
		if particle_bu.y != particle.y
			println(" old.y = " ,particle_bu.y, "new_y = " ,particle.y)
		end
	end
	
	if was_slide
		dx = particle.x - old_particle.x ;
		dy = particle.y - old_particle.y ;        

        if ( PERIODIC_Y && abs(dy) > HALF_BOARD_BORDER_Y )
			if dy > 0							
                dy = dy - BOARD_BORDER_Y ;
			else
				dy = dy + BOARD_BORDER_Y ;
			end
        end

        if ( PERIODIC_X && abs(dx) > HALF_BOARD_BORDER_X )
             if dx > 0							
                dx = dx - BOARD_BORDER_X ;
			else
				dx = dx + BOARD_BORDER_X ;
			end
        end	

		if abs(dy) > 6
			pr("abs(dy) > 6  : ", dy , " y_old=",old_particle.y, "  y_new=" , particle.y, "  x_old=",old_particle.x, "  x_new=" , particle.x)
		end
		
		if abs(dx) > 6
			pr("abs(dx) > 6: " , dx , "x_old=",old_particle.x, "x_new" , particle.x)
		end
		
		particle.v[1] = dx/DT;
		particle.v[2] = dy/DT
	end
		
end

function slide_to_borders_coordinate!(particle::particle_t, old_particle::particle_t)
	end_of_strip = BOARD_BORDER_X + CHANNEL_LENGTH;
	
    if particle.x <= BOARD_BORDER_X || particle.x >= end_of_strip
      return false
    end

	particle_bu = deepcopy(particle);
	for i = 1:length(CORNERS)
		if is_in_corner(particle,CORNERS[i]) 
			#println("in corner " , i ,  " : (x,y) = (",particle.x,",",particle.y,")")
			dist_to_corner =  dist_to_corner_center(particle,CORNERS[i]) 
			if dist_to_corner < CORNERS[i].r
				# slide to the closest point on the circle
				# Pnew = (P - O)/|P-O| * R     dist_to_corner = |P-0|
				particle.x = CORNERS[i].x + ( particle.x - CORNERS[i].x )/ dist_to_corner  * CORNERS[i].r ;
				particle.y = CORNERS[i].y + ( particle.y - CORNERS[i].y )/ dist_to_corner  * CORNERS[i].r ;
				
				# println("slide from corner: old(x,y) = ( ",particle_bu.x,",",particle_bu.y,") new(x,y) = (",particle.x,",",particle.y,")")
				return true;
			else
				return false;
			end
		end		
	end
	
    if particle.y < STRIPS_BORDERS[1,1] 
        if old_particle.x <= BOARD_BORDER_X
                particle.x =  BOARD_BORDER_X-0.1
        elseif old_particle.x >= end_of_strip
                particle.x =  end_of_strip+0.1
						
        else
                particle.y =  STRIPS_BORDERS[1,1]
				# println("slide A")
        end
		
        return true
    end

    for j=1:size(STRIPS_BORDERS,1)-1
         if particle.y <= STRIPS_BORDERS[j,2]
            return false
         end

         if  particle.y <= STRIPS_BORDERS[j+1,1] # => particle outside strip

            if old_particle.x <= BOARD_BORDER_X
                particle.x =  BOARD_BORDER_X - 0.1
            elseif old_particle.x >= end_of_strip
                particle.x =  end_of_strip+0.1
            else

                if STRIPS_BORDERS[j+1,1] - particle.y < particle.y - STRIPS_BORDERS[j,2]
                    particle.y = STRIPS_BORDERS[j+1,1]
					println("slide B")
                else
                    particle.y = STRIPS_BORDERS[j,2]
					# println("slide C")
                end

            end
				
            return true
         end
    end

    if particle.y > STRIPS_BORDERS[end,2]
          if old_particle.x <= BOARD_BORDER_X
                particle.x =  BOARD_BORDER_X-0.1
          elseif old_particle.x >= end_of_strip
                particle.x =  end_of_strip+0.1
          else
                particle.y =  STRIPS_BORDERS[end,2]
				# println("slide D")
          end

          return true
    end
	return false
end

function slide_to_gradual_borders_coordinates!(particle::particle_t, old_particle::particle_t)
   # if behind start line
   if particle.x <= BOARD_BORDER_X
      return false
    end

    y_curr_border  = STRIPS_BORDERS[1,1:2] +  [1 -1]*STRIP_SLOPE*(particle.x - BOARD_BORDER_X)
	#println("slide_to_gradual_borders: y_curr_border = ",y_curr_border)

    #region1 - under the strip
    if particle.y < y_curr_border[1]
        if old_particle.x <= BOARD_BORDER_X
                particle.x =  BOARD_BORDER_X-0.05
        else
                particle.y =  y_curr_border[1]
        end
		return true
    end

  # region3 - above the  strip
   if particle.y > y_curr_border[2]
        if old_particle.x <= BOARD_BORDER_X
              particle.x =  BOARD_BORDER_X - 0.05
        else
              particle.y = y_curr_border[2];
        end
		return true
    end

	return false
end

function dist_to_corner_center(p::particle_t,corner::corner_t)
	return sqrt( ( p.x - corner.x)^2 + ( p.y - corner.y)^2 ) 
end

function  is_in_corner(p::particle_t,corner::corner_t)
	# if (P - O1) and (P - O2) has the opposite sign -> P is between O1 and O2
	
	# x1 = corner.x ; x2 = corner.xsign*corner.r + corner.x;
	# y1 = corner.y ; y2 = corner.ysign*corner.r + corner.y;
	# result2 =  (( p.x > x1 &&  p.x < x2 ) ||  ( p.x < x1 &&  p.x > x2 )) &&
			# (( p.y > y1 &&  p.y < y2 ) ||  ( p.y < y1 &&  p.y > y2 ))
	
	result = ( ( p.x - corner.x )* ( p.x - (corner.xsign*corner.r + corner.x) ) < 0 && 
	          ( p.y - corner.y )* ( p.y - (corner.ysign*corner.r + corner.y) ) < 0     );
			  # if result
				# println(" (corner.xsign*corner.r + corner.x) = " ,  (corner.xsign*corner.r + corner.x));
				# println(" (corner.ysign*corner.r + corner.y) = " ,  (corner.ysign*corner.r + corner.y));
				# println("particle(x,y) = (" , p.x, "," , p.y,")")
				# exit()
			  # end		
			
	# if result !=result2
		# println("result1 !=result2 : " ,result1 , " != " ,result2 )
	# end			
	return 	result	
end

function  is_in_corner_force_influence(p::particle_t,corner::corner_t)
	# if (P - O1) and (P - O2) has the opposite sign -> P is between O1 and O2
	
	# x1 = corner.x ; x2 = corner.xsign*corner.r + corner.x;
	# y1 = corner.y ; y2 = corner.ysign*corner.r + corner.y;
	# result2 =  (( p.x > x1 &&  p.x < x2 ) ||  ( p.x < x1 &&  p.x > x2 )) &&
			# (( p.y > y1 &&  p.y < y2 ) ||  ( p.y < y1 &&  p.y > y2 ))
	
	result = ( ( p.x - corner.x )* ( p.x - (corner.xsign*corner.d + corner.x) ) < 0 && 
	          ( p.y - corner.y )* ( p.y - (corner.ysign*corner.d + corner.y) ) < 0     );
			  # if result
				# println(" (corner.xsign*corner.r + corner.x) = " ,  (corner.xsign*corner.r + corner.x));
				# println(" (corner.ysign*corner.r + corner.y) = " ,  (corner.ysign*corner.r + corner.y));
				# println("particle(x,y) = (" , p.x, "," , p.y,")")
				# exit()
			  # end		
			
	# if result !=result2
		# println("result1 !=result2 : " ,result1 , " != " ,result2 )
	# end			
	return 	result	
end


function init_arrays_1D_to_zero(NumOfArrays,arrLength)
    a = Array(Array{Float64,1},NumOfArrays);
    for i=1:NumOfArrays
        a[i]= zeros(arrLength);
    end

    return tuple(a...)
end

function init_arrays_2D_to_zero(NumOfArrays,arrLength)
    a = Array(Array{Float64,2},NumOfArrays);
    for i=1:NumOfArrays
        a[i]= zeros(2,arrLength);
    end

    return tuple(a...)
end

