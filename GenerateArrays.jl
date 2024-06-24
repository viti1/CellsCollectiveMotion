# generating initial arrays (r , v , eta )
include("update_map.jl")

function generate_particles_locations() # Fil Rectangle
        # -- generate location
        particles = particle_t[]
        x = 0.1*dX + ( 0:dX:dX*(X_LEN-1) )
        y = 0.1*dY + ( 0:dY:dY*(Y_LEN-1) )
        rho = 9/((dX^2+dY^2)*pi)
        xarr = Float64[]; yarr = Float64[];
        for i= 1 : X_LEN
            for j = 1:Y_LEN
                p = particle_t(x[i]+rand()*dX*0.8 , y[j]+rand()*dX*0.8  ,[0.0;0.0],[0.0;0.0],Array(Int64,0),Array(Float64,0),rho)
                push!(particles,p);
            end
        end

        # -- find neibors ( check 8 nbrs)
        for k = 1:length(particles)
                listt = [k+Y_LEN-1,k+Y_LEN,k+Y_LEN+1, k-1,k+1, k-Y_LEN-1,k-Y_LEN,k-Y_LEN+1 ]
                if ( k % Y_LEN == 0 ) # first row
                        push!(listt,k+2*Y_LEN-1)
                elseif ( k+1 ) % Y_LEN == 0  # last row
                        push!(listt,k-2*Y_LEN+1)
                end

                for nbr in listt
                     if nbr>0 && nbr<=length(particles)
                        check_neighbor!(particles,k,nbr);
                     end
                end
        end

        return particles
end

function generate_particles_locations1() #circle
        # -- generate location
        particles = particle_t[]

        for x = 0:dX:CR
             y = 0.0 ;
             while ( x^2+(y - 0.5*dY)^2 <= CR^2)
                if x^2+( y + 0.5*dY)^2 > CR^2  # next y
                        y1 = sqrt(CR^2 - x^2)
                else
                        y1 = y
                end

                if x == CR-dX
                        x1 = sqrt(CR^2 - y1^2)
                else
                        x1 = x
                end

               

                push!(particles, particle_t(x1 + CX , y1 + CY  ,[0.0,0.0],[0.0,0.0],Int[] ,Float64[],0.0 ) )

                if y!=0
                    push!(particles, particle_t(x1 + CX, -y1 + CY ,[0.0,0.0],[0.0,0.0],Int[] ,Float64[] ,0.0) )
                end
                if x!=0
                    push!(particles, particle_t(-x1 + CX, y1 + CY  ,[0.0,0.0],[0.0,0.0],Int[] ,Float64[] ,0.0) )
                        if y!=0
                            push!(particles, particle_t(-x1 + CX, -y1 + CY ,[0.0,0.0],[0.0,0.0],Int[] ,Float64[] ,0.0) )
                        end
                end
                y+=dY
             end
        end
        NUM_OF_PARTICLES = length(particles)

        # -- find neibors
        for k = 1:NUM_OF_PARTICLES
                for nbr = 1:NUM_OF_PARTICLES
                     check_neighbor!(particles,k,nbr);
                end
        end

        return particles
end

function generate_particles_locations3() #two sides
        # -- generate location
        GAP_X = dX*20
        particles = []
        x = 0.1*dX + ( 0:dX:dX*(X_LEN/2-1) )

        x = [x ; x + dX*(X_LEN/2) + GAP_X ]

        y = 0.1*dY + ( 0:dY:dY*(Y_LEN-1) )
        rho = 9/(dX*dY)
        for i= 1 : X_LEN
            for j = 1:Y_LEN
                p = particle_t(x[i]+rand()*dX*0.8 , y[j]+rand()*dX*0.8  ,[0.0,0.0],[0.0,0.0],Array(Int64,0),Array(Float64,0),rho)
                push!(particles,p);
            end
        end

        NUM_OF_PARTICLES = length(particles);
        # -- find neibors ( check 8 nbrs)

        for k = 1:NUM_OF_PARTICLES
                nbrs_to_check = [k+Y_LEN-1,k+Y_LEN,k+Y_LEN+1, k-1,k+1, k-Y_LEN-1,k-Y_LEN,k-Y_LEN+1 ]
                if ( k % Y_LEN == 0 ) # first row
                        push!(nbrs_to_check,k+2*Y_LEN-1)
                elseif ( k+1 ) % Y_LEN == 0  # last row
                        push!(nbrs_to_check,k-2*Y_LEN+1)
                end

                for nbr in nbrs_to_check
                     check_neighbor!(particles,k,nbr);
                end
        end

        return particles
end


function read_particles(path)
        particles = particle_t[]
        last = readdlm("$(path)last_pos.txt",'\t')'
        NUM_OF_PARTICLES = size(last,1)
        pr("NUM_OF_PARTICLES=",NUM_OF_PARTICLES)

        for k = 1:NUM_OF_PARTICLES
            push!(particles,particle_t(last[k,1],last[k,2],[last[k,3],last[k,4]],[last[k,5],last[k,6]],Int[],Float64[],last[k,7]))
        end

        #find neigbors
        for k = 1:NUM_OF_PARTICLES
                difx = abs(last[:,1] - last[k,1]);
                dify = abs(last[:,2] - last[k,2]);
                if  PERIODIC_Y
                        ff = dify .> HALF_BOARD_BORDER_Y;
                        dify[ff] =  BOARD_BORDER_Y - dify[ff]
                end
                if  PERIODIC_X
                        ff = difx .> HALF_BOARD_BORDER_X;
                        difx[ff] =  BOARD_BORDER_X - difx[ff]
                end
                listt = find( (difx.<=NEAR_NBR) .* (dify.<= NEAR_NBR) )

                for nbr in listt
                     check_neighbor!(particles,k,nbr);
                end
        end

        #check that there is no cell that are not "mutual nbrs"
        for i = 1:length(particles)
                rm_list = Int[]
                ind = 1
                for nb in particles[i].nbrs
                        if ( findfirst(particles[nb].nbrs .==i ) == 0 )
                              push!(rm_list,ind)
                              #println("$nb is nbr of $i, but not the opposite");
                        else
                            ind += 1;
                        end
                end

                for rm in rm_list
                    splice!(particles[i].nbrs,rm)
                    splice!(particles[i].angles,rm)
                end

                #if !isempty(rm_list) pr("new nbrs list :", particles[i].nbrs) end
        end

        return particles
end

function  initialize_eta(len)
        return randn(len,2)*dt_divtau* Noise_Coeff
end

function initialize_arrays(start_from_last,path)

        if start_from_last
                println(" -- load location ---- ")
                r_vec = read_particles(path)
        else
                println(" -- gen  location ---- ")
				if CIRCULAR_INITIAL_LOCATION
					r_vec = generate_particles_locations1();
				else
                	r_vec = generate_particles_locations();
				end
        end

        r_vec_new = deepcopy(r_vec);
        return r_vec, r_vec_new
end
