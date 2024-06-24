include("Printings.jl")

function find_coast_line(right_side::Bool, particles::Array{particle_t,1}, iteration::Int, cline_file::IOStream, prev_cl::Array{Int,1})
        if PERIODIC_X || BOX
             return (Int[] ,Int[])
        end

        coast_line = Int[]; not_cl=Int[];
        max_cl_size = Y_LEN*6;
        try_num = 1;
        side = right_side ? "right" : "left"
       # if CIRCULAR_INITIAL_LOCATION
        #    not_cl = Int[] 
        #else
            not_cl = find_small_detached_groups(particles);
        #end
        
        prev_coast_line = Int[]; # DEBUG
        while true
            i=2 ; coast_line = Int[];
            p = find_the_first_particle(particles,not_cl,right_side); # first one is just the maximum x coordinate that has at leat 3 n.n.
            if ( p == 0)
                    close_all_files()
                    error("i=$iteration: $(side)_cl: try_num=$try_num: Error in find the first particle!"); return coast_line;
            end
            #println("first p is $p")
            push!(coast_line,p)
            ( p, n_idx ) = find_the_second_particle( particles , p, right_side);

            push!(coast_line,p);
            (p ,n_idx) = find_next_neighbour( particles , coast_line,n_idx);

            while ( i < max_cl_size && length(find(coast_line.==p))<2 )
                    push!(coast_line,p);
                    (p ,n_idx) = find_next_neighbour(particles,coast_line,n_idx);
                    i+=1;
                    # check if its the end of the coast line
                    idx = findfirst(coast_line.==p)
                    if  i > 0.2*Y_LEN  && idx ==2 && coast_line[1]== coast_line[end]   # ( idx > 1 &&  coast_line[idx-1]==coast_line[end] ) || ( idx > 0 && coast_line[idx+1]==coast_line[end] )
                          pop!(coast_line);  i-=1;
                          break;
                    end
            end

            #check for errors
            if ( i == max_cl_size)
                    pr("max_cl_size = ",max_cl_size,"   YLEN = ",Y_LEN)
                    println("Error in coast line (max length exeeded) iteration ",iteration);
            end

            if  ~CIRCULAR_INITIAL_LOCATION && is_agglomeration(particles,coast_line)  &&  try_num < 6
                    inp = inpolygon(particles,coast_line);
                    pr("cl_$(side) try $(try_num): distinct group of particles size($(length(inp))) detected.")
                    pr("cline : ",coast_line', "group: ", inp')

                    for c in coast_line
                            if findfirst(inp,c) == 0
                                    error("$c is not included!!!")
                                    push!(inp,c);
                            end
                    end

                    #if ~check_if_closed_group(particles, inp,coast_line )
                            #error("not closed group");
                            #close_all_files()
                            #break;
                    #end

                    not_cl = [not_cl; inp];
                    try_num += 1;
            else
                    break;
            end
        end
        #if right_side print("right ") else print("left ") end
        #println("cl_size = ",length(coast_line)) ##DEBUG
        #delete_loops(coast_line);

        if iteration%PRINT_STEP==0 prarr(cline_file,coast_line) end
        return ( coast_line , not_cl) ;
end

function check_if_closed_group(particles::Array{particle_t,1}, group::Array{Int,1}, cl::Array{Int,1})
    # check all particles in group
    for inpk in group
            for nbr in particles[inpk].nbrs
                    if findfirst(group,nbr) == 0 # if nbr not in group
#                            if ~FIX
#                                close_all_files()
#                            else
#                                pr("cl : ",cl');
#                                close(clrID); close(cllID); close(xID); close(yID); #DEBUG
#                            end
                            pr("Warning: not closed group: $nbr is nbr in group , but does not belong to group")
                            return false
                    end
            end
    end

    # check all particles not in group
    for k=1:length(particles)
       if findfirst(group,k) == 0
              p = particles[k] ;
              for nbr in p.nbrs
                if findfirst(group,nbr) != 0 # if nbr in group
                        pr("Warning: i = $i: not closed group: $nbr is nbr of $k which is not in group , but $nbr belongs to group")
                        return false;
                end
              end
       end
    end
    return true

end

function find_small_detached_groups(particles::Array{particle_t,1});
    detached = Int[];
    for k=1:length(particles)
        if findfirst(detached,k)==0 #if group is not in the detached array yet
            if ( length(particles[k].nbrs) == 0 )
                push!(detached,k);
                pr("Small detached group detected (1): \n",k)
            elseif length(particles[k].nbrs) == 1
                    nbr = particles[k].nbrs[1];
                    if length( particles[nbr].nbrs ) ==1
                            if particles[nbr].nbrs[1] != k
                                    pr("particle to check = $k")
                                    pr("particle only nbr is $nbr")
                                    pr("nbr.nbrs[1] = ",nbr.nbrs[1])
                                    error("Error")
                            end
                            push!(detached,k);
                            push!(detached,nbr);
                            print("Small detached group detected (2) : ",[k ; nbr]')
                    end
            elseif length(particles[k].nbrs) == 2
                  nbr1 = particles[k].nbrs[1]
                  nbr2 = particles[k].nbrs[2]

                  if ( findfirst( particles[nbr1].nbrs , nbr2 ) != 0 &&
                       length(particles[nbr1].nbrs) == 2 &&
                       length(particles[nbr2].nbrs) == 2  ) ||
                     ( findfirst( particles[nbr1].nbrs , nbr2 ) == 0 &&
                       length(particles[nbr1].nbrs) == 1 &&
                       length(particles[nbr2].nbrs) == 1  )

                       detached = [ detached ; k ; nbr1 ; nbr2 ]
                       print("Small detached group detected (3) \n",[k ; nbr1 ; nbr2])
                  end
            end
        end
    end

    return detached  ;
end
## find the particle with maximum x coordinane, that has at least three nbr (eliminating dissattached groups)
function find_the_first_particle(particles::Array{particle_t,1},not_cl::Array{Int,1}, right_side::Bool)
        prev_cl = Int[]; #CHANGE DEBUG
        if right_side
            if isempty(prev_cl)
                max_x=0; idx=0;
                for i = 1:length(particles)
                    if  particles[i].x > max_x && length(particles[i].nbrs) >= 2 && findfirst(not_cl,i)==0  #find the particle with max x coord that has at least 2 n.n.
                           idx = i ;
                           max_x = particles[i].x;
                    end
                end
            else
                min_p = findmin(particles,prev_cl);
                curv = find_curvature(prev_cl);
            end
        else
             if isempty(prev_cl)
                 min_x=BOARD_BORDER_X; idx=0;
                 for i = 1:length(particles)
                     if  particles[i].x < min_x && length(particles[i].nbrs) >= 2 && findfirst(not_cl,i)==0  #find the particle with max x coord that has at least 2 n.n.
                           idx = i ;
                           min_x = particles[i].x;
                     end
                 end
             else

             end
        end

        return idx;
end

function find_the_second_particle(particles, idx::Int, right_side::Bool)
        min_angle = right_side ? -half_pi:half_pi

        j =  findfirst( particles[idx].angles .> min_angle ); #DEBUG
        if ( j == 0 ) ; j = 1;  end    # keeping in mind that the first in cline has at least two nbrs
                                       # if all the nbrs has negative angles -> the first angle is the smallest one
        return  particles[idx].nbrs[j], j
end


## idx current particle index
## n_idx - Input: prev particle index in nbrs array
##         Output :  particle index in nbrs array
function find_next_neighbour(particles::Array{particle_t,1}, cline::Array{Int,1} , n_idx::Int )
        p_idx = cline[end];
        to_print = false #p_idx==1758
        if to_print println("start $p_idx ") end

        # find prev_angle
        #j_temp = findfirst( particles[cline[end]].nbrs, cline[end-1] )
        #prev_angle = particles[cline[end]].angles[j_temp];
        prev_angle = particles[cline[end-1]].angles[n_idx];
        prev_angle += prev_angle > 0 ? -pi:pi

        if ( cline[end] == cline[end-1] ) error(); prev_angle = 0 ;  end  ##TBD: not suppose to be here?

        # find index of next_angle
        j = findfirst(particles[p_idx].angles .> prev_angle+0.00002 )
        if ( j == 0 )  j = 1 end
        if to_print; pr("prev_angle = $prev_angle ;  j= $j"); end

        # check if that index is crossing the coast line
        d_angle = particles[p_idx].angles[j] - prev_angle;
        if ( d_angle < 0 ) ; d_angle += 2*pi; end

        # 1. find the first particle that is not nbr
        k=length(cline)-1
        while ( k > 1 && sqr_dist(particles,p_idx,cline[k]) < SQR_NEAR_NBR ) k-=1 ; end
            # 2. check if there is an intersection between the ptcl that is nbr and the one that is not. ( both should not be real n.n. )
           while ( (k>1 && d_angle < CriticalAngle_Rad && sqr_dist(particles,particles[p_idx].nbrs[j],cline[end-1]) < SQR_NEAR_NBR )   )
                #just validity check
#                if !in_cline && !get_line_intersection(particles,p_idx,particles[p_idx].nbrs[j],cline[k+1],cline[k])
#                     println("Debug: interesting case, find_next_nbr in coast line p[",p_idx,"]"," j= ",j, " nbr[j]= ",particles[p_idx].nbrs[j],"line: from ",cline[k+1]," to ", cline[k]);
#                     #break
#                end
                if to_print println("checked $(particles[p_idx].nbrs[j]), going to next j (currj= $j)") end

                #increment j
                j = increment_nbr_index(particles[p_idx],j);
                if to_print println("checking $(particles[p_idx].nbrs[j]) (j=$j)") end

                #check that we haven't chacked all the nbrs already
                if ( particles[p_idx].nbrs[j]== cline[end-1]) ; println("find_next_neighbour: p[$p_idx] no suitble n.n."); break; end

                #calc new angle
                d_angle = particles[p_idx].angles[j] - prev_angle;
                if ( d_angle < 0 ) ; d_angle +=2*pi; end
           end

        if to_print println("end $p_idx. return: $(particles[p_idx].nbrs[j])") end

        return ( particles[p_idx].nbrs[j] , j );
end

# check if idx in cline, but its not the end of cline.
function is_in_cline(cline::Array{Int,1},p)
        idx = findfirst(cline.==p)
        if ( idx > 1 &&  cline[idx-1]==cline[end] ) || ( idx > 0 && cline[idx+1]==cline[end] )
              return false
        else
              return idx!=0
        end
end

function increment_nbr_index(particle::particle_t , n::Int)
        if ( length(particle.nbrs) == n )
                return 1 ;
        else
                return n+1;
        end
end

function decrement_nbr_index(particle::particle_t , n::Int)
        if ( n == 1 )
                if (length(particle.nbrs) == 0 ) ; println("Error Decrementing : particle has no neibours!\n"); end
                return length(particle.nbrs) ;
        else
                return n-1;
        end
end

function is_agglomeration(particles::Array{particle_t,1}, cl::Array{Int,1})
      dif = Array(Float64, size(cl));
      for i = 1:length(cl)-1
           dif[i] = particles[cl[i+1]].y - particles[cl[i]].y
      end
      dif[end] = particles[cl[1]].y - particles[cl[end]].y
      n = length( find( abs(dif) .> HALF_BOARD_BORDER_Y  ))

      return iseven(n)
end

function inpolygon( particles::Array{particle_t,1}, cl::Array{Int,1} )

          x = Array(Float64, length(particles));
          y = Array(Float64, length(particles));

          if length(cl) < 3
               return false & Array(Bool,length(particles))
          end

          for i=1:length(particles)
              x[i] = particles[i].x ;
              y[i] = particles[i].y ;
          end

          if PERIODIC_Y
              # Segments defined by 'border' crossings
              # find the minimum y value in the "upper" segment
              # needed for the case when part of the agglomaration is in upper side
              # and the second is in the down side
              dif =diff(y[ [cl,cl[1]] ]);
              idx = find( abs(dif) .> HALF_BOARD_BORDER_Y  )
              if !isempty(idx)
                    # check if the first segment is in upper part -> if maximum(y) is close to BOARD_BORDER

                    segment = cl[  (idx[1]+1) : idx[2]  ];
                    start_j = ( maximum(y[segment]) > (BOARD_BORDER_Y - NEAR_NBR) ) ? 1:2 ;

                    #find minimum
                    ymin = BOARD_BORDER_Y;

                    for j = start_j:2:length(idx)-1
                            segment = cl[  (idx[j]+1):idx[j+1]   ]
                            if minimum( y[segment] ) < ymin
                                 ymin = minimum(y[segment]);
                            end
                    end
                    if start_j == 2 # second segment is upper, then 0 segment is upper too
                        segment = [ cl[  (idx[end]+1):end  ] ; cl[1:idx[1]] ]
                        if minimum(y[segment]) < ymin
                                ymin = minimum(y[segment]);
                        end
                    end

                    # "move" all y under ymin, to BOARD_BORDER_Y
                    y[ y .< ymin ]  +=  BOARD_BORDER_Y;
                end
          end

          xv = x[[cl;cl[1]]];
          yv = y[[cl;cl[1]]];
          Nv = length(xv);

          # check only relevant particles
          mask = (x.>= minimum(xv)) & (x.<= maximum(xv)) & (y.>=minimum(yv)) & (y.<=maximum(yv));
          if !any(mask)
              error("mask is Empty!") #DEBUG
              return mask;
          end
          inbounds = find(mask);
          x = x[mask]';
          y = y[mask]';

          #vecorize
          Np = length(x);
          x = x[ones(Nv),:];
          y = y[ones(Nv),:];

        # Compute scale factors for eps that are based on the original vertex
        # locations. This ensures that the test points that lie on the boundary
        # will be evaluated using an appropriately scaled tolerance.
        # (m and mp1 will be reused for setting up adjacent vertices later on.)
        m = 1:Nv-1;
        mp1 = 2:Nv;
        avx = abs( 0.5*(  xv[m,:] + xv[mp1,:]) );
        avy = abs( 0.5*(  yv[m,:] + yv[mp1,:]) );
        scaleFactor = max(avx[m], avy[m]);
        scaleFactor = max(scaleFactor, avx[m,:].*avy[m,:] );
        # Translate the vertices so that the test points are
        # at the origin.
        xv = xv[:,ones(Np)] - x;
        yv = yv[:,ones(Np)] - y;

        # Compute the quadrant number for the vertices relative
        # to the test points.
        posX = xv .> 0.0;
        posY = yv .> 0.0;
        negX = !posX;
        negY = !posY;
        quad = (negX & posY) + 2*(negX & negY) + 3*(posX & negY);

        # Ignore crossings between distinct edge loops that are separated by NaNs
        nanidx = isnan(xv) | isnan(yv);
        quad[nanidx] = NaN;

        # Compute the sign() of the cross product and dot product
        # of adjacent vertices.
        theCrossProd = xv[m,:] .* yv[mp1,:] - xv[mp1,:] .* yv[m,:];
        signCrossProduct = sign(theCrossProd);

        dotProduct = xv[m,:] .* xv[mp1,:] + yv[m,:] .* yv[mp1,:];

        # Compute the vertex quadrant changes for each test point.
        diffQuad = diff(quad);

        # Fix up the quadrant differences.  Replace 3 by -1 and -3 by 1.
        # Any quadrant difference with an absolute value of 2 should have
        # the same sign as the cross product.
        idx = (abs(diffQuad) .== 3);
        diffQuad[idx] = -diffQuad[idx]/3;
        idx = abs(diffQuad) .== 2;
        diffQuad[idx] = 2*signCrossProduct[idx];

        # Find the inside points.
        # Ignore crossings between distinct loops that are separated by NaNs
        nanidx = isnan(diffQuad);
        diffQuad[nanidx] = 0;
        inn = (sum(diffQuad,1) .!= 0); inn = inn[:];


        # Find the points on the polygon.  If the cross product is 0 and
        # the dot product is nonpositive anywhere, then the corresponding
        # point must be on the contour.
        on = any((signCrossProduct .== 0) & (dotProduct .<= 0),1); on = on[:];

        inn = inn | on;
        mask[inbounds[!inn]] = false;
        return find(mask.==true)
 end

 function get_line_intersection( particles::Array{particle_t,1}, i0::Int , i1::Int , i2::Int , i3::Int )
     p0_x = particles[i0].x ;  p0_y = particles[i0].y;
     p1_x = particles[i1].x ;  p1_y = particles[i1].y;
     p2_x = particles[i2].x ;  p2_y = particles[i2].y;
     p3_x = particles[i3].x ;  p3_y = particles[i3].y;v::Array{Float64,2},eta::Array{Float64,2}

     if PERIODIC_Y
             if ( p0_y-p1_y > HALF_BOARD_BORDER_Y)
                     p1_y += BOARD_BORDER_Y;
             elseif ( p1_y-p0_y > HALF_BOARD_BORDER_X )
                     p0_y +=BOARD_BORDER_Y;
             end

             if ( p2_y-p3_y > HALF_BOARD_BORDER_Y )
                     p3_y += BOARD_BORDER_Y;
             elseif ( p3_y-p2_y > HALF_BOARD_BORDER_Y  )
                     p2_y +=BOARD_BORDER_Y;
             end
     end

     s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
     s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

     s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
     t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

     if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
         # intersection detected - cal the intersection coords
             #    i_x = p0_x + (t * s1_x);
             #    i_y = p0_y + (t * s1_y);
             #if PERIODIC_Y && i_y > BOARD_BORDER_Y
             #      i_y -= BOARD_BORDER_Y;
             #elseif PERIODIC_X && i_x > BOARD_BORDER_X
             #       i_x -= BOARD_BORDER_Y;
             #end
         return true;
     end

     return false; # No collision
 end

function delete_loops(cline::Array{Int,1} )
        i=1;
        while i<=length(cline)
            # check if loop
            n = findfirst(cline[(i+1):end],cline[i]);
            if n!=0 && iseven(n)
               if n < ( length(cline) - n ) #choose the shorter loop
                       if cline[i:i+n/2-1] == flipud(cline[i+n/2+1:i+n]) # check if a polyndrom
                               splice!(cline,i+1:i+n)
                       end
                       i+=1;
               else
                       poli = [ cline[i+n+1:end] ; cline[1:i] ];
                       if poli[1:floor(end/2)] == flipud(poli[ceil(end/2)+1:end]) # check if a polyndrom
                           splice!(cline,i+n+1:length(cline))
                           splice!(cline,1:i)
                           i=1
                       else
                           i+=1;
                       end
               end
            else
               i+=1
            end
        end
#		println(cline')
end
