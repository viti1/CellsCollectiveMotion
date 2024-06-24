function sqr_dist(particles,idx1, idx2)

        dy = abs(particles[idx1].y - particles[idx2].y) ;
        dx = abs(particles[idx1].x - particles[idx2].x) ;

        if ( PERIODIC_Y && dy > HALF_BOARD_BORDER_Y )
                dy = BOARD_BORDER_Y - dy ;
        end

        if ( PERIODIC_X && dx > HALF_BOARD_BORDER_X )
                dx = BOARD_BORDER_X - dx ;
        end

        return dx*dx+dy*dy;
end

function decrement_nbr_index(pnbrs,j);
        if j == 1 || j == 0  # 0 means the current j is after the last index
           return length(pnbrs) ;
        else
           return j-1
        end
end

function  increment_nbr_index(pnbrs, j )
        if ( length(pnbrs) == j )
                return 1 ;
        else
                return j+1;
        end
end

function my_abs( x )
        if x >= 0
             return x
        else
             return -x
        end
end

function check_neighbor!(particles::Array{particle_t,1}, p_idx::Int , n_idx::Int)
        if  p_idx==n_idx #|| findfirst(particles[p_idx].nbrs, n_idx) != 0
                return;
        end

        overlap1 = false;
        overlap2 = false;
        pnbrs = particles[p_idx].nbrs;
        pangles = particles[p_idx].angles;
        to_print =  false #p_idx == 1341 && ( n_idx == 1279  )  ;

        dx = particles[n_idx].x - particles[p_idx].x;
        dy = particles[n_idx].y - particles[p_idx].y;

        if ( PERIODIC_Y && abs( dy ) >   HALF_BOARD_BORDER_Y  )
                dy = dy > 0 ? -( BOARD_BORDER_Y - abs(dy))  :  BOARD_BORDER_Y - abs(dy)  ;
        end

        if ( PERIODIC_X && abs( dx ) >  HALF_BOARD_BORDER_X  )
                dx = dx > 0 ? -( BOARD_BORDER_X - abs(dx))  :  BOARD_BORDER_X - abs(dx)  ;
        end

        if dx*dx+dy*dy > SQR_NEAR_NBR
            return;
        end

        the_angle = atan2( dy , dx );

        if  isempty( pnbrs ) # p_idx has no nbrs yet
                push!(pnbrs, n_idx);  push!(pangles, the_angle);
                return;
        end

        # == compare overlap with prev nbr
        r_pow2 = dx*dx + dy*dy;
        j_prev = findfirst( pangles .>= the_angle );
        if j_prev==0; j_prev = 1; end;
        if pnbrs[j_prev] == n_idx ;  return;  end  #particles already exist in nbr array

        if ( to_print ) println("before prev  particles[$p_idx].nbrs: \n", particles[p_idx].nbrs)  end
        k = 1;
        while ( k < 10 )
            j_prev = decrement_nbr_index(pnbrs,j_prev);
            if length(pnbrs)==0; break; end

            d_angle =  the_angle - pangles[j_prev];
            if  ( d_angle < 0 ); d_angle = 2*pi + d_angle; end

            if ( d_angle < half_pi )
                ( overlap, r_nbr_pow2 ) = check_overlap(particles,p_idx,n_idx,j_prev,r_pow2,to_print) ;
                 overlap1 = overlap1 || overlap ; #TBD

                if ( overlap && r_nbr_pow2 < r_pow2 )
                        if to_print   println("p[$p_idx]: $n_idx is overlapped by $(pnbrs[j_prev]) (prev)") end
                        j_next = increment_nbr_index(pnbrs,j_prev)
                        ( next_overlap, r_nbr_pow2 ) = check_overlap(particles,p_idx,n_idx,j_next,r_pow2,to_print) ;
                        if to_print   println("$n_idx is checked on overlapping with $(pnbrs[j_next]). result =  $next_overlap") end
                        if next_overlap && r_nbr_pow2 > r_pow2
#                            println("p[$p_idx], $n_idx is overlapped by $(pnbrs[j_prev]), but also overlapping $(pnbrs[j_next]) - hance it is removed")
                            splice!(pnbrs,j_next); splice!(pangles,j_next);
                        end
                        return; # existing particle hides the new one
                elseif !overlap # no overlap
                     break;
                end
            else # no overlap for sure
                 break
            end

            #continue if there was an overlap -> looking for the next one
            splice!(pnbrs,j_prev); splice!(pangles,j_prev);
            if to_print println("p[$p_idx]: $(pnbrs[j_prev]) is overlapped(and deleted) by $n_idx (prev)"); end
        end

        if ( k==10 ) println("Error KOKO1"); end;
        if ( to_print ) println("before next  particles[$p_idx].nbrs: ", particles[p_idx].nbrs')  end


        # ==  compare overlap with next nbr
        j_next = findfirst( pangles .>= the_angle );
        if j_next ==0 ; j_next=1; end;

        if isempty(pnbrs)
               d_angle = 2*pi; #this will cause to skip the following loop
        else
               d_angle = pangles[j_next] - the_angle;
               if ( d_angle < 0 ) ; d_angle = 2*pi + d_angle ; end
        end

        k=1
        while ( d_angle < half_pi &&  k<10)
                ( overlap, r_nbr_pow2 ) = check_overlap(particles,p_idx,n_idx,j_next,r_pow2,to_print) ;
                overlap2 = overlap || overlap2; #TBD is there was an overlap at least one time its overlap

                if ( overlap && r_nbr_pow2 < r_pow2 )
                    if (to_print) println("p[$p_idx]: $n_idx is overlapped by $(pnbrs[j_next]) (next)") end
                    return; #  existing particle hides the new one

                elseif overlap
                    if (to_print) println("p[$p_idx]: $(pnbrs[j_next]) is overlapped by $n_idx (next)") end
                    splice!(pnbrs,j_next); splice!(pangles,j_next);

                    if ( j_next >  length(pnbrs) ) ; j_next = 1; end
                    if  length(pnbrs) == 0 break; end

                    d_angle = pangles[j_next] - the_angle;
                    if ( d_angle < 0 ) ; d_angle = 2*pi + d_angle ; end
                else # no overlap
                        break;
                end
                k+=1;
                #continue if there was an overlap -> looking for the next one
        end

        if ( k==10 ) println("Error p[$p_idx] KOKO2"); end;
        if ( to_print ) println("before insert particles[$p_idx].nbrs: \n", particles[p_idx].nbrs)  end

        # == insert the particle ==
        if ( !isempty(pangles) && the_angle > pangles[end] ) # n_idx should follow the last nbr in the array
             push!(pnbrs,n_idx) ; push!(pangles,the_angle)
        else
             insert!(pnbrs,j_next,n_idx);
             insert!(pangles,j_next,the_angle);
        end

        if length(pnbrs)!=length(pangles)
             println("Error in check nbr p[$p_idx] , n_idx = $n_idx: length(pnbrs)!=length(pangles) ");
        end

        if ( to_print ) print("at end particles[$p_idx].nbrs: ", particles[p_idx].nbrs') end
end


function check_overlap(particles,p_idx,n_idx,je,r_pow2,to_print)
        dist_pow2  = sqr_dist(particles, n_idx , particles[p_idx].nbrs[je] );
        r_nbr_pow2 = sqr_dist(particles, p_idx , particles[p_idx].nbrs[je] );
        r_max_pow2 = max(r_nbr_pow2,r_pow2);
        r_min_pow2 = min(r_nbr_pow2,r_pow2);
        dist =	sqrt(dist_pow2 );
        r_min = sqrt(r_min_pow2);
        cos_beta = ( r_min_pow2 + dist_pow2 - r_max_pow2 ) /  ( 2 * r_min * dist );

        if (to_print) println("beta = ",acos(cos_beta), "  CRITICAL_ANG = ", acos(COS_CRITICAL_ANGLE) ); end
        return ( cos_beta < COS_CRITICAL_ANGLE , r_nbr_pow2 ) ;
end

function add_to_list( listt::Array{Int,1} , idx::Int )
        if findfirst(listt, idx ) == 0 ; push!(listt,idx); end
end

function update_map(r_new::Array{particle_t,1},r_cur::Array{particle_t,1},cl::Array{Int,1},curv::Array{Float64,1},iteration::Int)
        listt = Int[];
        ADDITIONAL_CHECK = iteration%50 == 0;
        high_curv = cl[curv.<MIN_CURV_FOR_DIVIDE];
        # -----------  division ----------------------
        if DIVIDE
            for i = 1:length(r_new)
                 if i in high_curv continue end
                 rho = r_new[i].rho
                 div_time = MIN_DIV*(1+(rho/RHO1)^4)
                 prob = DT/div_time;
                 if rand() <  prob
                     if find_min_distance(r_new,i) > NEW_CELL_DIST
                        if PRINT_TO_STDOUT @printf("dividing cell %d :  div_time = %.1f ; prob = %lf\n",i,div_time,prob) end
                        divide_cell(r_new,r_cur,i,iteration*DT)
                     end
                 end
            end
        end

        if ADDITIONAL_CHECK || length(r_cur) < length(r_new)
                #pr("additional nbrs check")
                x = Array(Float64,length(r_new))
                y = Array(Float64,length(r_new))
                j=1;
                for p in r_new
                   x[j]  = p.x
                   y[j]  = p.y
                   j +=1;
                end
        end

        for i = 1:length(r_new)
                listt = Array(Int,0);
                # prepare list of nbrs
                if ADDITIONAL_CHECK || i>length(r_cur)
                        difx = abs(x - x[i]);
                        dify = abs(y - y[i]);
                        if  PERIODIC_Y
                                ff = dify .> HALF_BOARD_BORDER_Y;
                                dify[ff] =  BOARD_BORDER_Y - dify[ff]
                        end
                        if  PERIODIC_X
                                ff = difx .> HALF_BOARD_BORDER_X;
                                difx[ff] =  BOARD_BORDER_X - difx[ff]
                        end

                        listt = find( (difx.<=NEAR_NBR) .* (dify.<= NEAR_NBR) )
                else
                        for n_idx in r_cur[i].nbrs
                            add_to_list(listt,n_idx);
                            for nn_idx in r_cur[n_idx].nbrs
                                 add_to_list(listt,nn_idx); #check neighbours neighbours
                            end
                        end
                end

                r_new[i].nbrs = Int64[]
                r_new[i].angles = Float64[]
                for nbr in listt
                      check_neighbor!(r_new,i,nbr);
                end

                if i > length(r_cur) # only for fix
                   for n in r_new[i].nbrs
                        check_neighbor!(r_new,n,i)
                   end
                end
        end

        # ----------- check non-mutual nbrs ------------------

        for i = 1:length(r_new)
                rm_list = Int[]
                ind = 1
                for nb in r_new[i].nbrs
                        if ( findfirst(r_new[nb].nbrs .==i ) == 0 )
                              push!(rm_list,ind)
                              #println("$nb is nbr of $i, but not the opposite");
                        else
                            ind += 1;
                        end
                end


                for rm in rm_list
                    splice!(r_new[i].nbrs,rm)
                    splice!(r_new[i].angles,rm)
                end
                #if !isempty(rm_list) println("new nbrs list :", r_new[i].nbrs) end
        end

        # --------------- calc rho ---------------------------
        for i = 1:length(r_new)
              sqr_max = 0
              for nb in r_new[i].nbrs
                     dist2 = sqr_dist(r_new,i,nb)
                     if sqr_max < dist2
                        sqr_max= dist2;
                     end
              end
              if sqr_max==0 sqr_max=SQR_NEAR_NBR end
              r_new[i].rho = (length(r_new[i].nbrs)+1)/(sqr_max*pi);
        end

end

function find_min_distance(particles::Array{particle_t,1},i::Int)
    mind = (NEAR_NBR+1)^2;
    for nb in particles[i].nbrs
          dist2 = sqr_dist(particles,i,nb)
          if dist2 < mind
             mind = dist2;
          end
    end
    return sqrt(mind)
end

function divide_cell(particles::Array{particle_t,1}, particles_prev::Array{particle_t,1}, idx::Int , time)
        p = deepcopy(particles[idx]);
        teta = rand()*2*pi;
        particles[idx].x = particles[idx].x + (NEW_CELL_DIST/2)*cos(teta)
        particles[idx].y = particles[idx].y + (NEW_CELL_DIST/2)*sin(teta)
        if particles[idx].y > BOARD_BORDER_Y && PERIODIC_Y
              particles[idx].y =  particles[idx].y - BOARD_BORDER_Y
        elseif particles[idx].y < 0 && PERIODIC_Y
              particles[idx].y =  particles[idx].y + BOARD_BORDER_Y
        end

        if particles[idx].x > BOARD_BORDER_X && PERIODIC_X
              particles[idx].x=  particles[idx].x - BOARD_BORDER_X
        elseif particles[idx].x < 0 && PERIODIC_X
              particles[idx].x =  particles[idx].x + BOARD_BORDER_X
        end


        p.x = p.x - 2.5*cos(teta);
        p.y = p.y - 2.5*sin(teta);
        if p.y > BOARD_BORDER_Y && PERIODIC_Y
              p.y =  p.y - BOARD_BORDER_Y
        elseif particles[idx].y < 0 && PERIODIC_Y
              p.y =  p.y + BOARD_BORDER_Y
        end

        if p.x > BOARD_BORDER_X && PERIODIC_X
              p.x =  p.x - BOARD_BORDER_X
        elseif particles[idx].x < 0 && PERIODIC_X
              p.x =  p.x + BOARD_BORDER_X
        end

        push!(particles,p);
        push!(particles_prev,p); #just to enlarge the array

        # find new nbrs
        idx_new = length(particles)
        nbrs_list = [ particles[idx].nbrs ; idx ]
        particles[idx_new].nbrs = Int[]
        particles[idx_new].angles = Float64[]

        push!( particles_prev[idx].nbrs ,idx_new )
        for nb in particles_prev[idx].nbrs
                push!(particles_prev[nb].nbrs,idx_new);
                push!(particles_prev[nb].angles,10*pi); # just to keep the angles array in the same length
        end
        particles_prev[idx_new].nbrs = [ particles_prev[idx].nbrs ; idx ]
end


function delete_particles!(particles::Array{particle_t,1},new_p::Array{particle_t,1},
                           cl1::Array{Int,1},
                           cl2::Array{Int,1},
                           group_to_delete::Array{Int,1})
        if DELETE_DEATACHED && !isempty(group_to_delete)
            list = Array(Int,length(particles));
            list[:] = 1:length(particles);
            sorted_group = sort(group_to_delete, rev=true)
            prev_j = 0;

            print("deleting group:",group_to_delete')

            #deletting
            for j in sorted_group
                if prev_j != j # prevent double occurences
                      splice!(particles,j);
                      splice!(new_p,j);
                      splice!(list,j); #for future reindexing
                end
                prev_j = j;
            end

            #reindexing of nbrs
            for k =1:length(particles)
                reindex!(particles[k].nbrs,particles[k].angles,list)
            end

            #reindexing of clines
            reindex!(cl1,list);
            reindex!(cl2,list);
        end
end

function reindex!(arr::Array{Int,1}, ll::Array{Int,1}) #used for cline
        j=1
        while j <= length(arr)
            n = findfirst(ll,arr[j])
            if n == 0 #arr[j] is from deleted group
                error("Reindex: arr[$j]=$(arr[j]) is from deleted group")
            else
                arr[j] = n;
                j+=1
            end
        end
end

function reindex!(arr::Array{Int,1}, angles_arr::Array{Float64,1}, ll::Array{Int,1})
        j=1
        while j <= length(arr)
            n = findfirst(ll,arr[j])
            arr[j] = n;
            if n == 0 #arr[j] is from deleted group
                splice!(arr,j)
                splice!(angles_arr,j)
            else
                j+=1
            end
        end
end

