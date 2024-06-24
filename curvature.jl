function calc_curvature(particles, cl)
# Input: particles  - array of particle struct (among others, contain x and y fields)
#        cl  - contour line particles indices
# -----------------------------------------------------------------------------------

        #check
        if length(cl) < 4
                return (zeros(length(cl),2) , zeros(length(cl)) , zeros(length(cl)) )
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

        # smoothing
        Hmod_orig = H_abs.*H_sign;
        Hmod_orig[Hmod_orig .> F_H_MAX] = F_H_MAX;
        Hmod_orig[Hmod_orig .< -F_H_MAX] = -F_H_MAX;
        Hmod = weighted_gaussian_smooth(H_abs.*H_sign,ds);
        H = [Hmod Hmod] .* [H_sign H_sign] .* H_normal

        # second derivative of H and its smoothing
        dd_Hmod_orig  = ds_double_derivative(Hmod, ds);
        dd_Hmod = weighted_gaussian_smooth(dd_Hmod_orig[:],ds);
        ddH = [H_sign H_sign].* [ dd_Hmod dd_Hmod ] .* H_normal ;

        nanidx = findfirst(isnan(Hmod));
        if nanidx!=0 pr("H_mod[$nanidx] is NaN ! cl[$nanidx] = ",cl[nanidx]) end
        nanidx = findfirst(isnan(dd_Hmod));
        if nanidx!=0 pr("dd_Hmod[$nanidx] is NaN ! cl[$nanidx] = ",cl[nanidx]) end
        return (H, Hmod, Hmod_orig , ddH, dd_Hmod, dd_Hmod_orig , t_vec, ds)
end

function ds_double_derivative( arr::Array{Float64}, ds::Array{Float64,1} )
    first_der  = [ diff(arr) ,  arr[1,:]-arr[end,:] ]
    for s = 1:size(arr,2) ; first_der[:,s] = first_der[:,s]./ds ; end

    ds_avg = 0.5*( [ ds[end]; ds[1:end-1] ] + ds );
    second_der = [ first_der[1,:]-first_der[end,:] ; diff(first_der)];

    for s = 1:size(arr,2)
           second_der[:,s] = second_der[:,s]./ ds_avg
    end
    return second_der
end


function weighted_gaussian_smooth(orig_y::Array{Float64,1},ds::Array{Float64,1})
        len = length(orig_y);
        if len!= length(ds) error("y and ds vectors must be the same length") end
        orig_s = Array(Float64,len);
        orig_s[1] = ds[1];
        for i=2:len
             orig_s[i] = orig_s[i-1] + ds[i];
        end
        L = orig_s[end];
        ret_y = Array(Float64,len);
        for i=1:len
             ids = abs(orig_s-orig_s[i]);
             L_minus_ids = L - ids;
             ids_idx = ids .> L_minus_ids;
             ids[ids_idx] = L_minus_ids[ids_idx]
             gaussian = my_Gaussain(ids,0,SMOOTH_GAUSSIAN_WIDTH)
             ret_y[i] = sum(gaussian.*orig_y)/sum(gaussian);
        end

        return ret_y
end

function my_Gaussain ( x::Array{Float64,1}, mu::Number, sigma::Number)
    norm = 1/sigma/sqrt(2*pi);
    y = norm*exp (-((x - mu).^2)/(2*sigma^2));
    return y;
end
