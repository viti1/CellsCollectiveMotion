function smooth(H,ds)
        if WEIGHTED_SMOOTH
            #inhomogenic_smooth(H,ds,0.05)
            weighted_gaussian_smooth(H,ds)
        else
            filter_smooth(H,0.4)
        end
end

function inhomogenic_smooth(H::Array{Float64},ds::Array{Float64,1},fraction::Float64)
    if length(H)!=length(ds) error("length(H)!=length(ds)") end
    S_eq = linspace(0.0,sum(ds),int(ceil(sum(ds)/minimum(ds)))) #equal distances
    S_orig = Array(Float64,length(ds)+1) #original distances
    S_orig[1]=0.0;
    for i=2:length(ds)+1
         S_orig[i] = S_orig[i-1]+ds[i-1]
    end
    H_eq = convert_points(H,S_orig,S_eq)
    smoothed_eq = filter_smooth(H_eq,fraction);
    ret = convert_points(smoothed_eq,S_eq,S_orig)
end

function convert_points(y_orig, x_orig, x_new)
        if length(y_orig) != length(x_orig)-1 pr("Warning smooth: length(y_orig)($(length(y_orig))) != length(x_orig)-1 ($(length(x_orig)-1))") end
        y_new = Array(Float64,length(x_new)-1);
        m_orig = diff([y_orig;y_orig[1]]);

        for k = 1:length(x_new)-1
             io = findfirst( x_orig .> x_new[k] ) - 1
             y_new[k] = y_orig[io] + m_orig[io]*(x_new[k] - x_orig[io])
        end
        return y_new
end

function filter_smooth(arr::Array{Float64},fraction::Float64)
        arr_size = size(arr);
        len = length(arr);
        half_len = floor(len/2)+1;
        arr = [ arr[half_len:end] ; arr[1:end] ; arr[1:half_len]];
        start_idx =  len - half_len + 2 ;
        sizef = length(arr)
        filter = sqrt_filter(sizef,0,floor(sizef*fraction/2))   ;
        ret = ifft(fft(arr).*filter);
        ret = ret[ len-half_len+2 : 2*len-half_len+1 ];
        ret = real(ret) ;
end

function gauss_smooth(arr::Array{Float64})
        halfWidth = 15.0 ;
        gaussFilter = my_Gaussain(-halfWidth:halfWidth, 0.0, 2 ) ;

        ret=Array(Float64,size(arr))
        for k = 1:size(arr,2)
            smoothedVector = conv([flipud(arr[:,k]) ; arr[2:end,k]], gaussFilter) ;
            ret[:,k] = smoothedVector[(end+1)/2:end-halfWidth];
        end
        return ret
end

function my_Gaussain( window, mu::Number, sigma::Number)
    y = exp(-((window - mu).^2)/(2*sigma^2)) ;
    y = y / sum(y);
    return y;
end

function my_Gaussain( x::Array{Float64,1}, mu::Number, sigma::Number)
    norm = 1/sigma/sqrt(2*pi);
    y = norm*exp(-((x - mu).^2)/(2*sigma^2));
    return y;
end


function sqrt_filter(FiltSize,f_low, f_high )
        if f_low > f_high
            error("Error: f_high must be larger than f_low");
        elseif f_low < 0 || f_high < 0
            error("Error: negative frequency is illegal");
        end

        u_low  =  round(f_low  ) + 1 ;
        u_high =  round(f_high ) + 1 ;
        u_max  =  floor(FiltSize/2) + 1 ;

        if u_high > u_max
            println("f_high is out of boundaries of this image frequencies, setting it to max frequency $(u_max/du_inv)");
            u_high = u_max;
        end

        filt = zeros(u_max);
        filt[u_low:u_high] = ones(u_high - u_low + 1);

        if isodd(FiltSize) # odd FiltSize
            filt = [filt , flipud(filt[2:end])];
        else # even FiltSize
            filt = [filt , flipud(filt[2:end-1])];
        end

        if length(filt) !=  FiltSize
            error("Error generating Filter");
        end

        return filt
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




# ---------------------- old stuff --------------------------------------------
#function find_five_points(A,clsize,i,ds)

#    prev = (i==1) ? clsize : i-1
#    prev_prev = ( i < 3 ) ? clsize+i-2 : i-2
#    next = (i==clsize) ? 1 : i+1
#    next_next = (i > clsize-2) ? 2-(clsize-i) : i+2

#    ret = Array(Float64,5,2)
#    ret[1,:] = A[prev_prev,:]
#    ret[2,:] = A[prev,:]
#    ret[3,:] = A[i,:]
#    ret[4,:] = A[next,:]
#    ret[5,:] = A[next_next,:]

#    for k=1:4
#        if ret[k+1,2]-ret[k,2] .< - HALF_BOARD_BORDER_Y
#              ret[k+1,2] = ret[k+1,2] + BOARD_BORDER_Y
#        elseif ret[k+1,2]-ret[k,2] > HALF_BOARD_BORDER_Y
#              ret[k+1,2] = ret[k+1,2] - BOARD_BORDER_Y
#        end
#    end

#    s = Array(Float64,5,1);
#    s[1] = -ds[prev]-ds[prev_prev]
#    s[2] = -ds[prev]
#    s[3] = 0
#    s[4] = ds[i]
#    s[5] = ds[next]+ds[i]
#    return (ret,s)
#end

#function spline_curv(curv,ds)
#        clsize = size(curv,1)
#        second_der = Array(Float64,size(curv))
#        forth_der = Array(Float64,size(curv))

#        for i=1:clsize
#                #construct a matrix
#                (y,x) = find_five_points(curv,clsize,i,ds);
#                mat = [ ones(5,1) x x.^2 x.^3 x.^4]
#                a_coef = inv(mat) * y ;
#                second_der[i,:] = 2*a_coef[3,:]
#                forth_der[i,:] = 24*a_coef[5,:]
#        end
#        return ( second_der , forth_der );
#end

#function smooth!(H, H_modulus,np)
#       TH = 0.009
#       cl_size = length(H_modulus);
#       h=0
#       flag = true
#       while flag && h<30
#           idx = find( abs(diff( [ H[:,1]; H[1,1] ] )).> TH  )
#           H_new = ([H[2:end,:];H[1,:]] + H[:,:]+ [ H[end,:]; H[1:end-1,:] ] )/3
#           H[idx,:] = H_new[idx,:]
#           flag = findfirst(abs(diff( [ H[:,1]; H[1,1] ] )).> TH ) != 0 ;
#           h+=1
#       end

#       H_sign = ( (H[:,1].*np[:,2]-H[:,2].*np[:,1]) .> 0 )*2 - 1 ;
#       H_modulus[:] = sqrt(H[:,1].^2+H[:,2].^2).*H_sign ;
#end

#function spline(H,ds)
#    clsize = size(H,1)
#    ddH = Array(Float64,size(H))
#    for i=1:clsize
#            #construct a matrix
#            (y,x) = find_five_points(H,clsize,i);
#            mat = [ ones(5,1) x.^2 x.^3 x.^4]
#            a_coef = inv(mat) * y ;
#            ddH[i,:] = 2*a_coef[3,:]
#    end
#    return ddH;
#end

