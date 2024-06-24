function inargs(input_param::ASCIIString)
        if length(ARGS) > 1
           for s in ARGS[2:end]
                if beginswith(s,input_param)
                   if ( s[1:search(s,'=')-1] == input_param || s == input_param )
                        return true
                   end
                end
           end
        end
        return false
end

function extract_param(input_param::ASCIIString)
        for s in ARGS[2:end]
             if s[1:search(s,'=')-1] == input_param

                    if contains(s,"[")
                        argument = s[search(s,'=')+2:end-1];
                        ret=float(split(argument,","))
                    else

                        ret = float(s[search(s,'=')+1:end])
                    end
                    println("input parameter: $input_param = ",ret)
                    return ret
             end
        end
        error("parameter $input_param cannot be extracted")
end

function extract_bool_param(input_param::ASCIIString)
        for s in ARGS[2:end]
             if s[1:search(s,'=')-1] == input_param
                    val = s[search(s,'=')+1:end]
                    println("input parameter: $input_param = ",val)
                    if val == "true"
                        return true
                    elseif val == "false"
                        return false
                    else
                        error("Wrong Parameter $input_param=$val ." )
                    end
             end
        end
        error("parameter $input_param cannot be extracted")
end

function extract_string_param(input_param::ASCIIString)
        for s in ARGS[2:end]
             if s[1:search(s,'=')-1] == input_param
                    val = s[search(s,'=')+1:end]
                    println("input parameter: $input_param = ",val)
                    return val
             end
        end
        error("parameter $input_param cannot be extracted")
end



function checkinargsiflong(islong_original::Bool)
        if length(ARGS) > 1
           for s in ARGS[2:end]
                if s=="long"
                   return true
                elseif s=="short"
                   return false
                end
           end
        end
        return islong_original
end

#BETA = checkinargs(BETA,"beta");
#( inargs , BETA_param ) = checkargs(ARGS,"beta")
#if BETA_inargs; BETA=BETA_param end


function beta2sigma(beta::Number)
        beta = float(beta)
        btosig =
        [ 0   0.157
          10  0.45
          20  0.6
          30  0.75
          60  1
          90  1.35 #1.4 was bit too large
          130 1.8
        ]

        if beta < 0
                error("Beta (=$beta) must be positive")
        elseif beta > maximum(btosig[:,1])
                error("Beta (=$beta) mut be less than ",maximum(btosig[:,1]))
        end

        ind = findfirst(btosig[:,1].>beta) - 1  ;
        d_beta = beta - btosig[ind,1];
        slope  = ( btosig[ind+1,2] - btosig[ind,2] )/( btosig[ind+1,1] - btosig[ind,1] )
        return btosig[ind,2] + d_beta*slope;
end

function max2kappa(new_slope::Number)
    base_slope = 6* 61.42/0.01;
    base_k = 3.1e7

    return base_k*new_slope/base_slope
end
