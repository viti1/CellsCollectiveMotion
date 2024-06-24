
# printing to files
function pr(str...)
   println(str...)
end

function prarr(io::IOStream,arr::Array)
   println(io,join(arr," "))
end

function prarr(arr::Array)
   println(join(arr," "))
end

function arrtostr(Arr::Array)
    return join(Arr," ")
end

function arrtostrn(Arr::Array)
     return  join(Arr," ")*"\n";
end

function print_particles_array(particles::Array{particle_t,1} ,idx::Int, files )
        for k = 1:length(particles) #DEBUG
            @printf(files[1] , "%lf\t", particles[k].x)
            @printf(files[2] , "%lf\t", particles[k].y)
            @printf(files[3] , "%lf\t", particles[k].v[1])
            @printf(files[4] , "%lf\t", particles[k].v[2])
            @printf(files[5] , "%lf\t", particles[k].eta[1])
            @printf(files[6] , "%lf\t", particles[k].eta[2])
            @printf(files[7], "%lf\t", particles[k].rho)
        end
        for j=1:7
            pr(files[j]);
        end


        if idx > length(particles)
                println("Warning: Printing neighbors of p=%idx is not possible since num of particls is ",length(particles))
        elseif idx!=0
                print(files[8],[idx particles[idx].nbrs'])
        end

end

# print the particles array as 7 files
function print_last( particles::Array{particle_t,1} )
         pos = Array(Float64,length(particles),7)

         for k = 1:length(particles)
            pos[k,1] = particles[k].x
            pos[k,2] = particles[k].y
            pos[k,3] = particles[k].v[1]
            pos[k,4] = particles[k].v[2]
            pos[k,5] = particles[k].eta[1]
            pos[k,6] = particles[k].eta[2]
            pos[k,7] = particles[k].rho
         end

         writedlm("$(directory)poscur.txt",pos');
end

function print_cl_forces(files::Array{IOStream,1},
                Fres::Array{Float64,2},
                Fcell::Array{Float64,2},
                Fcord::Array{Float64,2},
                H::Array{Float64,2},
                Hmod::Array{Float64,1} ,
                ddH::Array{Float64,2},
                ds::Array{Float64,1},
                ddHmod::Array{Float64} ,
                Hmod_orig::Array{Float64} ,
                ddHmod_orig::Array{Float64} ,
                )
        if !PRINT_FORCE
                println(files[1],join(Hmod',"\t"))
                return
        end
        println(files[1],join(Fres[1,:],"\t"))
        println(files[2],join(Fres[2,:],"\t"))
        println(files[3],join(Fcell[1,:],"\t "))
        println(files[4],join(Fcell[2,:],"\t"))
        println(files[5],join(Fcord[1,:],"\t"))
        println(files[6],join(Fcord[2,:],"\t"))
        println(files[7],join(H[:,1]',"\t"))
        println(files[8],join(H[:,2]',"\t"))
        println(files[9],join(Hmod',"\t"))
        println(files[10],join(ddH[:,1]',"\t"))
        println(files[11],join(ddH[:,2]',"\t"))
        println(files[12],join(ds',"\t"))

        println(files[13],join(ddHmod',"\t"))
        println(files[14],join(Hmod_orig',"\t"))
        println(files[15],join(ddHmod_orig',"\t"))
end

function print_cl_forces(files::Array{IOStream,1},
                Fres::Array{Float64,2},
                Fcell::Array{Float64,2},
                Fcord::Array{Float64,2},
                H::Array{Float64,2},
                Hmod::Array{Float64,1} ,
                ddH::Array{Float64,2},
                )
    if !PRINT_FORCE
        println(files[1],join(Hmod',"\t"))
        return
    end

    print(files[1],arrtostrn(Fres[1,:]))
    print(files[2],arrtostrn(Fres[2,:]))
    print(files[3],arrtostrn(Fcell[1,:]))
    print(files[4],arrtostrn(Fcell[2,:]))
    print(files[5],arrtostrn(Fcord[1,:]))
    print(files[6],arrtostrn(Fcord[2,:]))
    print(files[7],arrtostrn(H[:,1]'))
    print(files[8],arrtostrn(H[:,2]'))
    print(files[9],arrtostrn(Hmod'))
    print(files[10],arrtostrn(ddH[:,1]'))
    print(files[11],arrtostrn(ddH[:,2]'))
end

function print_cl_forces_of_particle( idx::Int, cline,
                                Fres, Fcell, Fcord,
                                side_string::ASCIIString)
        idx_cl = findfirst(cline,idx);
        if idx_cl!=0
            pr(" ",side_string," side cline ")
            println("\tFres:   " , Fres[:,idx_cl]')
            println("\tFcell:  " , Fcell[:,idx_cl]')
            println("\tFcords: " , Fcord[:,idx_cl]')
        end
end

function print_bulk_forces_of_particle( files::Array{IOStream,1},
                            f_friction::Array{Float64} ,
                            f_vicsek::Array{Float64} ,
                            f_interaction::Array{Float64},
                            f_noise::Array{Float64},
							f_wall::Array{Float64}	)						
    if !PRINT_FORCE; return; end
    print(files[1], f_friction[1]    ,'\t')
    print(files[2], f_friction[2]    ,'\t')
    print(files[3], f_vicsek[1]      ,'\t')
    print(files[4], f_vicsek[2]      ,'\t')
    print(files[5], f_interaction[1] ,'\t')
    print(files[6], f_interaction[2] ,'\t')
    print(files[7], f_noise[1]       ,'\t')
    print(files[8], f_noise[2]       ,'\t')
	print(files[9], f_wall[1]       ,'\t')
    print(files[10] , f_wall[2]       ,'\t')

end

function print_bulk_forces_of_particle(  f_friction::Array{Float64} ,
                            f_vicsek::Array{Float64},
                            f_interaction::Array{Float64},
                            f_noise::Array{Float64},
							f_wall::Array{Float64}
                         )

        println("~~~ The Forces are ~~    ");
        println("\tFriction   : ", f_friction' )
        println("\tVicsek     : " , f_vicsek')
        println("\tInteraction: ", f_interaction')
        println("\tNoise      : ", f_noise')
		println("\tWall      : ", f_wall')
end

function print_new_line_to_files(files::Array{IOStream,1})
    for fl in files
            pr(fl)
    end
end

