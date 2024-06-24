include("definitions.jl")
#if ~isdefined(:DEF_INCLUDED) require("definitions.jl") end
include("files.jl")
include("Printings.jl")
include("GenerateArrays.jl")
include("movement.jl")
include("coast_line.jl")
#  -------------------   MAIN -------------------------------------------------------------
        tic() ; iter=1;
        cline_r = Int[]; cline_l = Int[]; to_delete2 = Int[];
        (frc_files,pos_files,cl_files) = open_all_files();

        if iter%PRINT_STEP==0 && PRINT_EVERY_I; print("----- Iteration 1 -") end
        ( r_cur, r_new ) = initialize_arrays(start_from_last,directory);

        if iter%PRINT_STEP==0; print_particles_array(r_cur, P1 ,pos_files); end
        print_last(r_cur);

        (cline_r, to_delete1)  = find_coast_line(true,r_cur,iter,cl_files[1],cline_r);
        if TWO_SIDED; (cline_l, to_delete2) = find_coast_line(false, r_cur,iter,cl_files[2],cline_l); end

        iter=2
        while ( iter <= NUM_OF_CYCLES )
            if iter%PRINT_STEP==0 && PRINT_EVERY_I pr("----- Iteration $(iter) ($(iter/PRINT_STEP))---- DIR = $(directory[11:end]) ") end

            delete_particles!(r_cur,r_new,cline_r,cline_l,[to_delete1; to_delete2]) #detected in find_coast_line
            cells_movement!(r_cur, r_new,
                             cline_r, cline_l,
                             frc_files..., iter)
                             # particles might be added (proliferate) after movement (in update_map )
            if iter%PRINT_STEP==0 print_particles_array(r_new, P1 , pos_files) end
            print_last(r_new); # to have the last iteration if something goes wrong

             ( cline_r , to_delete1 ) = find_coast_line(true, r_new,iter,cl_files[1],cline_r);
             if TWO_SIDED; ( cline_l, to_delete2) = find_coast_line(false, r_new,iter,cl_files[2], cline_l); end

            iter+=1;
            ( r_cur , r_new ) = (r_new , r_cur )
        end

        close_all_files();
	toc()
