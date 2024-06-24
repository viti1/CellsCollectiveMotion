#declare all files... (make the names global variables)

function open_pos_files()
        x_file      = open("$(directory)x.txt"      , "w") #1
        y_file      = open("$(directory)y.txt"      , "w") #2
        Vx_file     = open("$(directory)Vx.txt"    , "w")  #3
        Vy_file     = open("$(directory)Vy.txt"    , "w")  #4
        eta_x_file  = open("$(directory)etax.txt"  , "w") #5
        eta_y_file  = open("$(directory)etay.txt"  , "w") #6
        rho_file    = open("$(directory)rho.txt"  , "w")   #7

        nbr_file    = open("$(directory)nbr.txt"    , "w") #8

        return      [ x_file,  y_file,
                      Vx_file, Vy_file,
                      eta_x_file, eta_y_file,
                      rho_file,
                      nbr_file ]
end

function open_cline_files()
        cline_r_file  = open("$(directory)cline_r.txt"  , "w")
        if TWO_SIDED
            cline_l_file  = open("$(directory)cline_l.txt"  , "w")
            return [cline_r_file, cline_l_file]
        else
            return [ cline_r_file ];
        end

end

function open_force_files()
        if !PRINT_FORCE
              Hmod_file_r = open("$(directory)Hmod_r.txt"    , "w")
              if TWO_SIDED
                     Hmod_file_l = open("$(directory)Hmod_l.txt"    , "w")
                     return ( IOStream[], [Hmod_file_r] , [Hmod_file_l] )
              else
                     return ( IOStream[], [Hmod_file_r], IOStream[]);
              end
        end

        #  == all forces==
        Ffric_file_x    = open("$(directory)Ffricx.txt"    , "w")  #1
        Ffric_file_y    = open("$(directory)Ffricy.txt"    , "w")  #2
        Fvicsek_file_x  = open("$(directory)Fvicsekx.txt"  , "w")  #3
        Fvicsek_file_y  = open("$(directory)Fvicseky.txt"  , "w")  #4
        FU_file_x       = open("$(directory)FUx.txt"       , "w")  #5
        FU_file_y       = open("$(directory)FUy.txt"       , "w")  #6
        Fnoise_file_x   = open("$(directory)Fnoisex.txt"   , "w")  #7
        Fnoise_file_y   = open("$(directory)Fnoisey.txt"   , "w")  #8

        frc_files_bulk = [ Ffric_file_x , Ffric_file_y,      # 1,2
                          Fvicsek_file_x, Fvicsek_file_y,    # 3,4
                          FU_file_x, FU_file_y,             # 5,6
                          Fnoise_file_x, Fnoise_file_y     # 7,8
                        ]

        # == cline forces - right side ==


        Fresx_file_r = open("$(directory)Fresx_r.txt"    , "w")
        Fresy_file_r = open("$(directory)Fresy_r.txt"    , "w")

        Fcellx_file_r = open("$(directory)Fcellx_r.txt"    , "w")
        Fcelly_file_r = open("$(directory)Fcelly_r.txt"    , "w")

        Fcordx_file_r = open("$(directory)Fcordx_r.txt"    , "w")
        Fcordy_file_r = open("$(directory)Fcordy_r.txt"    , "w")

        Hx_file_r = open("$(directory)Hx_r.txt"    , "w")
        Hy_file_r = open("$(directory)Hy_r.txt"    , "w")
        Hmod_file_r = open("$(directory)Hmod_r.txt"    , "w")
        ddHx_file_r = open("$(directory)ddHx_r.txt"    , "w")
        ddHy_file_r = open("$(directory)ddHy_r.txt"    , "w")
        ds_file_r = open("$(directory)ds_r.txt"    , "w")

           #  -for debugging -
        ddHmod_file_r = open("$(directory)ddH_mod_r.txt"    , "w") #
        Hmod_orig_r = open("$(directory)Hmod_orig_r.txt"    , "w") #
        ddHmod_orig_r = open("$(directory)ddHmod_orig_r.txt"    , "w") #
           # -----

        frc_files_r = [
                Fresx_file_r,  Fresy_file_r,  # 1,2
                Fcellx_file_r, Fcelly_file_r, # 3,4
                Fcordx_file_r, Fcordy_file_r, # 5,6

                Hx_file_r, Hy_file_r,   # 7,8
                Hmod_file_r,            # 9
                ddHx_file_r,   ddHy_file_r,   # 10,11

                ds_file_r, #12

                ddHmod_file_r, # 13 - for debugging
                Hmod_orig_r,   # 14 - for debugging
                ddHmod_orig_r  # 15 - for debugging
                ]

        # == cline forces - left side ==
       if TWO_SIDED
          Fresx_file_l = open("$(directory)Fresx_l.txt"    , "w")
          Fresy_file_l = open("$(directory)Fresy_l.txt"    , "w")

          Fcellx_file_l = open("$(directory)Fcellx_l.txt"    , "w")
          Fcelly_file_l = open("$(directory)Fcelly_l.txt"    , "w")

          Fcordx_file_l = open("$(directory)Fcordx_l.txt"    , "w")
          Fcordy_file_l = open("$(directory)Fcordy_l.txt"    , "w")

          Hx_file_l = open("$(directory)Hx_l.txt"    , "w")
          Hy_file_l = open("$(directory)Hy_l.txt"    , "w")
          Hmod_file_l = open("$(directory)Hmod_l.txt"    , "w")
          ddHx_file_l = open("$(directory)ddHx_l.txt"    , "w")
          ddHy_file_l = open("$(directory)ddHy_l.txt"    , "w")
          ds_file_l = open("$(directory)ds_l.txt"    , "w")

            # - for debugging -
          ddHmod_file_l = open("$(directory)ddH_mod_l.txt"    , "w") #
          Hmod_orig_l = open("$(directory)Hmod_orig_l.txt"    , "w") #
          ddHmod_orig_l = open("$(directory)ddHmod_orig_l.txt"    , "w") #
            # ----

          frc_files_l = [
                  Fresx_file_l,  Fresy_file_l,   # 1,2
                  Fcellx_file_l, Fcelly_file_l,  # 3,4
                  Fcordx_file_l, Fcordy_file_l,  # 5,6

                  Hx_file_l, Hy_file_l,   # 7,8
                  Hmod_file_l,            # 9
                  ddHx_file_l,   ddHy_file_l,   # 10,11

                  ds_file_l, #12

                  ddHmod_file_l, # 13 - for debugging
                  Hmod_orig_l,   # 14 - for debugging
                  ddHmod_orig_l  # 15 - for debugging
          ]

           return ( frc_files_bulk,  frc_files_r , frc_files_l )
       else
           return ( frc_files_bulk,  frc_files_r, IOStream[])
       end
end

function open_all_files()
        return ( open_force_files() , open_pos_files(), open_cline_files() )
end

function close_all_files()
        all_files = [frc_files...;pos_files;cl_files]

        for fl in all_files
              close(fl)
        end
end
