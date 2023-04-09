
#--------------------------------------------------#
#                     ODE function                 #
#--------------------------------------------------#

function func_2D!(dU,U,p,t)

    Dc,DuS,DvS,DuB,DvB,χ,γ,η,dx,dy,Nx,Ny,param_reac_S,param_reac_B,αS,βS,αB,αB = p 


    ############################ SURFACE ############################

    for i in 1:Nx
        r_uS, r_vS = reaction_S(U[i,1,1],U[i,1,2],param_reac_S)
        ip1 = min(i+1,Nx)
        im1 = max(i-1,1)
        ΔuS = (U[ip1,1,1] + U[im1,1,1] - 2*U[i,1,1])/dx^2
        ΔvS = (U[ip1,1,2] + U[im1,1,2] - 2*U[i,1,2])/dx^2
        dU[i,1,1] = DuS * ΔuS + r_uS + η*αS*(U[i,Ny+1,1]-U[i,1,1])
        dU[i,1,2] = DvS * ΔvS + r_vS + η*βS*(U[i,Ny+1,2]-U[i,1,2])
    end

    #################################################################


    ############################ BULK ###############################

    @inbounds for j in 2:(Ny+1)
        @simd for i in 1:Nx

            imax = min(i,Nx)
            ip1 = min(i+1,Nx)
            ip2 = min(i+2,Nx)
            im1 = max(i-1,1)
            im2 = max(i-2,1)
            jmax = min(j,Ny+1)
            jp1 = min(j+1,Ny+1)
            jp2 = min(j+2,Ny+1)
            jm1 = max(j-1,2)
            jm2 = max(j-2,2)

            ########### Flux ###############################
            
            ### X ###
            
            dx_uB_im = (taxis(U[i,j,1],U[i,j,2],param_reac_B) - taxis(U[im1,j,1],U[im1,j,1],param_reac_B))/dx
            # if dx_uB_im>0
            #     ## Third order
            #     c_im = 1.0/6.0 * (-U[im2,j,3] + 5*U[im1,j,3] + 2*U[imax,j,3]) 
            # else
            #     ## Third order
            #     c_im = 1.0/6.0 * (2*U[im1,j,3] + 5*U[imax,j,3] - U[ip1,j,3])
            # end
            c_im = IfElse.ifelse(dx_uB_im>0.0,1.0/6.0 * (-U[im2,j,3] + 5*U[im1,j,3] + 2*U[imax,j,3]),1.0/6.0 * (2*U[im1,j,3] + 5*U[imax,j,3] - U[ip1,j,3]))
            flux_x_im = dx_uB_im * taxis_h(c_im,param_reac_B)
            
            dx_uB_ip = (taxis(U[ip1,j,1],U[ip1,j,2],param_reac_B) - taxis(U[i,j,1],U[i,j,2],param_reac_B))/dx
            # if dx_uB_ip>0
            #     ## Third order
            #     c_ip = 1.0/6.0 * (-U[im1,j,3] + 5*U[imax,j,3] + 2*U[ip1,j,3]) 
            # else
            #     ## Third order
            #     c_ip = 1.0/6.0 * (2*U[imax,j,3] + 5*U[ip1,j,3] - U[ip2,j,3])
            # end
            c_ip = IfElse.ifelse( dx_uB_ip>0.0,1.0/6.0 * (-U[im1,j,3] + 5*U[imax,j,3] + 2*U[ip1,j,3]),1.0/6.0 * (2*U[imax,j,3] + 5*U[ip1,j,3] - U[ip2,j,3]))
            flux_x_ip = dx_uB_ip * taxis_h(c_ip,param_reac_B)

            ### Y ###
            
            dy_uB_jm = (taxis(U[i,j,1],U[i,j,2],param_reac_B) - taxis(U[i,jm1,1],U[i,jm1,2],param_reac_B))/dy
            # if dy_uB_jm>0
            #     ## Third order
            #     c_jm = 1.0/6.0 * (-U[i,jm2,3] + 5*U[i,jm1,3] + 2*U[i,jmax,3]) 
            # else
            #     ## Third order
            #     c_jm = 1.0/6.0 * (2*U[i,jm1,3] + 5*U[i,jmax,3] - U[i,jp1,3])
            # end
            c_jm = IfElse.ifelse(dy_uB_jm>0.0,1.0/6.0 * (-U[i,jm2,3] + 5*U[i,jm1,3] + 2*U[i,jmax,3]),1.0/6.0 * (2*U[i,jm1,3] + 5*U[i,jmax,3] - U[i,jp1,3]))
            flux_y_jm = dy_uB_jm * taxis_h(c_jm,param_reac_B)
            
            dy_uB_jp = (taxis(U[i,jp1,1],U[i,jp1,2],param_reac_B) - taxis(U[i,j,1],U[i,j,2],param_reac_B))/dy
            # if dy_uB_jp>0
            #     ## Third order
            #     c_jp = 1.0/6.0 * (-U[i,jm1,3] + 5*U[i,jmax,3] + 2*U[i,jp1,3]) 
            # else
            #     ## Third order
            #     c_jp = 1.0/6.0 * (2*U[i,jmax,3] + 5*U[i,jp1,3] - U[i,jp2,3])
            # end
            c_jp = IfElse.ifelse(dy_uB_jp>0.0,1.0/6.0 * (-U[i,jm1,3] + 5*U[i,jmax,3] + 2*U[i,jp1,3]),1.0/6.0 * (2*U[i,jmax,3] + 5*U[i,jp1,3] - U[i,jp2,3]))
            flux_y_jp = dy_uB_jp * taxis_h(c_jp,param_reac_B)
            ################################################

            ########## Laplace #############################
            Δc = (U[ip1,j,3] + U[im1,j,3] - 2*U[i,j,3])/dx^2 + (U[i,jp1,3] + U[i,jm1,3] - 2*U[i,j,3])/dy^2

            b_u_BS = U[i,Ny+1,1] + dy*η*αB/DuB .* (U[i,1,1] - U[i,Ny+1,1])
            uB_jp = j==(Ny+1) ? b_u_BS : U[i,j+1,1]
            ΔuB = (U[ip1,j,1] + U[im1,j,1] - 2*U[i,j,1])/dx^2 + (uB_jp + U[i,jm1,1] - 2*U[i,j,1])/dy^2

            b_v_BS = U[i,Ny+1,2] + dy*η*βB/DvB .* (U[i,1,2] - U[i,Ny+1,2])
            vB_jp = j==(Ny+1) ? b_v_BS : U[i,j+1,2]
            ΔvB = (U[ip1,j,2] + U[im1,j,2] - 2*U[i,j,2])/dx^2 + (vB_jp + U[i,jm1,2] - 2*U[i,j,2])/dy^2
            ################################################

            g = growth(U[i,j,3],param_reac_B)
            chemo = (flux_x_im - flux_x_ip)/dx + (flux_y_jm - flux_y_jp)/dy
            r_u,r_v = reaction_B(U[i,j,1],U[i,j,2],U[i,j,3],param_reac_B)
            dU[i,j,1] = DuB * ΔuB + r_u
            dU[i,j,2] = DvB * ΔvB + r_v
            dU[i,j,3] = Dc*Δc + χ*chemo + γ*g
        end
    end

end


function func_1D!(dU,U,p,t)

    DuS,DvS,DuB,DvB,Dc,χ,γ,η,dx,Nx,param_reac_S,param_reac_B,αS,βS,αB,βB = p 


    ############################ SURFACE ############################

    for i in 1:Nx
        ip1 = min(i+1,Nx)
        im1 = max(i-1,1)
        r_uS, r_vS = reaction_S(U[i,1],U[i,2],param_reac_S)
        ΔuS = (U[ip1,1] + U[im1,1] - 2*U[i,1])/dx^2
        ΔvS = (U[ip1,2] + U[im1,2] - 2*U[i,2])/dx^2
        dU[i,1] = DuS * ΔuS + r_uS + η*αS*(U[i,3]-U[i,1])
        dU[i,2] = DvS * ΔvS + r_vS + η*βS*(U[i,4]-U[i,2])
    end

    #################################################################


    ############################ BULK ###############################


    for i in 1:Nx
        imax = min(i,Nx)
        ip1 = min(i+1,Nx)
        ip2 = min(i+2,Nx)
        im1 = max(i-1,1)
        im2 = max(i-2,1)

        ########### Flux ###############################
            
        dx_uB_im = (taxis(U[i,3],U[i,4],param_reac_B) - taxis(U[im1,3],U[im1,4],param_reac_B))/dx

        c_im = IfElse.ifelse(dx_uB_im>0.0,1.0/6.0 * (-U[im2,5] + 5*U[im1,5] + 2*U[imax,5]),1.0/6.0 * (2*U[im1,5] + 5*U[imax,5] - U[ip1,5]))
        flux_x_im = dx_uB_im * taxis_h(c_im,param_reac_B)
        
        dx_uB_ip = (taxis(U[ip1,3],U[ip1,4],param_reac_B) - taxis(U[i,3],U[i,4],param_reac_B))/dx

        c_ip = IfElse.ifelse( dx_uB_ip>0.0,1.0/6.0 * (-U[im1,5] + 5*U[imax,5] + 2*U[ip1,5]),1.0/6.0 * (2*U[imax,5] + 5*U[ip1,5] - U[ip2,5]))
        flux_x_ip = dx_uB_ip * taxis_h(c_ip,param_reac_B)

        Δc = (U[ip1,5] + U[im1,5] - 2*U[i,5])/dx^2 
        g = growth(U[i,5],param_reac_B)
        chemo = (flux_x_im - flux_x_ip)/dx

        r_uB, r_vB = reaction_B(U[i,3],U[i,4],U[i,5],param_reac_B)
        ΔuB = (U[ip1,3] + U[im1,3] - 2*U[i,3])/dx^2
        ΔvB = (U[ip1,4] + U[im1,4] - 2*U[i,4])/dx^2
        dU[i,3] = DuB * ΔuB + r_uB + η*αB*(U[i,1]-U[i,3])
        dU[i,4] = DvB * ΔvB + r_vB + η*βB*(U[i,2]-U[i,4])
        dU[i,5] = Dc*Δc + χ*chemo + γ*g
    end

end
