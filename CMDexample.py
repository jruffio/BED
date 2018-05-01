__author__ = 'jruffio'

import numpy as np
from scipy.stats import norm
import os

outputDir = os.path.join("/home","sda","jruffio","BED","figures")

if 1:
    import astropy.io.fits as pyfits
    hdulist = pyfits.open("/home/sda/jruffio/NIRC2/besancon-p_Pup_0.230sq_deg.fits")
    table = hdulist[0].data
    hdulist.close()

    table = table[:,np.where(np.isfinite(table[7,:]))[0]]

    print(table.shape)
    N_stars = table.shape[1]
    Besancon_FOV = 0.23 # deg^2
    star_density = N_stars/Besancon_FOV

    mystar_Hmag = 5
    mystar_F444Wmag = 5
    mystar_H = 10**(-mystar_Hmag/2.5)
    mystar_F444W = 10**(-mystar_F444Wmag/2.5)

    myblob_Hcontmag = 12.5+0.8 +mystar_Hmag
    myblob_Hcont_sig = 10**(-myblob_Hcontmag/2.5)/5.
    myblob_Hcont = myblob_Hcont_sig
    myblob_F444Wcontmag = 12.5 +mystar_F444Wmag
    myblob_F444Wcont = 10**(-myblob_F444Wcontmag/2.5)
    myblob_F444Wcont_sig = 10**(-myblob_F444Wcontmag/2.5)/5.

    myblob2_Hcontmag = 15+0.8 +mystar_Hmag
    myblob2_Hcont_sig = 10**(-myblob2_Hcontmag/2.5)/5.
    myblob2_Hcont = 5*myblob2_Hcont_sig
    myblob2_F444Wcontmag = 15 +mystar_F444Wmag
    myblob2_F444Wcont = 10**(-myblob2_F444Wcontmag/2.5)
    myblob2_F444Wcont_sig = 10**(-myblob2_F444Wcontmag/2.5)/5.

    # plt.errorbar(np.log10(myblob_Hcont/myblob_F444Wcont))


    #0/ Besancon_mass
    #1/ Besancon_logg
    #2/ Besancon_teff
    #3/ Besancon_mh
    #4/ Besancon_rad_fixed
    #5/ Besancon_d
    #6/ Besancon_av
    #7/ 2MASS H magnitude
    #8/ 2MASS K magnitude
    #9/ JWST F444W magnitude
    print(table.shape)

    import matplotlib.pyplot as plt
    K_vec = 10**(-table[8,:]/2.5)
    F444W_vec =  10**(-table[9,:]/2.5)
    H_vec = 10**(-table[7,:]/2.5)
    Kmag_vec = table[8,:]
    Hmag_vec = table[7,:]
    F444Wmag_vec = table[9,:]
    mass_vec = table[0,:]


    # test color plot
    if 0:
        # plt.figure(2)
        # plt.hist(mass_vec,100)
        # plt.show()

        plt.figure(1,figsize=(6,6))
        Ns = 100
        color_bounds = [0,1.5]
        # mag_bounds = [5,30]
        mag_bounds = [0,5]
        xedges,yedges = np.linspace(color_bounds[0],color_bounds[1],Ns),np.linspace(mag_bounds[0],mag_bounds[1],Ns)
        H,xedges,yedges = np.histogram2d(Hmag_vec-F444Wmag_vec,mass_vec,bins=[xedges,yedges])
        x_centers = [(x1+x2)/2. for x1,x2 in zip(xedges[0:Ns-1],xedges[1:Ns])]
        y_centers = [(y1+y2)/2. for y1,y2 in zip(yedges[0:Ns-1],yedges[1:Ns])]
        ravel_H = np.ravel(H)
        ind = np.argsort(ravel_H)
        cum_ravel_H = np.zeros(np.shape(ravel_H))
        cum_ravel_H[ind] = np.cumsum(ravel_H[ind])
        cum_H = 1-np.reshape(cum_ravel_H/np.nanmax(cum_ravel_H),np.shape(H))
        image = np.clip(np.log10(H),0,np.inf).T
        image[np.where(cum_H.T>0.9973)] = np.nan
        where_star = np.where(mass_vec> 1)
        where_bd = np.where((mass_vec> 0)*(mass_vec< 1))
        print("N stars",np.size(Hmag_vec[where_star]))
        print("N bd",np.size(Hmag_vec[where_bd]))
        plt.scatter(Hmag_vec[where_star]-F444Wmag_vec[where_star],F444Wmag_vec[where_star],s=1,zorder=5,c="black",alpha=0.1)#,c=mass_vec
        plt.scatter(Hmag_vec[where_bd]-F444Wmag_vec[where_bd],F444Wmag_vec[where_bd],s=1,zorder=5,c="red",alpha=0.1)#,c=mass_vec
        plt.imshow(image,origin ="lower",
                   extent=[color_bounds[0],color_bounds[1],mag_bounds[0],mag_bounds[1]],
                   aspect="auto",zorder=10)#,alpha = 0.5,interpolation="spline16",
        plt.xlim(color_bounds)
        plt.ylim(mag_bounds)
        levels = [0.6827,0.9545,0.9973]
        xx,yy = np.meshgrid(x_centers,y_centers)
        CS = plt.contour(xx,yy,cum_H.T,levels = levels,linestyles="-",linewidths=[1],colors=("black",),zorder=15)
        plt.gca().invert_yaxis()
        plt.xlabel("H-F444W",fontsize=15)
        plt.ylabel("F444W",fontsize=15)
        plt.show()
        exit()

    # New color plot
    if 1:
        plt.figure(1,figsize=(7,7))
        Ns = 100
        color_bounds = [-0.5,2]# [0,1.5]
        mag_bounds = [5,30]
        xedges,yedges = np.linspace(color_bounds[0],color_bounds[1],Ns),np.linspace(mag_bounds[0],mag_bounds[1],Ns)

        Ns_hd = 1000
        xedges_hd,yedges_hd = np.linspace(color_bounds[0],color_bounds[1],Ns_hd),np.linspace(mag_bounds[0],mag_bounds[1],Ns_hd)
        dy = yedges_hd[1]-yedges_hd[0]
        dx = xedges_hd[1]-xedges_hd[0]
        x_centers_hd = [(x1+x2)/2. for x1,x2 in zip(xedges_hd[0:Ns_hd-1],xedges_hd[1:Ns_hd])]
        y_centers_hd = [(y1+y2)/2. for y1,y2 in zip(yedges_hd[0:Ns_hd-1],yedges_hd[1:Ns_hd])]
        xx_hd,yy_hd = np.meshgrid(x_centers_hd,y_centers_hd)

        where_star = np.where(mass_vec> 0.65)
        where_bd = np.where((mass_vec> 0)*(mass_vec< 0.65))
        print("N stars",np.size(Hmag_vec[where_star]))
        print("N bd",np.size(Hmag_vec[where_bd]))
        H_star_hd,xedges_hd,yedges_hd = np.histogram2d(Hmag_vec[where_star]-F444Wmag_vec[where_star],F444Wmag_vec[where_star],bins=[xedges_hd,yedges_hd])
        pdf_stars_hd = H_star_hd.T/(np.sum(H_star_hd.T)*dx*dy)
        H_bd_hd,xedges_hd,yedges_hd = np.histogram2d(Hmag_vec[where_bd]-F444Wmag_vec[where_bd],F444Wmag_vec[where_bd],bins=[xedges_hd,yedges_hd])
        pdf_bd_hd = H_bd_hd.T/(np.sum(H_bd_hd.T)*dx*dy)
        for where_obj,cmap,color in zip([where_star,where_bd],["viridis","hot"],["purple","red"]):
            H,xedges,yedges = np.histogram2d(Hmag_vec[where_obj]-F444Wmag_vec[where_obj],F444Wmag_vec[where_obj],bins=[xedges,yedges])
            x_centers = [(x1+x2)/2. for x1,x2 in zip(xedges[0:Ns-1],xedges[1:Ns])]
            y_centers = [(y1+y2)/2. for y1,y2 in zip(yedges[0:Ns-1],yedges[1:Ns])]
            ravel_H = np.ravel(H)
            ind = np.argsort(ravel_H)
            cum_ravel_H = np.zeros(np.shape(ravel_H))
            cum_ravel_H[ind] = np.cumsum(ravel_H[ind])
            cum_H = 1-np.reshape(cum_ravel_H/np.nanmax(cum_ravel_H),np.shape(H))
            image = np.clip(np.log10(H),0,np.inf).T
            image[np.where(cum_H.T>0.9545)] = np.nan
            # plt.scatter(Hmag_vec[where_obj]-F444Wmag_vec[where_obj],F444Wmag_vec[where_obj],s=2,zorder=5,c="black")#,c=mass_vec
            plt.imshow(image,origin ="lower",
                       extent=[color_bounds[0],color_bounds[1],mag_bounds[0],mag_bounds[1]],
                       aspect="auto",zorder=10,cmap=cmap,alpha = 0.5)#,alpha = 0.5,interpolation="spline16",,alpha = 0.75
            plt.xlim(color_bounds)
            plt.ylim(mag_bounds)
            levels = [0.6827]
            xx,yy = np.meshgrid(x_centers,y_centers)
            CS = plt.contour(xx,yy,cum_H.T,levels = levels,linestyles="-",linewidths=[2],colors=("white",),zorder=15)
            levels = [0.9545,0.9973]
            CS = plt.contour(xx,yy,cum_H.T,levels = levels,linestyles="-",linewidths=[2],colors=(color,),zorder=15)

        resel_area = np.pi*(3.5*14.166)**2 # mas^2
        Besancon_FOV = 0.23*(3.6e6)**2/resel_area # mas^2 ## 0.23 deg^2
        N_bd = np.size(Hmag_vec[where_bd])
        bd_density = N_bd/Besancon_FOV
        N_star = np.size(Hmag_vec[where_star])
        star_density = N_star/Besancon_FOV
        null_density = 1-(bd_density+star_density)

        mystar_Hmag = 5
        mystar_F444Wmag = 5


        Hcontmag_list = np.array([12.5,12.5,12.5]) +mystar_Hmag
        F444Wcontmag_list = np.array([10,12.5,15]) + mystar_F444Wmag
        Hcont_snr_list = np.array([30,3,-1])
        F444Wcont_snr_list = np.array([5,5,5])
        colors=["#ff9900","#006699","grey"]
        for k,(myblob_Hcontmag,myblob_F444Wcontmag,myblob_H_snr,myblob_F444W_snr,color) in enumerate(zip(Hcontmag_list,F444Wcontmag_list,Hcont_snr_list,F444Wcont_snr_list,colors)):
            print("candidate",k)
            myblob_Hcont_sig = 10**(-myblob_Hcontmag/2.5)/5.
            myblob_Hcont = myblob_Hcont_sig*myblob_H_snr
            print(myblob_Hcont_sig,myblob_Hcont,-2.5*np.log10(myblob_Hcont_sig),-2.5*np.log10(myblob_Hcont))
            myblob_F444Wcont_sig = 10**(-myblob_F444Wcontmag/2.5)/5.
            myblob_F444Wcont = myblob_F444Wcont_sig*myblob_F444W_snr
            levels = [0.99865,0.9545,0.84135]
            levels = np.array([norm.pdf(norm.ppf(l,0,1)) for l in levels])
            levels = levels*norm.pdf(np.zeros(3),loc=0,scale=myblob_F444Wcont_sig)*\
                     norm.pdf(np.zeros(3),loc=0,scale=myblob_Hcont_sig)/norm.pdf(np.zeros(3),0,1)
            like2d = norm.pdf(10**(-yy/2.5),loc=myblob_F444Wcont,scale=myblob_F444Wcont_sig)*\
                norm.pdf((10**(-xx/2.5))*(10**(-yy/2.5)),loc=myblob_Hcont,scale=myblob_Hcont_sig)
            like2d_hd = norm.pdf(10**(-yy_hd/2.5),loc=myblob_F444Wcont,scale=myblob_F444Wcont_sig)*\
                norm.pdf((10**(-xx_hd/2.5))*(10**(-yy_hd/2.5)),loc=myblob_Hcont,scale=myblob_Hcont_sig)
            image =np.log10(like2d)
            image[np.where(like2d<levels[0])] = np.nan
            # plt.imshow(-image,origin ="lower",
            #            extent=[color_bounds[0],color_bounds[1],mag_bounds[0],mag_bounds[1]],
            #            aspect="auto",zorder=10,cmap="hot")#,alpha = 0.5,interpolation="spline16",
            CS = plt.contour(xx_hd,yy_hd,like2d_hd,levels = levels,linestyles=[":","--","-"],linewidths=[2,2,2],colors=(color,),zorder=15)
            # CS = plt.contour(xx,yy,like2d,levels = [levels[0]],linestyles=":",linewidths=[1],colors=("black",),zorder=15)
            # CS = plt.contour(xx,yy,like2d,levels = [levels[1]],linestyles="--",linewidths=[1],colors=("black",),zorder=15)
            # CS = plt.contour(xx,yy,like2d,levels = [levels[2]],linestyles="-",linewidths=[1],colors=("black",),zorder=15)

            # plt.figure(3)
            # plt.subplot(1,3,1)
            # plt.imshow(like2d_hd)
            # plt.colorbar()
            # plt.subplot(1,3,2)
            # plt.imshow(H_star_hd.T)
            # plt.colorbar()
            # plt.subplot(1,3,3)
            # plt.imshow(H_bd_hd.T)
            # plt.colorbar()
            # plt.show()
            proba_star = np.sum(like2d_hd*pdf_stars_hd)*dx*dy*star_density
            proba_bd = np.sum(like2d_hd*pdf_bd_hd)*dx*dy*bd_density
            proba_null = norm.pdf(myblob_F444Wcont,loc=0,scale=myblob_F444Wcont_sig)*norm.pdf(myblob_Hcont,loc=0,scale=myblob_Hcont_sig)*null_density
            normalization = proba_star+proba_bd+proba_null
            print(star_density,bd_density,null_density)
            # print(proba_star,proba_bd,proba_null)
            print(proba_star/normalization,proba_bd/normalization,proba_null/normalization,proba_bd/proba_star)

        plt.gca().invert_yaxis()
        plt.gca().tick_params(axis='x', labelsize=20)
        plt.gca().tick_params(axis='y', labelsize=20)
        plt.xlabel("H-F444W",fontsize=20)
        plt.ylabel("F444W",fontsize=20)

        plt.savefig(os.path.join(outputDir,"starproba.pdf"),bbox_inches='tight')
        plt.savefig(os.path.join(outputDir,"starproba.png"),bbox_inches='tight')
        plt.show()


    # Old color plot
    if 0:
        plt.figure(1,figsize=(12,12))

        plt.subplot(2,2,1)
        I_F_vec = np.arange(-5*myblob_F444Wcont_sig+myblob_F444Wcont,5*myblob_F444Wcont_sig+myblob_F444Wcont,0.01*myblob_F444Wcont_sig)
        I_vec = np.logspace(-12,0,1000)
        I_log_vec = np.log10(I_vec)# np.linspace(-12,0,1000)
        R_vec = np.logspace(-3,3,num=1000)
        R_log_vec = np.linspace(-3,3,num=1000)
        integrand_H = norm.pdf(I_vec,loc=myblob_Hcont,scale=myblob_Hcont_sig)[:,None]*norm.pdf((I_vec[:,None]/R_vec[None,:]),loc=myblob_F444Wcont,scale=myblob_F444Wcont_sig)
        posterior = np.sum(integrand_H,axis=0)
        posterior=posterior/np.nanmax(posterior)
        # plt.plot(R_vec,posterior)
        plt.plot(R_log_vec,posterior)
        from  scipy.interpolate import interp1d
        f = interp1d(posterior,R_log_vec,bounds_error=False, fill_value=np.nan)
        R_log_upperlimit = f(0.01)
        # plt.plot([R_log_upperlimit,R_log_upperlimit],[0,0.01])
        plt.scatter([R_log_upperlimit],[0.01],s=50)
        plt.xlim([-1.,0])
        plt.ylim([0,1])
        plt.xlabel("Contrast ratio H/F444W")
        plt.ylabel("Posterior")
        # plt.gca().set_xscale("log")




        plt.subplot(2,2,3)
        Ns = 100
        xedges,yedges = np.linspace(-1,0,Ns),np.linspace(-12,0,Ns)
        H,xedges,yedges = np.histogram2d(np.log10(H_vec/F444W_vec*mystar_F444W/mystar_H),np.log10(F444W_vec/mystar_F444W),bins=[xedges,yedges])
        x_centers = [(x1+x2)/2. for x1,x2 in zip(xedges[0:Ns-1],xedges[1:Ns])]
        y_centers = [(y1+y2)/2. for y1,y2 in zip(yedges[0:Ns-1],yedges[1:Ns])]
        ravel_H = np.ravel(H)
        ind = np.argsort(ravel_H)
        cum_ravel_H = np.zeros(np.shape(ravel_H))
        cum_ravel_H[ind] = np.cumsum(ravel_H[ind])
        cum_H = 1-np.reshape(cum_ravel_H/np.nanmax(cum_ravel_H),np.shape(H))
        image = np.clip(np.log10(H),0,np.inf).T
        image[np.where(cum_H.T>0.9973)] = np.nan
        plt.scatter(np.log10(H_vec/F444W_vec*mystar_F444W/mystar_H),np.log10(F444W_vec/mystar_F444W),s=2,zorder=5,c="black")#,c=mass_vec
        plt.imshow(image,origin ="lower",extent=[-1,0,-12,0],aspect="auto",zorder=10)#,alpha = 0.5,interpolation="spline16",
        plt.xlim([-1,0])
        plt.ylim([-12,0])
        levels = [0.6827,0.9545,0.9973]
        xx,yy = np.meshgrid(x_centers,y_centers)
        CS = plt.contour(xx,yy,cum_H.T,levels = levels,linestyles="-",linewidths=[1],colors=("black",),zorder=15)
        plt.xlabel("Contrast ratio H/F444W")
        plt.ylabel("Planet-to-star flux (F444W Log contrast)")
        plt.errorbar(R_log_upperlimit,np.log10(myblob_F444Wcont),
                     xerr=[[10],[0]],
                     yerr=[[np.log10(myblob_F444Wcont)-np.log10(myblob_F444Wcont-myblob_F444Wcont_sig)],[np.log10(myblob_F444Wcont+myblob_F444Wcont_sig)-np.log10(myblob_F444Wcont)]],
                     linewidth=2,zorder=20,ecolor="red")
        ax0 = plt.gca()

        plt.subplot(2,2,2)
        levels = [0.99865,0.9545,0.84135]
        levels = np.array([norm.pdf(norm.ppf(l,0,1)) for l in levels])
        levels = levels*norm.pdf(np.zeros(3),loc=0,scale=myblob_F444Wcont_sig)*\
                 norm.pdf(np.zeros(3),loc=0,scale=myblob_Hcont_sig)/norm.pdf(np.zeros(3),0,1)
        like2d = norm.pdf(10**yy,loc=myblob_F444Wcont,scale=myblob_F444Wcont_sig)*\
            norm.pdf((10**xx)*(10**yy),loc=myblob_Hcont,scale=myblob_Hcont_sig)
        # print(like2d)
        # ravel_like2d = np.ravel(like2d)
        # ind = np.argsort(ravel_like2d)
        # cum_ravel_like2d = np.zeros(np.shape(ravel_like2d))
        # cum_ravel_like2d[ind] = np.cumsum(ravel_like2d[ind])
        # cum_like2d = 1-np.reshape(cum_ravel_like2d/np.nanmax(cum_ravel_like2d),np.shape(like2d))
        image =np.log10(like2d)
        image[np.where(like2d<levels[0])] = np.nan
        plt.imshow(-image,origin ="lower",extent=[-1,0,-12,0],aspect="auto",zorder=10,cmap="hot")#,alpha = 0.5,interpolation="spline16",
        plt.xlim([-1,0])
        # plt.xlim([-12,0])
        plt.ylim([-12,0])
        xx,yy = np.meshgrid(x_centers,y_centers)
        CS = plt.contour(xx,yy,like2d,levels = levels,linestyles="-",linewidths=[1],colors=("black",),zorder=15)
        plt.xlabel("Contrast ratio H/F444W")
        plt.ylabel("Planet-to-star flux (F444W Log contrast)")
        # from matplotlib.patches import Ellipse
        # myellipse = Ellipse(xy=[-1,np.log10(myblob_F444Wcont)],width=R_upperlimit+1,height=)
        # plt.gca().add_artist(myellipse)

        plt.sca(ax0)
        plt.imshow(-image,origin ="lower",extent=[-1,0,-12,0],aspect="auto",zorder=11,cmap="hot")#,alpha = 0.5,interpolation="spline16",

        plt.subplot(2,2,4)
        I_vec = np.logspace(-12,0,1000)
        I_log_vec = np.log10(I_vec)# np.linspace(-12,0,1000)
        dIlog = I_log_vec[1]-I_log_vec[0]
        # print(myblob_F444Wcont,myblob_F444Wcont_sig)
        like=norm.pdf(I_vec,loc=myblob_F444Wcont,scale=myblob_F444Wcont_sig)
        plt.plot(I_log_vec,like/np.nanmax(like),label="F444W")
        like=norm.pdf(I_vec,loc=myblob_Hcont,scale=myblob_Hcont_sig)
        plt.plot(I_log_vec,like/np.nanmax(like),label="H")
        plt.legend()
        plt.xlabel("Planet-to-star flux (Log contrast)")
        plt.ylabel("PDF")
        # # plt.plot(I_vec,norm.cdf(I_vec,loc=myblob_F444Wcont,scale=myblob_F444Wcont_sig))
        # post = np.cumsum(like)/(np.max(np.cumsum(like)))
        # plt.plot( I_log_vec,post)
        # plt.gca().set_xscale("log")

        plt.tight_layout()
        plt.show()





        exit()
        plt.scatter(H_vec/F444W_vec*mystar_F444W/mystar_H,F444W_vec/mystar_F444W,s=2,alpha=0.1,c=mass_vec)
        # plt.plot([1e-12,1e-0],[1e-12,1e-0])
        # plt.plot(K_vec,".")
        # plt.plot(F444W_vec/K_vec,".")
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.xlim([-0.1,1.4])
        plt.ylim([1e-12,1e-0])
        plt.figure(2)
        plt.scatter(Hmag_vec-F444Wmag_vec,F444Wmag_vec,s=2,alpha=0.1,c=mass_vec)
        plt.colorbar()
        plt.show()
        exit()


        # levels = [0.5,0.158655,0.0227501,0.00134990,0.0000316712,2.86652e-7,9.86588e-10] #0,1,2,3,4,5,6 sigmas
        # CS = plt.contour(xx,yy,tailposterior_map,levels = [levels[2]],linestyles="-.",linewidths=[3],colors=("white",))
        # fmt = {}
        # strs = [r"$5\sigma$",r"$4\sigma$",r"$3\sigma$"]#[r"$10^{-3}$",r"$10^{-5}$"]
        # for l,s in zip(CS.levels,strs):
        #     fmt[l] = s
        # plt.clabel(CS,inline=True,fontsize=15,fmt=fmt)#manual=[(30,-4),(50,-4),(70,-4)]
        # # plt.clabel(CS,inline=1,fontsize=10,fmt='%1.2f')#manual=[(0.5,interp1d(sep_cont_samples,contrast_5sig)(0.5))]