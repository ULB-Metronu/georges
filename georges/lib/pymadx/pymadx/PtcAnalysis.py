"""
Analysis utilities for PTC data.
"""

from . import Ptc as _Ptc
from .Data import Tfs as _Tfs
import numpy as _np
import csv
import time

class PtcAnalysis(object):
    """
    Deprecated.

    Optical function calculation for PTC tracking data.

    This has be reimplemented and replaced by C++ implementation
    in rebdsim.

    """
    def __init__(self,ptcInput = None, ptcOutput = None) : 

        # Load input rays 
        if type(ptcInput) == str: 
            self.ptcInput = _Ptc.LoadInrays(ptcInput)
        else : 
            self.ptcInput = ptcInput 

        # Load output rays 
        if type(ptcOutput) == str:
            self.ptcOutput = _Tfs(ptcOutput)
        else : 
            self.ptcOutput = ptcOutput
    
    def SamplerLoop(self):
        #for segment in ptcOutput.segments:
        for isampler in (0,self.ptcOutput.isegment,1):
            samplerData = self.ptcOutput.GetSegment(isampler)
            
            xrms  = samplerData.GetColumn('X').std()
            yrms  = samplerData.GetColumn('Y').std()
            pxrms = samplerData.GetColumn('PX').std()
            pyrms = samplerData.GetColumn('PY').std()

            print(isampler, xrms, pxrms, yrms, pyrms)

    def CalculateOpticalFunctions(self, output):
        """
        Calulates optical functions from a PTC output file
    
        output - the name of the output file
        """
        sampler  = []
        S        = []
        betx     = []
        bety     = []
        alphx    = []
        alphy    = []
        dispx    = []
        dispy    = []
        dispxp   = []
        dispyp   = []
        emittx   = []
        emitty   = []
        sigmax   = []
        sigmay   = []
        sigmaxp  = []
        sigmayp  = []
        sigmaxxp = []
        sigmayyp = []
        meanx    = []
        meany    = []
        meanxp   = []
        meanyp   = []
        Wgt      = []
        sigm_emx = []
        sigm_emy = []
        sigm_btx = []
        sigm_bty = []
        sigm_alx = []
        sigm_aly = []
        

        
        #for segment in ptcOutput.segments:
        for isampler in range(0,self.ptcOutput.isegment,1):
            #samplerData = self.ptcOutput.GetSegment(segment)
            samplerData = self.ptcOutput.GetSegment(isampler)
            print('segment:', (isampler+1) ,'/', self.ptcOutput.isegment)
        
            x  = samplerData.GetColumn('X')
            y  = samplerData.GetColumn('Y') 
            xp = samplerData.GetColumn('PX')
            yp = samplerData.GetColumn('PY')
            s  = samplerData.GetColumn('S')
            E0 = samplerData.GetColumn('E')  #this is the specified beam  energy
            PT = samplerData.GetColumn('PT') # this is defined as pt= deltaE/(p0*c)
            E=E0*(1+PT) #This is the energy with a spread

            #Calculate sums
            s_s    = _np.sum(s)
            x_s    = _np.sum(x)
            y_s    = _np.sum(y)
            xp_s   = _np.sum(xp)
            yp_s   = _np.sum(yp)
            E_s    = _np.sum(E)
            EE_s   = _np.sum(E*E)
            xx_s   = _np.sum(x*x)
            xxp_s  = _np.sum(x*xp)
            xpxp_s = _np.sum(xp*xp)
            xpE_s  = _np.sum(xp*E)
            xE_s   = _np.sum(x*E)
            yy_s   = _np.sum(y*y)
            yyp_s  = _np.sum(y*yp)
            ypyp_s = _np.sum(yp*yp)
            ypE_s  = _np.sum(yp*E)
            yE_s   = _np.sum(y*E)

            wgt   = len(x) #this is a valid for ptc output as weight is always 1

            #normalise
            s_s    /= wgt
            x_s    /= wgt
            y_s    /= wgt
            xp_s   /= wgt
            yp_s   /= wgt
            E_s    /= wgt
            EE_s   /= wgt
            xx_s   /= wgt
            xxp_s  /= wgt
            xpxp_s /= wgt
            xpE_s  /= wgt
            xE_s   /= wgt
            yy_s   /= wgt
            yyp_s  /= wgt
            ypyp_s /= wgt
            ypE_s  /= wgt
            yE_s   /= wgt
            
            x_sv    = _np.sum(x-_np.mean(x))
            xx_sv   = _np.sum((x-_np.mean(x))*(x-_np.mean(x)))
            y_sv    = _np.sum(y-_np.mean(y))
            yy_sv   = _np.sum((y-_np.mean(y))*(y-_np.mean(y)))

            xp_sv   = _np.sum(xp-_np.mean(xp))
            xpxp_sv = _np.sum((xp-_np.mean(xp))*(xp-_np.mean(xp)))
            yp_sv   = _np.sum(yp-_np.mean(yp))
            ypyp_sv = _np.sum((yp-_np.mean(yp))*(yp-_np.mean(yp)))
            
            #Calculate variances/sigmas
            variance_x = (xx_sv - (x_sv * x_sv)/wgt)/wgt
            variance_y = (yy_sv - (y_sv * y_sv)/wgt)/wgt
            sigma_x    = _np.sqrt(variance_x)
            sigma_y    = _np.sqrt(variance_y)
            variance_xp = (xpxp_sv - (xp_sv * xp_sv)/wgt)/wgt
            variance_yp = (ypyp_sv - (yp_sv * yp_sv)/wgt)/wgt
            sigma_xp   = _np.sqrt(variance_xp)
            sigma_yp   = _np.sqrt(variance_yp)
        
            #Calculate means
            x_m  = x_s/wgt
            y_m  = y_s/wgt
            xp_m = xp_s/wgt
            yp_m = yp_s/wgt
        
            #Calculate sigmas
            sigma_x_xp = xxp_s - xx_s * xp_s;
            sigma_y_yp = yyp_s - y_s * yp_s;
            
            #####################
            #Calculate moments for error calculations
            
            mx_1_1  = _np.sum((x-x_m)*(xp-xp_m))/wgt
            mx_0_2  = _np.sum((xp-xp_m)**2)/wgt
            mx_2_0  = _np.sum((x-x_m)**2)/wgt
            mx_2_2  = _np.sum((x-x_m)**2*(xp-xp_m)**2)/wgt
            mx_1_3  = _np.sum((x-x_m)*(xp-xp_m)**3)/wgt
            mx_3_1  = _np.sum((x-x_m)**3*(xp-xp_m))/wgt
            mx_4_0  = _np.sum((x-x_m)**4)/wgt
            mx_0_4  = _np.sum((xp-xp_m)**4)/wgt

            my_1_1  = _np.sum((y-y_m)*(yp-yp_m))/wgt
            my_0_2  = _np.sum((yp-yp_m)**2)/wgt
            my_2_0  = _np.sum((y-y_m)**2)/wgt
            my_2_2  = _np.sum((y-y_m)**2*(yp-yp_m)**2)/wgt
            my_1_3  = _np.sum((y-y_m)*(yp-yp_m)**3)/wgt
            my_3_1  = _np.sum((y-y_m)**3*(yp-yp_m))/wgt
            my_4_0  = _np.sum((y-y_m)**4)/wgt
            my_0_4  = _np.sum((yp-yp_m)**4)/wgt
            
            
            #caculate higher moment covariances
            cov_vv_x       = -((-3+wgt)*mx_2_0**2)/((-1+wgt)*wgt)+mx_4_0/wgt
            cov_vv_xp      = -((-3+wgt)*mx_0_2**2)/((-1+wgt)*wgt)+mx_0_4/wgt
            cov_cc_xxp_xxp = -((-2+wgt)*mx_1_1**2)/((-1+wgt)*wgt)+(mx_0_2*mx_2_0)/((-1+wgt)*wgt)+mx_2_2/wgt
            cov_vc_x_xxp   = -((-3+wgt)*mx_1_1*mx_2_0)/((-1+wgt)*wgt)+mx_3_1/wgt
            cov_vc_xp_xxp  = -((-3+wgt)*mx_1_1*mx_0_2)/((-1+wgt)*wgt)+mx_1_3/wgt
            cov_vv_x_xp    = (2*mx_1_1**2)/((-1+wgt)*wgt)-(mx_0_2*mx_2_0)/wgt+mx_2_2/wgt

            cov_vv_y       = -((-3+wgt)*my_2_0**2)/((-1+wgt)*wgt)+my_4_0/wgt
            cov_vv_yp      = -((-3+wgt)*my_0_2**2)/((-1+wgt)*wgt)+my_0_4/wgt
            cov_cc_yyp_yyp = -((-2+wgt)*my_1_1**2)/((-1+wgt)*wgt)+(my_0_2*my_2_0)/((-1+wgt)*wgt)+my_2_2/wgt
            cov_vc_y_yyp   = -((-3+wgt)*my_1_1*my_2_0)/((-1+wgt)*wgt)+my_3_1/wgt
            cov_vc_yp_yyp  = -((-3+wgt)*my_1_1*my_0_2)/((-1+wgt)*wgt)+my_1_3/wgt
            cov_vv_y_yp    = (2*my_1_1**2)/((-1+wgt)*wgt)-(my_0_2*my_2_0)/wgt+my_2_2/wgt
 
            #Calculate the moments using the sums
            xx_s   -= x_s*x_s
            xxp_s  -= x_s*xp_s
            xpxp_s -= xp_s*xp_s
            yy_s   -= y_s*y_s  
            yyp_s  -= y_s*yp_s
            ypyp_s -= yp_s*yp_s

            #Calculate emittance
            emitt_x = _np.sqrt(xx_s*xpxp_s - xxp_s*xxp_s)
            emitt_y = _np.sqrt(yy_s*ypyp_s - yyp_s*yyp_s)

            #Calculate optical functions
            beta_x  =  xx_s  / emitt_x;
            beta_y  =  yy_s  / emitt_y;
            alph_x  = -xxp_s / emitt_x;
            alph_y  = -yyp_s / emitt_y;
            
            #if there is no energy spread it is expected that the dispersion would evaluate to inf, hence ignore
            #'divide by 0' runtime warning
            with _np.errstate(divide='ignore'):
                disp_x  = (xE_s  - (x_s  * E_s)) / (EE_s - (E_s * E_s))
                disp_xp = (xpE_s - (xp_s * E_s)) / (EE_s - (E_s * E_s))
                disp_y  = (yE_s  - (y_s  * E_s)) / (EE_s - (E_s * E_s))
                disp_yp = (ypE_s - (yp_s * E_s)) / (EE_s - (E_s * E_s))
            

            #error calculation
            d_emx_d_xx   = xpxp_s/(2*emitt_x)
            d_emx_d_xxp  = -xxp_s/emitt_x
            d_emx_d_xpxp = xx_s/(2*emitt_x)
            
            var_emitt_x = (d_emx_d_xx)**2*cov_vv_x+(d_emx_d_xxp)**2*cov_cc_xxp_xxp+(d_emx_d_xpxp)**2*cov_vv_xp
            var_emitt_x += 2*(d_emx_d_xx*d_emx_d_xxp)*cov_vc_x_xxp+2*(d_emx_d_xx*d_emx_d_xpxp)*cov_vv_x_xp
            var_emitt_x += 2*(d_emx_d_xxp*d_emx_d_xpxp)*cov_vc_xp_xxp
            
            sigma_emitt_x = _np.sqrt(var_emitt_x)

            d_emy_d_yy   = ypyp_s/(2*emitt_y)
            d_emy_d_yyp  = -yyp_s/emitt_y
            d_emy_d_ypyp = yy_s/(2*emitt_y)

            var_emitt_y = (d_emy_d_yy)**2*cov_vv_y+(d_emy_d_yyp)**2*cov_cc_yyp_yyp+(d_emy_d_ypyp)**2*cov_vv_yp
            var_emitt_y += 2*(d_emy_d_yy*d_emy_d_yyp)*cov_vc_y_yyp+2*(d_emy_d_yy*d_emy_d_ypyp)*cov_vv_y_yp
            var_emitt_y += 2*(d_emy_d_yyp*d_emy_d_ypyp)*cov_vc_yp_yyp
            
            sigma_emitt_y = _np.sqrt(var_emitt_y)
         
            d_btx_d_xx   = (xx_s*xpxp_s-2*xxp_s**2)/(2*(xx_s*xpxp_s - xxp_s*xxp_s)**(3./2.))
            d_btx_d_xxp  = (xx_s*xxp_s)/((xx_s*xpxp_s - xxp_s*xxp_s)**(3./2.))
            d_btx_d_xpxp = -(xx_s**2)/(2*(xx_s*xpxp_s - xxp_s*xxp_s)**(3./2.))

            var_beta_x = (d_btx_d_xx)**2*cov_vv_x+(d_btx_d_xxp)**2*cov_cc_xxp_xxp+(d_btx_d_xpxp)**2*cov_vv_xp
            var_beta_x+= 2*(d_btx_d_xx*d_btx_d_xxp)*cov_vc_x_xxp+2*(d_btx_d_xx*d_btx_d_xpxp)*cov_vv_x_xp
            var_beta_x+= 2*(d_btx_d_xxp*d_btx_d_xpxp)*cov_vc_xp_xxp
            
            sigma_beta_x = _np.sqrt(var_beta_x)

            d_bty_d_yy   = (yy_s*ypyp_s-2*yyp_s**2)/(2*(yy_s*ypyp_s - yyp_s*yyp_s)**(3./2.))
            d_bty_d_yyp  = (yy_s*yyp_s)/((yy_s*ypyp_s - yyp_s*yyp_s)**(3./2.))
            d_bty_d_ypyp = -(yy_s**2)/(2*(yy_s*ypyp_s - yyp_s*yyp_s)**(3./2.))

            var_beta_y = (d_bty_d_yy)**2*cov_vv_y+(d_bty_d_yyp)**2*cov_cc_yyp_yyp+(d_bty_d_ypyp)**2*cov_vv_yp
            var_beta_y+= 2*(d_bty_d_yy*d_bty_d_yyp)*cov_vc_y_yyp+2*(d_bty_d_yy*d_bty_d_ypyp)*cov_vv_y_yp
            var_beta_y+= 2*(d_bty_d_yyp*d_bty_d_ypyp)*cov_vc_yp_yyp
            
            sigma_beta_y = _np.sqrt(var_beta_y)

            d_alx_d_xx   = -(xxp_s*xpxp_s)/(2*(xx_s*xpxp_s - xxp_s*xxp_s)**(3./2.))
            d_alx_d_xxp  = (xx_s*xpxp_s)/((xx_s*xpxp_s - xxp_s*xxp_s)**(3./2.))
            d_alx_d_xpxp = -(xx_s*xxp_s)/(2*(xx_s*xpxp_s - xxp_s*xxp_s)**(3./2.))

            var_alph_x = (d_alx_d_xx)**2*cov_vv_x+(d_alx_d_xxp)**2*cov_cc_xxp_xxp+(d_alx_d_xpxp)**2*cov_vv_xp
            var_alph_x+= 2*(d_alx_d_xx*d_alx_d_xxp)*cov_vc_x_xxp+2*(d_alx_d_xx*d_alx_d_xpxp)*cov_vv_x_xp
            var_alph_x+= 2*(d_alx_d_xxp*d_alx_d_xpxp)*cov_vc_xp_xxp
            
            sigma_alph_x = _np.sqrt(var_alph_x)

            
            d_aly_d_yy   = -(yyp_s*ypyp_s)/(2*(yy_s*ypyp_s - yyp_s*yyp_s)**(3./2.))
            d_aly_d_yyp  = (yy_s*ypyp_s)/((yy_s*ypyp_s - yyp_s*yyp_s)**(3./2.))
            d_aly_d_ypyp = -(yy_s*yyp_s)/(2*(yy_s*ypyp_s - yyp_s*yyp_s)**(3./2.))

            var_alph_y = (d_aly_d_yy)**2*cov_vv_y+(d_aly_d_yyp)**2*cov_cc_yyp_yyp+(d_aly_d_ypyp)**2*cov_vv_yp
            var_alph_y+= 2*(d_aly_d_yy*d_aly_d_yyp)*cov_vc_y_yyp+2*(d_aly_d_yy*d_aly_d_ypyp)*cov_vv_y_yp
            var_alph_y+= 2*(d_aly_d_yyp*d_aly_d_ypyp)*cov_vc_yp_yyp
            
            sigma_alph_y = _np.sqrt(var_alph_y)
            
            """
            print "mx_2_0= ",mx_2_0,", mx_0_2= ",mx_0_2,", mx_1_1= ",mx_1_1
            print "mx_4_0= ",mx_4_0,", mx_0_4= ",mx_0_4
            print "mx_3_1= ",mx_3_1,", mx_1_3= ",mx_1_3,", mx_2_2= ",mx_2_2
            print "my_2_0= ",my_2_0,", my_0_2= ",my_0_2,", my_1_1= ",my_1_1
            print "my_4_0= ",my_4_0,", my_0_4= ",my_0_4
            print "my_3_1= ",my_3_1,", my_1_3= ",my_1_3,", my_2_2= ",my_2_2

            print "cov_vv_x= ", cov_vv_x,", cov_vv_xp= ",cov_vv_xp,", cov_cc_xxp_xxp= ",cov_cc_xxp_xxp
            print "cov_vc_x_xxp= ",cov_vc_x_xxp,", cov_vc_xp_xxp= ",cov_vc_xp_xxp,", cov_vv_x_xp=",cov_vv_x_xp
            print "cov_vv_y= ",cov_vv_y,", cov_vv_yp= ",cov_vv_yp,", cov_cc_yyp_yyp= ",cov_cc_yyp_yyp
            print "cov_vc_y_yyp= ",cov_vc_y_yyp,", cov_vc_yp_yyp= ",cov_vc_yp_yyp,", cov_vv_y_yp=",cov_vv_y_yp
           
            print "dedx= ",d_emx_d_xx," ,dedxp= ",d_emx_d_xpxp,", dedxxp= ",d_emx_d_xxp
            print "dbdx= ",d_btx_d_xx,",dbdxp= ",d_btx_d_xpxp,", dbdxxp= ",d_btx_d_xxp
            print "dadx= ",d_alx_d_xx,",dadxp= ",d_alx_d_xpxp,", dadxxp= ",d_alx_d_xxp
            print "dedy= ",d_emy_d_yy," ,dedyp= ",d_emy_d_ypyp,", dedyyp= ",d_emy_d_yyp
            print "dbdy= ",d_bty_d_yy,",dbdyp= ",d_bty_d_ypyp,", dbdyyp= ",d_bty_d_yyp
            print "dady= ",d_aly_d_yy,",dadyp= ",d_aly_d_ypyp,", dadyyp= ",d_aly_d_yyp

            print "var_emitt_x= ",var_emitt_x,", var_beta_x= ",var_beta_x,", var_alph_x= ",var_alph_x
            print "var_emitt_y= ",var_emitt_y,", var_beta_y= ",var_beta_y,", var_alph_y= ",var_alph_y

            print "Emiitance_x= ",emitt_x," Sigma_Emittance_x = ",sigma_emitt_x, " Frac error= ", sigma_emitt_x/emitt_x 
            print "Beta_x= ",beta_x," Sigma_Beta_x = ",sigma_beta_x, " Frac error= ", sigma_beta_x/beta_x
            print "Alpha_x= ",alph_x," Sigma_Alpha_x = ",sigma_alph_x, " Frac error= ", sigma_alph_x/alph_x
            """

            sampler.append(isampler)
            S.append(s_s)            
            betx.append(beta_x)
            bety.append(beta_y)
            alphx.append(alph_x)
            alphy.append(alph_y)
            dispx.append(disp_x)
            dispy.append(disp_y)
            dispxp.append(disp_xp)
            dispyp.append(disp_yp)
            emittx.append(emitt_x)
            emitty.append(emitt_y)
            sigmax.append(sigma_x)
            sigmay.append(sigma_y)
            sigmaxp.append(sigma_xp)
            sigmayp.append(sigma_yp)
            sigmaxxp.append(sigma_x_xp)
            sigmayyp.append(sigma_y_yp)
            meanx.append(x_m)
            meany.append(y_m)
            meanxp.append(xp_m)
            meanyp.append(yp_m)
            Wgt.append(wgt)
            sigm_emx.append(sigma_emitt_x)
            sigm_emy.append(sigma_emitt_y)
            sigm_btx.append(sigma_beta_x)
            sigm_bty.append(sigma_beta_y)
            sigm_alx.append(sigma_alph_x)
            sigm_aly.append(sigma_alph_y)

        #prepare header    
        header = ['Segment','S[m]','Beta_x[m]','Beta_y[m]','Alph_x','Alph_y']
        header.extend(['Disp_x','Disp_xp','Disp_y','Disp_yp'])
        header.extend(['Emitt_x','Emitt_y'])
        header.extend(['Sigma_x[m]','Sigma_y[m]','Sigma_xp[rad]','Sigma_yp[rad]'])
        header.extend(['Sigma_x_xp[m*rad]','Sigma_y_yp[m*rad]'])
        header.extend(['Mean_x[m]','Mean_y[m]','Mean_xp[rad]','Mean_yp[rad]','Wgt'])
        header.extend(['Sigma_emitt_x','Sigma_emitt_y','Sigma_beta_x','Sigma_beta_y'])
        header.extend(['Sigma_alph_x','Sigma_alph_y'])


        #prepare optical function arrays for writing to file

        for i in range(len(S)):
            """
            S[i]=float("{0:.4e}".format(S[i]))
            betx[i]=float("{0:.4e}".format(betx[i]))
            bety[i]=float("{0:.4e}".format(bety[i]))
            alphx[i]=float("{0:.4e}".format(alphx[i]))
            alphy[i]=float("{0:.4e}".format(alphy[i]))
            dispx[i]=float("{0:.4e}".format(dispx[i]))
            dispy[i]=float("{0:.4e}".format(dispy[i]))
            dispxp[i]=float("{0:.4e}".format(dispxp[i]))
            dispyp[i]=float("{0:.4e}".format(dispyp[i]))
            emittx[i]=float("{0:.4e}".format(emittx[i]))
            emitty[i]=float("{0:.4e}".format(emitty[i]))
            sigmax[i]=float("{0:.4e}".format(sigmax[i]))
            sigmay[i]=float("{0:.4e}".format(sigmay[i]))
            sigmaxp[i]=float("{0:.4e}".format(sigmaxp[i]))
            sigmayp[i]=float("{0:.4e}".format(sigmayp[i]))
            sigmaxxp[i]=float("{0:.4e}".format(sigmaxxp[i]))
            sigmayyp[i]=float("{0:.4e}".format(sigmayyp[i]))
            meanx[i]=float("{0:.4e}".format(meanx[i]))
            meany[i]=float("{0:.4e}".format(meany[i]))
            meanxp[i]=float("{0:.4e}".format(meanxp[i]))
            meanyp[i]=float("{0:.4e}".format(meanyp[i]))
            Wgt[i]=float("{0:.4e}".format(Wgt[i]))
            """

        with open(output,'w') as ofile:        
            writer=csv.writer(ofile, delimiter='\t',lineterminator='\n',)
            timestamp = time.strftime("%Y/%m/%d-%H:%M:%S")
            writer.writerow(['# ','Optical functions from PTC output', timestamp])
            writer.writerow(header)
            for i in range(len(S)):
                row = [sampler[i],S[i],betx[i],bety[i],alphx[i],alphy[i],dispx[i],dispy[i]]
                row.extend([dispxp[i],dispyp[i],emittx[i],emitty[i]])
                row.extend([sigmax[i],sigmay[i],sigmaxp[i],sigmayp[i],sigmaxxp[i],sigmayyp[i]])
                row.extend([meanx[i],meany[i],meanxp[i],meanyp[i],Wgt[i]])
                row.extend([sigm_emx[i],sigm_emy[i],sigm_btx[i],sigm_bty[i]])
                row.extend([sigm_alx[i],sigm_aly[i]])
                writer.writerow(row)

            
        
