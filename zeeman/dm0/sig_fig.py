import math

def sig_fig(n,e):
   dec_place = int(math.floor(math.log10(e)))
   if dec_place >= 0:
       return str(int(round(n/(10**dec_place))*10**dec_place))+"\pm"+str(int(round(e/(10**dec_place),1)*10**dec_place))
   else:
#       return str(n)+' +/- '+str(e)
        return str(round(n,-dec_place))+"\pm"+str(round(e,-dec_place+1))

def sig_fig1(n,e):
   dec_place = int(math.floor(math.log10(e)))
   if dec_place >= 0:
       return str(int(round(n/(10**dec_place))*10**dec_place))+"\pm"+str(int(round(e/(10**dec_place),1)*10**dec_place))
   else:
#       return str(n)+' +/- '+str(e)
        return str(round(n,-dec_place))+"\pm"+str(round(e,-dec_place+1))
