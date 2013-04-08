import math

def sig_fig(n,e):
    dec_place = int(math.floor(math.log10(e)))
    if dec_place >= 0:
        val = str(int(round(n/(10**dec_place))*10**dec_place))
        err = str(int(round(e/(10**dec_place),1)*10**dec_place))
        return val+r'\pm'+err
    else:
        val = round(n,-dec_place)
        if err*dec_place < 3:
            err = round(e,-dec_place+1)
        else:
            err = round(e,-dec_place)
        return val + r'\pm' + err
 
