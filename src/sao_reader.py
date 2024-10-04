"""
SAO-file reading functions

Documentation: https://ulcar.uml.edu/~iag/SAO-4.htm
"""
import numpy as np
from scipy import interpolate
from calendar import timegm
from time import asctime,gmtime,strftime
#, mktime, ctime
from datetime import datetime

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

import os

################################################################################
def extract_critical_frequencies(scaled,index):
	foE = scaled[index][6] if scaled[index][6] < 9999 and not np.isnan(scaled[index][6]) else None
	return foE#, foF2

def Ne_prof(filename,flag=''):
    """
	output:
	- frequencies [MHz]
	- virtual heights [Km]
	- foE/E layer critical frequency [MHz]
	"""
    [tt,h,scaled,Ne,vhot,ftot] = read_sao_data(filename)
    gpath = filename[:-4]
    if not os.path.exists(gpath):
        os.mkdir(gpath)

    path = os.getcwd() + '/' + filename[:-4]

    n = len(tt)
    index = 0
    vhot = np.array(vhot)
    ftot = np.array(ftot)
	
    map_day = {}
    dates = []

    print(filename, n)
    for i in range(n):
        if not(np.isnan(tt[i])):
            if (i == n-1) or (flag != 'last'):
                date = asctime(gmtime(tt[i]))
                dates.append(date)
                # Display foE in the title, and handle the None case
                foE = extract_critical_frequencies(scaled,i)
                
                #-------------------------- Plot stuff -----------------------------#
                # if foE is not None:
                #     plt.title(date + f' | foE = {foE}')
                # else:
                #     plt.title(date + ' | foE = None')

                # ne_i = np.array(Ne[i]) * 1E-4
                # plt.ylabel('Height [km]')
                # plt.xlabel('Freq [$MHz$]')
                # plt.plot(ftot[i], vhot[i], 'ro-', lw=1.8, ms=2.5)
                # plt.grid()
                # npath = path + '/' + str(i) + '.png'
                # plt.savefig(npath)
                # plt.close()
				#-------------------------- ---------- -----------------------------#
                map_day[date] = (ftot[i], vhot[i], foE)

    return map_day, dates
################################################################################
"""

"""

def get_dimensions(lst):
    if isinstance(lst, list):
        return [len(lst)] + get_dimensions(lst[0]) if lst else [0]
    return []

def Ne_map(filename):

	[tt,h,scaled,Ne,vhot,ftot] = read_sao_data(filename)

	tt = np.array(tt)
	h = np.array(h)

	nt = len(tt)
	nh = len(h)

	t0 = 86400*int(tt[0]/86400)
	Ne_aux = np.zeros([nh,nt])
	t_aux = (tt - t0)/3600.

	for i in range(nt):
		Ne_aux[:,i] = np.array(Ne[i])*1E-4

	plt.figure()
	plt.pcolormesh(t_aux,h,Ne_aux)
	plt.title('Electron density - Date: '+strftime("%a %d %b %Y", gmtime(t0)))
#	plt.title('Electron density - Date: '+ctime(t0))
	plt.xlabel('UTC Time [h]')
	plt.ylabel('Height [km]')
	plt.xlim(0,24)
	plt.xticks(np.arange(6+1)*4)
	plt.ylim(100,1000)
	plt.clim(0,120)
	plt.grid()
	cb = plt.colorbar()
	cb.set_label('Ne x $10^4$ [$e/cm^3$]')
	plt.savefig('Ne_map.png')
	plt.close()


################################################################################
"""

"""

def read_sao_data(filename):

	num_hei = 65*3
	spacing = 15//3
	first_height = 100
	height = list(range(first_height,num_hei*spacing+first_height,spacing))

	block = 0
	Ne = []
	vhot = []
	ftot = []
#	num_freq = 0

	time=[]
	scaled=[]

	fid = open(filename)

	while 1:

			[IDFI,GCONST,SYSDES,IPREF,SCALED,IAF,DTT,IOTF2,OTF2,IOTHF2,OTHF2,IOAF2,IODF2,FTOF2,\
			 IOTF1,OTF1,IOTHF1,OTHF1,IOAF1,IODF1,FTOF1,IOTE,OTE,IOTHE,OTHE,IOAE,IODE,FTOE,IXTF2,\
			 XTF2,IXAF2,IXDF2,FTXF2,IXTF1,XTF1,IXAF1,IXDF1,FTXF1,IXTE,XTE,IXAE,IXDE,FTXE,MEDF,\
			 MEDE,MEDES,THF2,THF1,THE,QPCOEF,THVAL,IEDF,IOTSE,OTSE,IOASE,IODSE,FTOSE,IOTAE,OTAE,\
			 IOAAE,IODAE,FTOAE,HTAB,FTAB,NTAB,QL,DL,IEDFTP,IREAD,IERR] = read_sao(fid)

			if IERR==1: break

			block = block + 1
			IPREF=[0]+IPREF
			IDFI=[0]+IDFI

			yyyy = s2f(IPREF[3:7])
			doy = s2f(IPREF[7:10])
			mmm = s2f(IPREF[10:12])
			dom = s2f(IPREF[12:14])
			hh = s2f(IPREF[14:16])
			mm = s2f(IPREF[16:18])
			ss = s2f(IPREF[18:20])

			t=datetime(yyyy,mmm,dom,hh,mm,ss)
			#time.append(mktime(t.timetuple()))
			time.append(timegm(t.timetuple()))

			if IDFI[4]>0:
					SCALED=rep(SCALED,9999.)
					if len(SCALED)<50:
							SCALED=SCALED+[np.nan]*(50-len(SCALED))
					scaled.append(SCALED)
			else:
					scaled.append([np.nan]*50)

			auxheight=[i for i in height]
			if IDFI[51]>0:
					NTAB= rep(NTAB,9999.)
					[HTAB,hid] = uni(HTAB)
					NTAB = [NTAB[i] for i in hid]
#					FTAB= rep(FTAB,9999.)
#					[HTAB,hid] = uni(HTAB)
#					FTAB = [FTAB[i] for i in hid]
					Neb = cubspline(HTAB,NTAB,auxheight)
					if len(Neb)<num_hei:
							Neb=Neb+[np.nan]*(num_hei-len(Neb))
					Ne.append(Neb)
			else:
					Ne.append([np.nan]*num_hei)

			if IDFI[7]>0 or IDFI[12]>0 or IDFI[17]>0:
				VHOT = OTE+OTF1+OTF2
				FTOT = FTOE+FTOF1+FTOF2

				VHOT=rep(VHOT,9999.)

				vhot.append(VHOT)
				ftot.append(FTOT)
			else:
				vhot.append([np.nan])
				ftot.append([np.nan])

			maxv=max([len(i) for i in vhot])
			maxf=max([len(i) for i in ftot])
			vhot=[i+[np.nan]*(maxv-len(i)) for i in vhot]
			ftot=[i+[np.nan]*(maxf-len(i)) for i in ftot]

			IPREF=IPREF[1:]
			IDFI=IDFI[1:]

	fid.close()

	return [time,height,scaled,Ne,vhot,ftot]

################################################################################
"""

"""

def read_sao(IU):

	IERR = 1
	IREAD = 0

	IDFI = digscanf(IU,'%3d',[1,80])

	if IDFI=='flag_EOF': return ['']*69 + [1]

	IDFI=[0]+IDFI

	if IDFI[1]<1 or IDFI[1]>16:
		return
	if (IDFI[80]==0):
			FM1 = '%3d'
			FM2 = '%7f'
			FM3 = '%120c'
			FM4 = '%1c'
			FM5 = '%8f'
			FM6 = '%2d'
			FM7 = '%1d'
			FM8 = '%6f'
			FM9 = '%9f'
			FM10 = '%2d'

	if (IDFI[80]==1):
			FM1 = '%3d'
			FM2 = '%7f'
			FM3 = '%120c'
			FM4 = '%1c'
			FM5 = '%8f'
			FM6 = '%2d'
			FM7 = '%1d'
			FM8 = '%6f'
			FM9 = '%11f'
			FM10 = '%2d'

	if (IDFI[80]>=2):
			FM1 = '%8f'
			FM2 = '%7f'
			FM3 = '%120c'
			FM4 = '%1c'
			FM5 = '%8f'
			FM6 = '%2d'
			FM7 = '%1d'
			FM8 = '%8f'
			FM9 = '%11f'
			FM10 = '%3d'
			FM11 = '%8f'
			FM12 = '%20f'

	if (IDFI[1]>0):
			GCONST = digscanf(IU,FM2,[1,IDFI[1]])
	else:
			GCONST = []

	if (IDFI[2]>0):
			SYSDES =IU.readline()
			SYSDES =SYSDES[0:len(SYSDES)-2]
	else:
			SYSDES = ''

	if (IDFI[2]==2):
			OPMSG = IU.readline()
			OPMSG = OPMSG[0:len(SYSDES)-2]
	else:
			OPMSG = ''

	if (IDFI[3]>0):
			IPREF = digscanf(IU,FM4,[1,IDFI[3]])
	else:
			IPREF = []

	if (IDFI[4]>0):
			SCALED = digscanf(IU,FM5,[1,IDFI[4]])
	else:
			SCALED = []

	if (IDFI[5]>0):
			IAF = digscanf(IU,FM6,[1,IDFI[5]])
	else:
			IAF = []

	if (IDFI[6]>0):
			DTT = digscanf(IU,FM2,[1,IDFI[6]])
	else:
			DTT = []

	if (IDFI[7]>0):
			if (IDFI[80]>=2):
					OTF2 = digscanf(IU,FM1,[1,IDFI[7]])
					auxOTF2 = ['']*len(OTF2)
					for i in range(len(OTF2)): auxOTF2[i] = float(int(OTF2[i]))
					IOTF2 = auxOTF2
			else:
					IOTF2 = digscanf(IU,FM1,[1,IDFI[7]])
					auxIOTF2 = ['']*len(IOTF2)
					for i in range(len(IOTF2)): auxIOTF2[i] = float(IOTF2[i])
					OTF2 = auxIOTF2
	else:
			OTF2 = []
			IOTF2 = []

	if (IDFI[8]>0):
			if (IDFI[80]>=2):
					OTHF2 = digscanf(IU,FM1,[1,IDFI[8]])
					auxOTHF2 = ['']*len(OTHF2)
					for i in range(len(OTHF2)): auxOTHF2[i] = float(int(OTHF2[i]))
					IOTHF2 = auxOTHF2
			else:
					IOTHF2 = digscanf(IU,FM1,[1,IDFI[8]])
					auxIOTHF2 = ['']*len(IOTHF2)
					for i in range(len(IOTHF2)): auxIOTHF2[i] = float(IOTHF2[i])
					OTHF2 = auxIOTHF2

	else:
			OTHF2 = []
			IOTHF2 = []

	if (IDFI[9]>0):
			IOAF2 = digscanf(IU,FM10,[1,IDFI[9]])
	else:
			IOAF2 = []

	if (IDFI[10]>0):
			IODF2 = digscanf(IU,FM7,[1,IDFI[10]])
	else:
			IODF2 = []

	if (IDFI[11]>0):
			FTOF2 = digscanf(IU,FM8,[1,IDFI[11]])
	else:
			FTOF2 = []

	if (IDFI[12]>0):
			if (IDFI[80]>=2):
					OTF1 = digscanf(IU,FM1,[1,IDFI[12]])
					auxOTF1 = ['']*len(OTF1)
					for i in range(len(OTF1)): auxOTF1[i] = float(int(OTF1[i]))
					IOTF1 = auxOTF1
			else:
					IOTF1 = digscanf(IU,FM1,[1,IDFI[12]])
					auxIOTF1 = ['']*len(IOTF1)
					for i in range(len(IOTF1)): auxIOTF1[i] = float(IOTF1[i])
					OTF1 = auxIOTF1
	else:
			OTF1 = []
			IOTF1 = []

	if (IDFI[13]>0):
			if (IDFI[80]>=2):
					OTHF1 = digscanf(IU,FM1,[1,IDFI[13]])
					auxOTHF1 = ['']*len(OTHF1)
					for i in range(len(OTHF1)): auxOTHF1[i] = float(int(OTHF1[i]))
					IOTHF1 = auxOTHF1
			else:
					IOTHF1 = digscanf(IU,FM1,[1,IDFI[13]])
					auxIOTHF1 = ['']*len(IOTHF1)
					for i in range(len(IOTHF1)): auxIOTHF1[i] = float(IOTHF1[i])
					OTHF1 = auxIOTHF1

	else:
			OTHF1 = []
			IOTHF1 = []

	if (IDFI[14]>0):
			IOAF1 = digscanf(IU,FM10,[1,IDFI[14]])
	else:
			IOAF1 = []

	if (IDFI[15]>0):
			IODF1 = digscanf(IU,FM7,[1,IDFI[15]])
	else:
			IODF1 = []

	if (IDFI[16]>0):
			FTOF1 = digscanf(IU,FM8,[1,IDFI[16]])
	else:
			FTOF1 = []

	if (IDFI[17]>0):
			if (IDFI[80]>=2):
					OTE = digscanf(IU,FM1,[1,IDFI[17]])
					auxOTE = ['']*len(OTE)
					for i in range(len(OTE)): auxOTE[i] = float(int(OTE[i]))
					IOTE = auxOTE
			else:
					IOTE = digscanf(IU,FM1,[1,IDFI[17]])
					auxIOTE = ['']*len(IOTE)
					for i in range(len(IOTE)): auxIOTE[i] = float(IOTE[i])
					OTE = auxIOTE
	else:
			OTE = []
			IOTE = []


	if (IDFI[18]>0):
			if (IDFI[80]>=2):
					OTHE = digscanf(IU,FM1,[1,IDFI[18]])
					auxOTHE = ['']*len(OTHE)
					for i in range(len(OTHE)): auxOTHE[i] = float(int(OTHE[i]))
					IOTHE = auxOTHE
			else:
					IOTHE = digscanf(IU,FM1,[1,IDFI[18]])
					auxIOTHE = ['']*len(IOTHE)
					for i in range(len(IOTHE)): auxIOTHE[i] = float(IOTHE[i])
					OTHE = auxIOTHE
	else:
			OTHE = []
			IOTHE = []

	if (IDFI[19]>0):
			IOAE = digscanf(IU,FM10,[1,IDFI[19]])
	else:
			IOAE = []

	if (IDFI[20]>0):
			IODE = digscanf(IU,FM7,[1,IDFI[20]])
	else:
			IODE = []

	if (IDFI[21]>0):
			FTOE = digscanf(IU,FM8,[1,IDFI[21]])
	else:
			FTOE = []

	if (IDFI[22]>0):
			if (IDFI[80]>=2):
					XTF2 = digscanf(IU,FM1,[1,IDFI[22]])
					auxXTF2 = ['']*len(XTF2)
					for i in range(len(XTF2)): auxXTF2[i] = float(int(XTF2[i]))
					IXTF2 = auxXTF2
			else:
					IXTF2 = digscanf(IU,FM1,[1,IDFI[22]])
					auxIXTF2 = ['']*len(IXTF2)
					for i in range(len(IXTF2)): auxIXTF2[i] = float(IXTF2[i])
					XTF2 = auxIXTF2
	else:
			XTF2 = []
			IXTF2 = []

	if (IDFI[23]>0):
			IXAF2 = digscanf(IU,FM10,[1,IDFI[23]])
	else:
			IXAF2 = []

	if (IDFI[24]>0):
			IXDF2 = digscanf(IU,FM7,[1,IDFI[24]])
	else:
			IXDF2 = []

	if (IDFI[25]>0):
			FTXF2 = digscanf(IU,FM8,[1,IDFI[25]])
	else:
			FTXF2 = []

	if (IDFI[26]>0):
			if (IDFI[80]>=2):
					XTF1 = digscanf(IU,FM1,[1,IDFI[26]])
					auxXTF1 = ['']*len(XTF1)
					for i in range(len(XTF1)): auxXTF1[i] = float(int(XTF1[i]))
					IXTF1 = auxXTF1
			else:
					IXTF1 = digscanf(IU,FM1,[1,IDFI[26]])
					auxIXTF1 = ['']*len(IXTF1)
					for i in range(len(IXTF1)): auxIXTF1[i] = float(IXTF1[i])
					XTF1 = auxIXTF1
	else:
			XTF1 = []
			IXTF1 = []

	if (IDFI[27]>0):
			IXAF1 = digscanf(IU,FM10,[1,IDFI[27]])
	else:
			IXAF1 = []

	if (IDFI[28]>0):
			IXDF1 = digscanf(IU,FM7,[1,IDFI[28]])
	else:
			IXDF1 = []

	if (IDFI[29]>0):
			FTXF1 = digscanf(IU,FM8,[1,IDFI[29]])
	else:
			FTXF1 = []

	if (IDFI[30]>0):
			if (IDFI[80]>=2):
					XTE = digscanf(IU,FM1,[1,IDFI[30]])
					auxXTE = ['']*len(XTE)
					for i in range(len(XTE)): auxXTE[i] = float(int(XTE[i]))
					IXTE = auxXTE
			else:
					IXTE = digscanf(IU,FM1,[1,IDFI[30]])
					auxIXTE = ['']*len(IXTE)
					for i in range(len(IXTE)): auxIXTE[i] = float(IXTE[i])
					XTE = auxIXTE
	else:
			XTE = []
			IXTE = []

	if (IDFI[31]>0):
			IXAE = digscanf(IU,FM10,[1,IDFI[31]])
	else:
			IXAE = []

	if (IDFI[32]>0):
			IXDE = digscanf(IU,FM7,[1,IDFI[32]])
	else:
			IXDE = []

	if (IDFI[33]>0):
			FTXE = digscanf(IU,FM8,[1,IDFI[33]])
	else:
			FTXE = []

	if (IDFI[34]>0):
			MEDF = digscanf(IU,FM6,[1,IDFI[34]])
	else:
			MEDF = []

	if (IDFI[35]>0):
			MEDE = digscanf(IU,FM6,[1,IDFI[35]])
	else:
			MEDE = []

	if (IDFI[36]>0):
			MEDES = digscanf(IU,FM6,[1,IDFI[36]])
	else:
			MEDES = []

	if (IDFI[37]>0):
			THF2 = digscanf(IU,FM9,[1,IDFI[37]])
	else:
			THF2 = []

	if (IDFI[38]>0):
			THF1 = digscanf(IU,FM9,[1,IDFI[38]])
	else:
			THF1 = []

	if (IDFI[39]>0):
			THE = digscanf(IU,FM9,[1,IDFI[39]])
	else:
			THE = []

	if (IDFI[40]>0):
			if (IDFI[80]<2):
					THVAL = digscanf(IU,FM9,[1,IDFI[40]])
					QPCOEF = []
			else:
					QPCOEF = digscanf(IU,FM12,[1,IDFI[40]])
					THVAL = []
	else:
			THVAL = []
			QPCOEF = []

	if (IDFI[41]>0):
			IEDF = digscanf(IU,FM7,[1,IDFI[41]])
	else:
			IEDF = []

	if (IDFI[42]>0):
			THVAL = digscanf(IU,FM9,[1,IDFI[42]])
	else:
			THVAL = []

	if (IDFI[43]>0):
			if (IDFI[80]>=2):
					OTSE = digscanf(IU,FM1,[1,IDFI[43]])
					auxOTSE = ['']*len(OTSE)
					for i in range(len(OTSE)): auxOTSE[i] = float(int(OTSE[i]))
					IOTSE = auxOTSE
			else:
					IOTSE = digscanf(IU,FM1,[1,IDFI[43]])
					auxIOTSE = ['']*len(IOTSE)
					for i in range(len(IOTSE)): auxIOTSE[i] = float(IOTSE[i])
					OTSE = auxIOTSE
	else:
			OTSE = []
			IOTSE = []

	if (IDFI[44]>0):
			IOASE = digscanf(IU,FM10,[1,IDFI[44]])
	else:
			IOASE = []

	if (IDFI[45]>0):
			IODSE = digscanf(IU,FM7,[1,IDFI[45]])
	else:
			IODSE = []

	if (IDFI[46]>0):
			FTOSE = digscanf(IU,FM8,[1,IDFI[46]])
	else:
			FTOSE = []

	if (IDFI[47]>0):
			if (IDFI[80]>=2):
					OTAE = digscanf(IU,FM1,[1,IDFI[47]])
					IOTAE=IOTSE
			else:
					IOTAE = digscanf(IU,FM1,[1,IDFI[47]])
					OTAE = OTSE
	else:
			OTAE = []
			IOTAE = []

	if (IDFI[48]>0):
			IOAAE = digscanf(IU,FM10,[1,IDFI[48]])
	else:
			IOAAE = []

	if (IDFI[49]>0):
			IODAE = digscanf(IU,FM7,[1,IDFI[49]])
	else:
			IODAE = []

	if (IDFI[50]>0):
			FTOAE = digscanf(IU,FM8,[1,IDFI[50]])
	else:
			FTOAE = []

	if (IDFI[51]>0):
			HTAB = digscanf(IU,FM8,[1,IDFI[51]])
			FTAB = digscanf(IU,FM8,[1,IDFI[52]])
			NTAB = digscanf(IU,FM11,[1,IDFI[53]])
	else:
			HTAB = []
			FTAB = []
			NTAB = []

	if (IDFI[54]>0):
			QL = digscanf(IU,FM4,[1,IDFI[54]])
	else:
			QL = []

	if (IDFI[55]>0):
			DL = digscanf(IU,FM4,[1,IDFI[55]])
	else:
			DL = []

	if (IDFI[56]>0):
			IEDFTP = digscanf(IU,FM7,[1,IDFI[56]])
	else:
			IEDFTP = []

	IERR = 0
	IREAD = 1

	del(IDFI[0])

	return  [IDFI,GCONST,SYSDES,IPREF,SCALED,IAF,DTT,IOTF2,OTF2,IOTHF2,\
			OTHF2,IOAF2,IODF2,FTOF2,IOTF1,OTF1,IOTHF1,OTHF1,IOAF1,\
			IODF1,FTOF1,IOTE,OTE,IOTHE,OTHE,IOAE,IODE,FTOE,IXTF2,\
			XTF2,IXAF2,IXDF2,FTXF2,IXTF1,XTF1,IXAF1,IXDF1,FTXF1,IXTE,\
			XTE,IXAE,IXDE,FTXE,MEDF,MEDE,MEDES,THF2,THF1,THE,QPCOEF,\
			THVAL,IEDF,IOTSE,OTSE,IOASE,IODSE,FTOSE,IOTAE,OTAE,IOAAE,\
			IODAE,FTOAE,HTAB,FTAB,NTAB,QL,DL,IEDFTP,IREAD,IERR]

################################################################################
"""

"""

import numpy as np

def digscanf(fid, format, size):
    total = np.prod(size)
    datalen = int(np.floor(float(format[1:len(format)-1])))
    linelen = int(datalen * np.floor(120 / datalen))
    numline = int(np.ceil(datalen * total / linelen))
    aux = [''] * linelen * numline

    for i in range(numline):
        tmp = fid.readline()
        if tmp == '':
            break
        tmp = tmp.rstrip('\n')
        lim = min(linelen, len(tmp))
        aux[(linelen * i):(linelen * i + lim)] = tmp[0:lim]

    data = int(datalen)
    A = [[] for _ in range(total)]

    if tmp == '':
        return 'flag_EOF'

    for i in range(0, len(aux), data):
        temp = ''
        if i // data >= total:
            break
        for j in range(data):
            temp += aux[i + j]
        A[i // data] = temp

    ft = format[-1]

    if ft == 'f':
        A = [float(x) for x in A]
    elif ft == 'd':
        A = [int(x) for x in A]

    return A


################################################################################
"""

"""

def cubspline(x,y,xi):

	z=sorted(zip(x,y))
	len_z=len(z)
	x=[z[i][0] for i in range(len_z)]
	y=[z[i][1] for i in range(len_z)]

	x_min=min(x)
	x_max=max(x)

	for i in xi:
		if i<x_min or i>x_max:
			xi[xi.index(i)]=np.nan

	tck=interpolate.splrep(x,y,s=0)
	yi=interpolate.splev(xi,tck,der=0)

	return yi

################################################################################
"""

"""

def uni(x):
	aux=list(set(x))
	x.reverse()
	lenx=len(x)
	len_aux=len(aux)
	ind=[[]]*len_aux

	for i in range(len_aux):
		ind[i]=lenx-x.index(aux[i])-1

	x.reverse()
	return [aux,ind]

################################################################################
"""

"""

def s2f(x):
	tmp=''
	for i in x: tmp=tmp+i
	return int(tmp)

################################################################################
"""

"""

def rep(x,a):
	for i in range(len(x)):
		if x[i]==a: x[i]=np.nan
	return x
