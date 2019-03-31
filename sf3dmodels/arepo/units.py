#********************
#AREPO units (in cgs)
#********************
ulength = 3.0856e20 # 100 pc in cm
umass = 1.991e33 # 1 solar mass in grams
uvel = 1.0e5 # 1e3 m/s in cm/s
utime = ulength/uvel
udensity = umass/ulength/ulength/ulength
uenergy= umass*uvel*uvel
ucolumn = umass/ulength/ulength
umag = umass**0.5 / ulength**1.5 * uvel
uMyr=utime/(3600.*24.*365.25*1.e6)
uparsec=ulength/3.0856e18
kb = 1.3806485e-16

