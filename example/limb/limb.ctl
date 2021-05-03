# ======================================================================
# Forward model...
# ======================================================================

# Table directory...
TBLBASE = ./boxcar

# Emitters...
NG = 5
EMITTER[0] = CO2
EMITTER[1] = H2O
EMITTER[2] = O3
EMITTER[3] = F11
EMITTER[4] = CCl4

# Channels...
ND = 2
NU[0] = 792.0000
NU[1] = 832.0000

# use the GPU: 0:never, 1:always, -1:if possible
USEGPU = 0
