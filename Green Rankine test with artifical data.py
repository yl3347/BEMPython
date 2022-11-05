#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math 
NEXT_NODE = np.array([1,2,3,0]) #Minus one for indexing

'''
Note that Fortran array index different than Python M(1:3) => M[0:3]
'''
# Unneeded part:
# ! Inputs
#     REAL(KIND=PRE), DIMENSION(3),    INTENT(IN) :: M
#     REAL(KIND=PRE), DIMENSION(4, 3), INTENT(IN) :: Face_nodes
#     REAL(KIND=PRE), DIMENSION(3),    INTENT(IN) :: Face_center, Face_normal
#     REAL(KIND=PRE),                  INTENT(IN) :: Face_area, Face_radius
# 
# ! Outputs
#     REAL(KIND=PRE),               INTENT(OUT) :: S0
#     REAL(KIND=PRE), DIMENSION(3), INTENT(OUT) :: VS0

def COMPUTE_INTEGRAL_OF_RANKINE_SOURCE(M,Face_nodes, Face_center, Face_normal, Face_area, Face_radius,S0, VS0):
    L = 0                                               # INTEGER                         :: L
    RO, GZ, DK, GY = float(0),float(0),float(0),float(0)# REAL(KIND=PRE)                  :: RO, GZ, DK, GY
    RR = np.zeros(4)                                    # REAL(KIND=PRE), DIMENSION(4)    :: RR
    DRX = np.zeros((3,4))                               # REAL(KIND=PRE), DIMENSION(3, 4) :: DRX
    ANT = DNT = ANL = DNL = ALDEN = AT = float(0)        # REAL(KIND=PRE)                  :: ANT, DNT, ANL, DNL, ALDEN, AT
    PJ = GYX = ANTX = ANLX = DNTX = np.zeros(3)         # REAL(KIND=PRE), DIMENSION(3)    :: PJ, GYX, ANTX, ANLX, DNTX

    # Distance from center of mass of the face to M.
    RO = np.linalg.norm(M[0:3] - Face_center[0:3])      # RO = NORM2(M(1:3) - Face_center(1:3)) 

    if RO > 7 * Face_radius:                            # IF (RO > 7*Face_radius) THEN
        # Asymptotic value if face far away from M
        S0 = Face_area/RO                               # S0       = Face_area/RO
        VS0[0:3] = (Face_center[0:3] - M) * S0 / RO ** 2# VS0(1:3) = (Face_center(1:3) - M)*S0/RO**2
        
    else:
        # Called Z in [Del]
        GZ = np.dot((M[0:3] - Face_center[0:3]), Face_normal[0:3])  # GZ = DOT_PRODUCT(M(1:3) - Face_center(1:3), Face_normal(1:3))

        for L in range(4):                                     # DO CONCURRENT (L = 1:4)
            # Distance from vertices of Face to M.
            RR[L] = np.linalg.norm(M[0:3] - Face_nodes[L,0:3])     # RR(L) = NORM2(M(1:3) - Face_nodes(L, 1:3))
            # Normed vector from vertices of Face to M.
            DRX[:,L] = (M[0:3] - Face_nodes[L,0:3])/RR[L]          # DRX(:, L) = (M(1:3) - Face_nodes(L, 1:3))/RR(L)
                                                                    # END DO

        S0 = VS0[:] = 0

        for L in range(4):
            # Distance between two consecutive points, called d_k in [Del]
            DK = np.linalg.norm(Face_nodes[NEXT_NODE[L],:] - Face_nodes[L,:])       # DK = NORM2(Face_nodes(NEXT_NODE(L), :) - Face_nodes(L, :))

            if DK >= float(1e-3) * Face_radius :        # IF (DK >= REAL(1e-3, PRE)*Face_radius) THEN\
                # Normed vector from one corner to the next
                PJ[:] = (Face_nodes[NEXT_NODE[L],:] - Face_nodes[L,:]) / DK          # PJ(:) = (Face_nodes(NEXT_NODE(L), :) - Face_nodes(L, :))/DK

                # The following GYX(1:3) are called (a,b,c) in [Del]
                GYX[0] = Face_normal[1] * PJ[2] - Face_normal[2] * PJ[1]            # GYX(1) = Face_normal(2)*PJ(3) - Face_normal(3)*PJ(2)
                GYX[1] = Face_normal[2] * PJ[0] - Face_normal[0] * PJ[2]            # GYX(2) = Face_normal(3)*PJ(1) - Face_normal(1)*PJ(3)
                GYX[2] = Face_normal[0] * PJ[1] - Face_normal[1] * PJ[0]            # GYX(3) = Face_normal(1)*PJ(2) - Face_normal(2)*PJ(1)
                # Called Y_k in  [Del]
                GY = np.dot( M - Face_nodes[L,:] , GYX)                             # GY = DOT_PRODUCT(M - Face_nodes(L, :), GYX)  

                # Called N^t_k in [Del]
                ANT = 2 * GY * DK                                                   # ANT = 2*GY*DK 
                # Called D^t_k in [Del]
                DNT = (RR[NEXT_NODE[L]]+RR[L])**2 - DK * DK + 2 * np.abs(GZ) * (RR[NEXT_NODE[L]] + RR[L])   # DNT = (RR(NEXT_NODE(L))+RR(L))**2 - DK*DK + 2*ABS(GZ)*(RR(NEXT_NODE(L))+RR(L)) 
                # Called N^l_k in [Del]
                ANL = RR[NEXT_NODE[L]] + RR[L] + DK                                 # ANL = RR(NEXT_NODE(L)) + RR(L) + DK
                # Called D^l_k in [Del]
                DNL = RR[NEXT_NODE[L]] + RR[L] - DK                                 # DNL = RR(NEXT_NODE(L)) + RR(L) - DK 
                # Called D^l_k in [Del]
                #using abs value for testing purpose
                ALDEN = math.log(np.abs(ANL/DNL))                                            # ALDEN = LOG(ANL/DNL)

                if np.abs(GZ) >= float(1e-4) * Face_radius :                        # IF (ABS(GZ) >= REAL(1e-4, PRE)*Face_radius) THEN
                    AT = np.arctan(ANT/DNT)                                         # AT = ATAN(ANT/DNT)
                else:                                                               # ELSE
                    AT = 0.0                                                        # AT = 0.
                                                                                    # ENDIF
                # Called N^l_k_{x,y,z} in [Del]
                ANLX[:] = DRX[:,NEXT_NODE[L]] + DRX[:,L]                            # ANLX(:) = DRX(:, NEXT_NODE(L)) + DRX(:, L)
                
                # Called N^t_k_{x,y,z} in [Del]
                ANTX[:] = 2 * DK * GYX[:]                                           # ANTX(:) = 2*DK*GYX(:)
                # Called D^t_k_{x,y,z} in [Del]
                #What is ONE? Replacing with 1 for testing now
                DNTX[:] = 2 *(RR[NEXT_NODE[L]] + RR[L]) + np.abs(GZ) * ANLX[:] + 2 * np.copysign(1,GZ) * (RR[NEXT_NODE[L]] + RR[L]) * Face_normal[:] #DNTX(:) = 2*(RR(NEXT_NODE(L)) + RR(L) + ABS(GZ))*ANLX(:) + 2*SIGN(ONE, GZ)*(RR(NEXT_NODE(L)) + RR(L))*Face_normal(:)

                if np.abs(GY) < 1e-5 :
                    # Edge case where the singularity is on the boundary of the face (GY = 0, ALDEN = infty).
                    # This case seems to only occur when computating the free surface elevation,
                    # so no fix has been implemented for VS0, which is not needed then.

                    S0 = S0 - 2 * AT * np.abs(GZ)                                   # S0 = S0 - 2*AT*ABS(GZ)

                else:
                    # general case

                    S0 = S0 + GY * ALDEN - 2 * AT * np.abs(GZ)                      # S0 = S0 + GY*ALDEN - 2*AT*ABS(GZ)
                                                                                    # END IF
                temp=VS0
                VS0[:] = temp[:] + ALDEN * GYX[:] - 2 * np.copysign(1,GZ) * AT * Face_normal[:] + GY * (DNL -  ANL) / (ANL*DNL)*ANLX[:] - 2 * np.abs(GZ) * (ANTX[:] * DNT - DNTX[:] * ANT) / (ANT**2 + DNT**2)
                # VS0(:) = VS0(:) + ALDEN*GYX(:) - 2*SIGN(ONE, GZ)*AT*Face_normal(:) + GY*(DNL-ANL)/(ANL*DNL)*ANLX(:) - 2*ABS(GZ)*(ANTX(:)*DNT - DNTX(:)*ANT)/(ANT*ANT+DNT*DNT)
            # END IF
        # END DO
    # END IF


# In[3]:


import numpy as np
import math

#PURE SUBROUTINE COMPUTE_ASYMPTOTIC_RANKINE_SOURCE(M,Face_center, Face_area,S0, VS0)
def COMPUTE_ASYMPTOTIC_RANKINE_SOURCE(M,Face_center, Face_area,S0, VS0):
#     Unneeded part from Fortran:
#     ! Inputs
#     REAL(KIND=PRE), DIMENSION(3), INTENT(IN) :: M
#     REAL(KIND=PRE), DIMENSION(3), INTENT(IN) :: Face_center
#     REAL(KIND=PRE),               INTENT(IN) :: Face_area

#     ! Outputs
#     REAL(KIND=PRE),               INTENT(OUT) :: S0
#     REAL(KIND=PRE), DIMENSION(3), INTENT(OUT) :: VS0

#     Local Variable
    RO = 0 # REAL(KIND=PRE) :: RO
    
#     Distance from center of mass of the face to M.
    RO = np.linalg.norm(M[0:3] - Face_center[0:3]) # RO = NORM2(M(1:3) - Face_center(1:3))
    
    if RO > 1e-7 :# IF (RO > REAL(1e-7, KIND=PRE)) THEN
#         Asymptotic value if face far away from M
        S0= Face_area/RO #S0       = Face_area/RO
        VS0[0:3] = (Face_center[0:3]-M)*S0/RO**2 # VS0(1:3) = (Face_center(1:3) - M)*S0/RO**2
    else:
#         Singularity...
        S0 = 0 # S0 = ZERO
        VS0[0:3] = 0 # VS0(1:3) = ZERO
    
        
    


# In[15]:


get_ipython().run_cell_magic('time', '', 'M=np.array([3,5,7])\nFN=np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])\nFC=np.array([2,4,6])\nFNorm=np.array([4,8,10])\nFA= 5\nFR=2.5\nS0=1\nVS0=np.array([2,5,9])\nCOMPUTE_ASYMPTOTIC_RANKINE_SOURCE(M,FC,FA,S0,VS0)\nprint(VS0)\nCOMPUTE_INTEGRAL_OF_RANKINE_SOURCE(M,FN,FC,FNorm,FA,FR,S0,VS0)\nprint(VS0)')


# In[ ]:




