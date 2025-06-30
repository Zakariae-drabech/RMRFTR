import numpy as np
import math

from scipy import stats

import scipy.linalg as la
import scipy.signal as ss 
from time import time

class RMRFTR:
    """
    Recursive Markov Random Fields Trend Restore (RMRFTR)     
    
    Parameters:
    -----------
    mu : float
        Smoothing parameter
    h0 : float  
        threshold  for detection of an abrupt jump
    MAX_nbr_CPs : int, optional
        Maximum number of change points to detect
    discontinue : bool, default=True
        Whether to allow discontinuous jumps at change points
    min_segment : int, default=3
        Minimum segment length between change points     
    """

    def __init__(self, mu, h0, MAX_nbr_CPs=None, discontinue=True,min_segment=3):
        # Input validation
        if mu <= 0:
            raise ValueError("mu must be positive")
        if h0 <= 0:
            raise ValueError("h0 must be positive")
        if min_segment < 3:
            raise ValueError("min_segment must be at least 3")
            
        self.mu =mu
        self.mu4 = mu**4
        self.h0 = h0
        self.MAX_nbr_CPs = MAX_nbr_CPs if MAX_nbr_CPs is not None else 100
        self.discontinue = discontinue
        self.min_segment=min_segment

        self.beta = mu * h0**2 / (2*np.sqrt(2)) 

        self.ChangePoints=[]
        self.y_in=None
        self.x_mrf=None

    def _A_mat_lap(self,n,l=None):
        if l is None:
            l=np.zeros(n)
        mu4=self.mu4
        D=np.zeros(n) 
        D[2:-2]+=1+(6-l[3:-1]-l[1:-3]-4*l[2:-2])*mu4
        D[0]=1+(1-l[0])*mu4
        D[1]=1+(5-4*l[1]-l[2])*mu4
        D[-1]=1+(1-l[-1])*mu4
        D[-2]=1+(5-4*l[-1]-l[-2])*mu4
        U1=np.zeros(n-1)
        U1[1:-1]=-(4-2*l[2:-1]-2*l[1:-2])*mu4
        U1[0]=-(2-2*l[1])*mu4
        U1[-1]=-(2-2*l[-1])*mu4
        U2=np.zeros(n-2)+(1-l[1:-1])*mu4
        A=np.vstack((np.append([0,0],U2),np.append(0,U1),D,np.append(U1,0),np.append(U2,[0,0])))
        #A=np.diag(U2,-2)+np.diag(U2,2)+np.diag(U1,-1)+np.diag(U1,1)+np.diag(D,0)
        return A
    def _computeEnergy_cont(self,y):
        n=len(y)
        mu4=self.mu4
        beta=self.beta/2
        l=np.zeros(n)
        E=np.full(n, np.inf)
        r=self.min_segment
        for i in range(r,n-r):
            l_rup=l.copy()
            l_rup[i]=1
            A = self._A_mat_lap(n,l_rup)
            x=la.solve_banded ((2,2),A,y)
            E[i]=np.linalg.norm(x-y)**2+mu4*np.sum(((x[2:]-2*x[1:-1]+x[0:-2])**2)*(1-l_rup[1:-1]))+beta
        return E
    def _computeEnergy_dis(self,y):
        """ Brute force Solver  (Algorithm 1)"""
        n=len(y)
        mu4=self.mu4
        beta=self.beta
        l=np.zeros(n)
        E=np.full(n, np.inf)
        r=self.min_segment
        for i in range(r,n-r):
            l_rup=l.copy()
            l_rup[i-1]=1
            l_rup[i]=1
            Al = self._A_mat_lap(i)
            Ar = self._A_mat_lap(n-i)
            xl=la.solve_banded((2,2),Al,y[:i])
            xr=la.solve_banded((2,2),Ar,y[i:])
            x=np.append(xl,xr)

            E[i]=np.linalg.norm(x-y)**2+mu4*np.sum(((x[2:]-2*x[1:-1]+x[0:-2])**2)*(1-l_rup[1:-1]))+beta
        return E

    def RecursiveMRFTR(self,y,a=0,b=None):
        """  
        Algoritm 2
        
        Parameters:
        -----------
        y : array-like
            Input signal segment
        a : int  
            Starting index offset
        b : int, optional
            Ending index offset
        """
        n=len(y)
        mu4=self.mu4
        A  = self._A_mat_lap(n)
        x=la.solve_banded((2,2),A,y)
        E_prec = np.linalg.norm(x-y)**2+mu4*np.sum(((x[2:]-2*x[1:-1]+x[0:-2])**2))
        if self.discontinue:
            E=self._computeEnergy_dis(y)
        else:
            E=self._computeEnergy_cont(y)
        ik=np.argmin(E)
        Ek=np.min(E)
        min_segment=self.min_segment+2 # make sure that the segment size is at least 3+2 to have a CP
        if Ek<E_prec and len(self.ChangePoints)<self.MAX_nbr_CPs:
            cp=ik+a
            self.ChangePoints.append(cp)
            yl=y[:ik]
            yr=y[ik:]
            if len(yl)>=min_segment:
                self.RecursiveMRFTR(yl,a=a,b=cp)
            if len(yr)>=min_segment:
                    self.RecursiveMRFTR(yr,a=cp,b=b)
 
    def restore(self,y_in):
        """
        Restore signal using detected change points
        
        Parameters:
        -----------
        y_in : array-like
            Input noisy signal
            
        Returns:
        --------
        x_mrf : array 
            Restored signal
        """
        n=len(y_in)
        l=np.zeros(n)
        _ChangePoints_=self.ChangePoints.copy()
        if self.discontinue:
            _ChangePoints_+=[ik-1 for ik in _ChangePoints_]
        l[_ChangePoints_]=1
        A = self._A_mat_lap(n,l=l)
        self.x_mrf=la.solve_banded((2,2),A,y_in)
        self.y_in=y_in
    def restore_withKnownCPs(self, y_in, ChangePoints):
        n=len(y_in)
        l=np.zeros(n)
        if self.discontinue:
            ChangePoints+=[ik-1 for ik in ChangePoints]
        l[ChangePoints]=1
        A = self._A_mat_lap(n,l=l)
        self.x_mrf=la.solve_banded((2,2), A, y_in)
        self.y_in=y_in
    def error(self,x_exact):
        """
        Compute mean squared error against ground truth
        
        Parameters:
        -----------
        x_exact : array-like
            Ground truth signal
            
        Returns:
        --------
        mse : float
            Mean squared error
        """
        return np.mean((x_exact-self.x_mrf)**2)
    def run(self, y):
        """
        detect change points and restore signal
        
        Parameters:
        -----------
        y : array-like
            Input noisy signal
        """
        # Reset state
        self.ChangePoints = []
        
        # Detect change points
        self.RecursiveMRFTR(y)
        
        # Restore signal
        self.restore(y)