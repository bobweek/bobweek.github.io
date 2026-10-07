"""coev_likelihood.py -- the maximum-likelihood machinery of Week & Nuismer (2019),
transliterated line for line from functions.R in the measuring.coevolution
repository (R ^ -> Python **, sqrt -> np.sqrt; nothing else changed).

ml_sol   analytic ML estimates of A1, A2, B1, B2 and the offset k from the
         five moments (mu1, mu2, V1, V2, C) and the background parameters
mu       equilibrium means given the selection strengths
SIGMA    equilibrium covariance matrix given the selection strengths
"""
import numpy as np

def ml_sol(mu1,mu2,V1,V2,C,n1,n2,th1,th2,G1,G2):
    A1 = (
        (-(G2*(n1)**2*(C - V1)*V1*((C)**2*(th1 - mu1) + (V1)**2*(-th2 + mu2) + C*V1*(-th1 + th2 - mu1 + mu2))) + 
           G1*n1*n2*V1*(V1*(V2)**2*(th1 - 2*th2 + 3*mu1 - 2*mu2) + (C)**3*(-th2 + mu2) + (C)**2*(V1*(-th2 + mu2) + V2*(th1 + th2 - 3*mu1 + mu2)) - C*V2*(V2*(th1 - mu1) + V1*(th1 - 3*th2 + mu1 + mu2))) + 
           C*np.sqrt((n1)**2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)**2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                            (G1*n2*V1*(V2)**2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)**3*(-th2 + mu2) + 
                               (C)**2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))**2)) - 
           V1*np.sqrt((n1)**2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)**2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                             (G1*n2*V1*(V2)**2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)**3*(-th2 + mu2) + 
                                (C)**2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))**2)))/
    (4.*G1*(n1)**2*n2*V1*(-(C)**2 + V1*V2)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)))
    )
    A2 = (
        (G1*n1*n2*V1*(C - V2)*((V2)**2*(-th1 + mu1) + (C)**2*(th2 - mu2) + C*V2*(th1 - th2 + mu1 - mu2)) + 
           G2*(n1)**2*V1*((C)**3*(th1 - mu1) - (C)**2*(V2*(-th1 + mu1) + V1*(th1 + th2 + mu1 - 3*mu2)) + (V1)**2*V2*(2*th1 - th2 + 2*mu1 - 3*mu2) + C*V1*(V1*(th2 - mu2) + V2*(-3*th1 + th2 + mu1 + mu2))) - 
           C*np.sqrt((n1)**2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)**2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                            (G1*n2*V1*(V2)**2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)**3*(-th2 + mu2) + 
                               (C)**2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))**2)) + 
           V2*np.sqrt((n1)**2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)**2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                             (G1*n2*V1*(V2)**2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)**3*(-th2 + mu2) + 
                                (C)**2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))**2)))/
    (4.*G2*(n1)**2*n2*V1*(-(C)**2 + V1*V2)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)))
    )
    B1 = (
        (G2*(n1)**2*V1*((C)**2*(th1 - mu1) + (V1)**2*(-th2 + mu2) + C*V1*(-th1 + th2 - mu1 + mu2)) + 
           G1*n1*n2*(V1*(V2)**2*(-th1 + mu1) + C*V1*V2*(-th1 + th2 - mu1 + mu2) + (C)**2*(2*V2*(th1 - mu1) + V1*(-th2 + mu2))) - 
           np.sqrt((n1)**2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)**2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                          (G1*n2*V1*(V2)**2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)**3*(-th2 + mu2) + 
                             (C)**2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))**2)))/
    (4.*G1*(n1)**2*n2*((C)**2 - V1*V2)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)))
    )
    B2 = (
        (G1*n1*n2*V1*V2*((V2)**2*(-th1 + mu1) + (C)**2*(th2 - mu2) + C*V2*(th1 - th2 + mu1 - mu2)) + 
           G2*(n1)**2*V1*((C)**2*(V2*(-th1 + mu1) + 2*V1*(th2 - mu2)) + C*V1*V2*(th1 - th2 + mu1 - mu2) + (V1)**2*V2*(-th2 + mu2)) - 
           V2*np.sqrt((n1)**2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)**2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                             (G1*n2*V1*(V2)**2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)**3*(-th2 + mu2) + 
                                (C)**2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))**2)))/
    (4.*G2*(n1)**2*n2*V1*(-(C)**2 + V1*V2)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)))
    )
    k = (
        (-(G1*n1*n2*V1*(C - V2)*(V2*(-th1 + mu1) + C*(th2 - mu2))) - G2*(n1)**2*(C - V1)*V1*(C*(th1 - mu1) + V1*(-th2 + mu2)) + 
          np.sqrt((n1)**2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)**2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                         (G1*n2*V1*(V2)**2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)**3*(-th2 + mu2) + 
                            (C)**2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))**2)))/(2.*C*n1*V1*(G2*n1*V1 + G1*n2*V2))
    )
    return dict(A1=float(A1), A2=float(A2), B1=float(B1), B2=float(B2), k=float(k))


def mu(A1,A2,B1,B2,k,n1,n2,th1,th2):
    mu1 = (A1*(A2 + B2)*th1 + B1*(2*B2*k + A2*(th2 + k)))/(A2*B1 + A1*(A2 + B2))
    mu2 = (A1*B2*th1 + A1*A2*th2 + A2*B1*th2 + (A1 + 2*B1)*B2*k)/(A2*B1 + A1*(A2 + B2))
    return np.array([mu1, mu2])


def SIGMA(A1,A2,B1,B2,k,n1,n2,th1,th2,G1,G2):
    COV = (
        (B1*(A1 + B1)*G1*n1 + B2*(A2 + B2)*G2*n2)/
    (2.*(A2*B1 + A1*(A2 + B2))*((A1 + B1)*G1 + (A2 + B2)*G2)*n1*n2)
    )
    VAR1 = (
        ((B1)**2*G1*n1 + A2*B1*G1*n2 + (A2 + B2)*(A1*G1 + (A2 + B2)*G2)*n2)/
    (2.*(A2*B1 + A1*(A2 + B2))*((A1 + B1)*G1 + (A2 + B2)*G2)*n1*n2)
    )
    VAR2 = (
        (((A1 + B1)**2*G1 + (A2*B1 + A1*(A2 + B2))*G2)*n1 + (B2)**2*G2*n2)/
    (2.*(A2*B1 + A1*(A2 + B2))*((A1 + B1)*G1 + (A2 + B2)*G2)*n1*n2)
    )
    return np.array([[VAR1, COV], [COV, VAR2]])
