"""Tests for the functions in the model module.
"""
import numpy as np
from model import model as m
from model import ahdata as dC


def test_acm():
    """Test for acm fn in model module.
    """
    d=dC.dalecData(201)
    assert m.acm(0.01, d.p17, d.p11, d, 200) < 1e-3
    

def test_fitpoly():
    """Test for fitpolynomial function in model module.
    """
    assert m.fitpolynomial(0., 0.) == 0.
    assert m.fitpolynomial(0., 1.) == -0.188459767342504


def test_tempterm():
    """Test for temp_term fn in the model module.
    """
    d=dC.dalecData(10)
    assert m.temp_term(0., d, 0) == 1.    
    
    
def test_phionset():
    """Test for phi_onset fn in model module.
    """
    d=dC.dalecData(10)
    assert m.phi_onset(0., 0.001, d, 0) == 0.
    
    
def test_phifall():
    """Test for phi_fall fn in model module.
    """
    d=dC.dalecData(10)
    assert m.phi_fall(0., 0.01, 1.1, d, 0) == 0.
    
    
def test_dalecv2():
    """Test for dalecv2 and dalecv2_inuput fn.
    """
    d=dC.dalecData(10)
    p = d.pvals
    p[0:6] = 0.
    p[1] = 1e-12
    dalecout = m.dalecv2_input(p, d, 0)
    assert np.allclose(p, dalecout, 1e-12, 1e-12)


def test_lindalecv2(lam=1e-4):
    """Test for lin_dalecv2 fn.
    """
    datx = dC.dalecData(10)
    datdx = dC.dalecData(10)
    pvalx = datx.pvals
    pvaldx = datdx.pvals*0.1
    mxdx = m.dalecv2_input(pvalx+lam*pvaldx, datx, 0)
    mx = m.dalecv2_input(pvalx, datx, 0)
    mat0 = m.lin_dalecv2(pvalx, datx, 0)[1]
    print abs(np.linalg.norm(mxdx-mx) / \
             np.linalg.norm(np.dot(mat0, lam*pvaldx)) - 1)  
    assert abs(np.linalg.norm(mxdx-mx) / \
             np.linalg.norm(np.dot(mat0, lam*pvaldx)) - 1) < 1e-8
             
             
def test_linmodev(lam=1e-4):
    """Test for linmod_list and linmod_evolve fns.
    """
    datx = dC.dalecData(10)
    datdx = dC.dalecData(10)
    pvalx = datx.pvals
    pvaldx = datdx.pvals*0.1   
    mxdx = m.mod_list(pvalx+lam*pvaldx, datx, 0, 10)
    mx = m.mod_list(pvalx, datx, 0, 10)
    matlist = m.linmod_list(pvalx, datx, 0, 10)[1]
    linmodev = m.linmod_evolve(lam*pvaldx, matlist, datx, 0, 10)
    print abs(np.linalg.norm(mxdx[10]-mx[10]) / \
             np.linalg.norm(linmodev[10]) - 1)  
    assert abs(np.linalg.norm(mxdx[10]-mx[10]) / \
             np.linalg.norm(linmodev[10]) - 1) < 1e-8
             
             
def test_linmodevfac(lam=1e-4):
    """Test for linmod_list and linmod_evolvefac fns.
    """
    datx = dC.dalecData(10)
    datdx = dC.dalecData(10)
    pvalx = datx.pvals
    pvaldx = datdx.pvals*0.1   
    mxdx = m.mod_list(pvalx+lam*pvaldx, datx, 0, 10)
    mx = m.mod_list(pvalx, datx, 0, 10)
    matlist = m.linmod_list(pvalx, datx, 0, 10)[1]
    linmodev = m.linmod_evolvefac(lam*pvaldx, matlist, datx, 0, 10)
    print abs(np.linalg.norm(mxdx[10]-mx[10]) / \
             np.linalg.norm(linmodev[10]) - 1)  
    assert abs(np.linalg.norm(mxdx[10]-mx[10]) / \
             np.linalg.norm(linmodev[10]) - 1) < 1e-8

