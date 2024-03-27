import numpy as np
from scipy.optimize import curve_fit,minimize
import matplotlib.pyplot as plt


def process_data(cfg, data):
    return []

def ppm_to_mol(c_ppm, airpres_hpa, totalvol_ml, airvol_ml, eqtemp_c, R):
    """Convert ppm to mol assuming ideal gas.

    Parameters
    ----------
    c_ppm : float
        Concentration of ppm
    airpres_hpa : float
        Air pressure in hPa
    totalvol_ml : float
        Total volume of the vial in mL
    airvol_ml : float
        Volume of the headspace (gas) in mL
    eqtemp_c : float
        Equilibrium teperature in degC
    R : float
        Ideal gas constant in TODO

    Returns
    -------
    c_mol, c_pa
    c_mol: float
        Gas concentration in mols
    c_pa: float
        Gas partial pressure in Pascal
    """
    c_pa = c_ppm/10000*airpres_hpa
    c_mol = c_pa*(totalvol_ml-airvol_ml)/1000/1000/(R*(eqtemp_c+273.15))
    return c_mol, c_pa

def equibrium(ce_mol, ce_pa, ca_mol, eqtemp_c, hcp25, dlnHcpd1_T, airvol_ml):
    """Calculates concentration in the vial assuming gas equilimbrium between the gas and liquid phase.

    Parameters
    ----------
    ce_mol: float
        
    ce_pa: float

    ca_mol: float
    eqtemp_c : float
        Equilibrium teperature in degC
    hcp25: float
        Henry's coefficient at 25 degC and 1 atm
    dlnHcpd1_T: float
    airvol_ml : float
        Volume of the headspace (gas) in mL

    Returns
    -------
    TODO
    """
    hcp_t = hcp(eqtemp_c, hcp25, dlnHcpd1_T)
    caq_mol = hcp_t*ce_pa*airvol_ml/1000
    ntotal = ce_mol + caq_mol
    c_molm3 = (ntotal - ca_mol)/(airvol_ml/1000/1000)
    c_umolL = c_molm3*1000
    return c_umolL

def hcp(airtemp_c, hcp25, dlnHcpd1_T):
    hcp_t = hcp25*np.exp(dlnHcpd1_T*(1/(airtemp_c+273.15)-1/298.15))/1000
    return hcp_t

def reaction_constant(eqtemp_c):
    eqtemp_k = eqtemp_c + 273.15
    pkw = -(6.0875-4470.99/eqtemp_k-0.01706*eqtemp_k)
    pk1 = 0.000011*eqtemp_c**2-0.012*eqtemp_c+6.58
    pk2 = 0.00009*eqtemp_c**2-0.0137*eqtemp_c+10.62
    return pkw, pk1, pk2

def bottle_equilibruim_dic(pHeq, eqtemp_c):
    pkw, pk1, pk2 = reaction_constant(eqtemp_c)
    a0 = 1/(1+10**(-pk1+pHeq)+10**(-(pk1+pk2)+2*pHeq))
    a1 = 1/(10**(-pHeq+pk1)+1+10**(-pk2+pHeq))
    a2 = 1/(10**(-2*pHeq+pk1+pk2)+10**(-pHeq+pk2)+1)
    return pkw, a0, a1, a2, a0+a1+a2

def carbon_balance(pHeq, ce_mol, eqtemp_c, alk, hcp25, dlnHcpd1_T, airvol_ml):
    pkw, a0, a1, a2, _ = bottle_equilibruim_dic(pHeq, eqtemp_c)
    dic_molm3 = 1000*(alk*0.001 - 10**(-pkw+pHeq)+10**(-pHeq))/(a1+2*a2)
    h2co3_molm3 = a0*dic_molm3
    hco3_molm3 = a1*dic_molm3
    co3_molm3 = a2*dic_molm3
    pco2_pa = h2co3_molm3/(hcp(eqtemp_c, hcp25, dlnHcpd1_T)*1000)
    dic_mol = dic_molm3 * airvol_ml*1E-6+ce_mol
    return pco2_pa, dic_mol, dic_molm3, h2co3_molm3

def opt_fuc_bottle(pHeq, ce_pa, ce_mol, eqtemp_c, alk, hcp25, dlnHcpd1_T, airvol_ml):
    pco2_pa, _, _, _= carbon_balance(pHeq, ce_mol, eqtemp_c, alk, hcp25, dlnHcpd1_T, airvol_ml)
    return (pco2_pa-ce_pa)**2*1E6

def opt_fuc_insitu(pHeq, wsdic_mol, ca_mol, ce_mol, eqtemp_c, alk, hcp25, dlnHcpd1_T, airvol_ml):
    _, _, dic_molm3, _ = carbon_balance(pHeq, ce_mol, eqtemp_c, alk, hcp25, dlnHcpd1_T, airvol_ml)
    wsdic_molm3 = (wsdic_mol-ca_mol)/(airvol_ml*1E-6)
    return (dic_molm3-wsdic_molm3)**2*1E9

def alkalinity_correction(c_ppm, alk, pH, airpres_hpa, volumen_ml, temp_c, constant):
    ce_ppm, ca_ppm = c_ppm
    hcp25, dlnHcpd1_T, R = constant
    ctdtemp_c, eqtemp_c = temp_c
    totalvol_ml, airvol_ml = volumen_ml
    ce_mol, ce_pa = ppm_to_mol(ce_ppm, airpres_hpa, totalvol_ml, airvol_ml, eqtemp_c, R)
    ca_mol, _ = ppm_to_mol(ca_ppm, airpres_hpa, totalvol_ml, airvol_ml, eqtemp_c, R)
    opt_bottle = minimize(opt_fuc_bottle, pH, args=(ce_pa, ce_mol, eqtemp_c, alk, hcp25, dlnHcpd1_T, airvol_ml), method='Nelder-Mead')
    _, dic_mol, _, _ = carbon_balance(opt_bottle.x[0], ce_mol, eqtemp_c, alk, hcp25, dlnHcpd1_T, airvol_ml)
    opt_insitu = minimize(opt_fuc_insitu, pH, args=(dic_mol, ca_mol, ce_mol, ctdtemp_c, alk, hcp25, dlnHcpd1_T, airvol_ml), method='Nelder-Mead')
    _, _, _, co2_molm3 = carbon_balance(opt_insitu.x[0], ce_mol, ctdtemp_c, alk, hcp25, dlnHcpd1_T, airvol_ml)
    return co2_molm3*1000

def main():
    ce_ppm = 300
    ca_ppm = 402.9
    airpres_hpa = 1036
    totalvol_ml = 1160
    airvol_ml = 600
    eqtemp_c = 7
    ctdtemp_c = 5
    R = 8.314
    dlnHcpd1_T = 2392.86
    hcp25 = 3.3E-4
    alk = np.arange(3,4,.1)
    pH = 8.1

    constant = (hcp25, dlnHcpd1_T, R)
    co2_mmolm3 = []
    for alk_i in alk:
        co2_mmolm3.append(alkalinity_correction((ce_ppm, ca_ppm), alk_i, pH, airpres_hpa, (totalvol_ml, airvol_ml), (ctdtemp_c, eqtemp_c), constant))
    co2_mmolm3 = np.array(co2_mmolm3)
    ce_mol, ce_pa = ppm_to_mol(ce_ppm, airpres_hpa, totalvol_ml, airvol_ml, eqtemp_c, R)
    ca_mol, _ = ppm_to_mol(ca_ppm, airpres_hpa, totalvol_ml, airvol_ml, eqtemp_c, R)
    c_umolL = equibrium(ce_mol, ce_pa, ca_mol, eqtemp_c, hcp25, dlnHcpd1_T, airvol_ml)
    plt.plot((alk[round(len(alk)/2)]-alk)*50, (co2_mmolm3[round(len(alk)/2)]-co2_mmolm3)/co2_mmolm3[round(len(alk)/2)]*1000)
    plt.show()
    #print(c_umolL)



if __name__ == "__main__":
    main()
