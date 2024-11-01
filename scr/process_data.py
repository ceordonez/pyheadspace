
import numpy as np
from scr.read_data import read_constant

def process_data(cfg, data):
    """Process data.

    Parameters
    ----------
    cfg : TODO
    data : TODO

    Returns
    -------
    TODO

    """
    constant = read_constant('utils/constant.yml')
    for varname in cfg['varnames']:
        variables = cfg['varnames'][varname]
        colname = '_'.join([varname, 'mmolm3'])
        c_mmolm3, d13c_permil = equilibrium(cfg, varname, data, variables, constant)
        data[colname] = c_mmolm3
        colname = '_'.join(['d13', varname, 'permil'])
        data[colname] = d13c_permil
    return data


def ppm_to_mol(cfg, c_ppm, temp_c, data, constant):
    """Transform ppm concentration to mol.

    Parameters
    ----------
    cfg : TODO
    c_ppm : TODO
    temp_c : TODO
    data : TODO

    Returns
    -------
    TODO

    """
    c_pa = c_ppm/10000*data['AirP_hPa']
    c_mol = c_pa*(cfg['hs_vol'])/1000/1000/(constant['R']*(temp_c + 273.15))
    return c_mol, c_pa

def hcp(varname, temp_c, sal_psu, hcp25, dlnHcpd1_T):
    hcp_t = hcp25*np.exp(dlnHcpd1_T*(1/(temp_c + 273.15)-1/298.15))/1000
    if varname == 'ch4':
        hcp_t = hcpch4_sal(hcp_t, sal_psu, temp_c)
    elif varname == 'co2':
        hcp_t = hcpco2_sal(hcp_t, sal_psu, temp_c)
    return hcp_t

def equilibrium(cfg, varname, data, variables, constant):
    """TODO: Docstring for equilibrium.

    Parameters
    ----------
    cfg : TODO
    data : TODO
    variables : TODO
    constant : TODO

    Returns
    -------
    TODO

    """
    c_ppm = data[variables[0]]
    ce_mol, ce_pa = ppm_to_mol(cfg, c_ppm, data['Te_degC'], data, constant)
    catm_ppm = data[variables[1]]
    ca_mol, _ = ppm_to_mol(cfg, catm_ppm, data['AirT_degC'], data, constant)
    hcp_molLPa = hcp(varname, data['Te_degC'], data['sal_psu'], constant[variables[2]], constant[variables[3]])
    caq_mol = hcp_molLPa*ce_pa*(cfg['tot_vol'] - cfg['hs_vol'])/1000
    ntotal = ce_mol + caq_mol
    c_molm3 = (ntotal - ca_mol)/((cfg['tot_vol'] - cfg['hs_vol'])/1000/1000)
    c_mmolm3 = c_molm3*1000
    d13c_permil = dC_equilibrium(cfg, data, variables, constant, caq_mol)
    return c_mmolm3, d13c_permil

def dC_equilibrium(cfg, data, variables, constant, caq_mol):
    """TODO: Docstring for equilibrium.

    Parameters
    ----------
    cfg : TODO
    data : TODO
    variables : TODO
    constant : TODO

    Returns
    -------
    TODO

    """
    ce_ppm = data[variables[0]]
    dce_permil = data[variables[4]]
    ce_mol, _ = ppm_to_mol(cfg, ce_ppm, data['Te_degC'], data, constant)
    e13c_mol, e12c_mol = c13(constant, ce_mol, dce_permil)

    ca_ppm = data[variables[1]]
    dca_permil = data[variables[5]]
    ca_mol, _ =ppm_to_mol(cfg, ca_ppm, data['AirT_degC'], data, constant)
    a13c_mol, a12c_mol = c13(constant, ca_mol, dca_permil)

    aq13c12c = constant['RVDPD']*(data[variables[4]]/1000+1)
    aq12c_mol = caq_mol/(1 + aq13c12c)
    aq13c_mol = caq_mol - aq12c_mol
    tot13c_mol = e13c_mol + aq13c_mol
    tot12c_mol = e12c_mol + aq12c_mol
    w13c_mol = tot13c_mol - a13c_mol
    w12c_mol = tot12c_mol - a12c_mol
    w13c12c = w13c_mol/w12c_mol
    dc_permil = (w13c12c/constant['RVDPD']-1)*1000
    return dc_permil



def c13(constant, c_mol, dc_permil):
    """TODO: Docstring for c13.

    Parameters
    ----------
    cfg : TODO
    variables : TODO
    constant : TODO
    data : TODO

    Returns
    -------
    TODO

    """
    e13c12c = constant['RVDPD']*(dc_permil/1000+1)
    e12c = c_mol/(1 + e13c12c)
    e13c = c_mol - e12c
    return e13c, e12c

def hcpco2_sal(hcp_molLPa, sal_psu, temp_c):
    """Correct Henrys coefficient for CO2 by salinity (Lee et al., 2020).

    Parameters
    ----------
    hcp : TODO
    sal_psu : TODO
    temp_c : TODO

    Returns
    -------
    TODO
    """
    temp_k = temp_c + 273.15
    ksalt = 0.11572 - 6.0293E-4*temp_k +3.5817E-6*temp_k**2 - 3.772E-9*temp_k**3
    hcpsalt_molLPa = hcp_molLPa*np.exp(ksalt*sal_psu*1E-3)
    return hcpsalt_molLPa

def hcpch4_sal(hcp_molLPa, sal_psu, temp_c):
    """Correct Henrys coefficient for CH4 by salinity (Lee et al., 2020).

    Parameters
    ----------
    hcp : TODO
    sal_psu : TODO
    temp_c : TODO

    Returns
    -------
    TODO
    """
    #A1 = -68.8862
    #A2 = 101.4956
    #A3 = 28.7314
    #B1 = -0.076146
    #B2 = 0.043970
    #B3 = -0.0068612

    #temp_k = temp_c + 273.15
    #lnb = A1 + A2*(100/temp_k) + A3*np.log(temp_k/100) + sal_psu/100*(B1 + B2*temp_k/100 + (B3*temp_k/100)**2)

    ksalt = 3.38828 - 0.0318765*temp_c + 1.22003E-4*temp_c**2 - 2.31891E-7*temp_c**3 +2.22938E-10*temp_c**4 - 8.83764E-14*temp_c**5
    hcpsalt_molLPa = hcp_molLPa*np.exp(ksalt*sal_psu*1E-3)

    return hcpsalt_molLPa


