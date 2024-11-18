import numpy as np

from scr.read_data import read_constant


def process_data(cfg, data):
    """Process data.

    Parameters
    ----------
    cfg : dict
        Configuration file
    data : pandas.DataFrame
        Input data.

    Returns
    -------
    data : pandas.Dataframe
        Results including CO2 and CH4 concentrations in mmolm3 and isotopic corrected by headspace method.

    """
    constant = read_constant("utils/constant.yml")
    for varname in cfg["varnames"]:
        variables = cfg["varnames"][varname]
        colname = "_".join([varname, "mmolm3"])
        c_mmolm3, d13c_permil = equilibrium(cfg, varname, data, variables, constant)
        data[colname] = c_mmolm3
        colname = "_".join(["d13", varname, "permil"])
        data[colname] = d13c_permil
    return data


def equilibrium(cfg, varname, data, variables, constant):
    """Make equilibrium for the headspace method.

    Parameters
    ----------
    cfg : dict
        Configuration file.
    data : pandas.DataFrame
        Input data.
    variables : list
        Variables to be used to calculate equilibrium (defined in cfg file).
    constant : dict
        Constant data from utils/constant.yml.

    Returns
    -------
    c_mmolm3, d13c_permil : pandas.Series, pandas.Series
        c_mmolm3 is concentration in the water in mmolm3
        d13c_permil is the isotopic composition in the water in permille

    """
    c_ppm = data[variables[0]]
    # mols in the headspace after shaking
    ce_mol, ce_pa = ppm_to_mol(cfg, c_ppm, data["Te_degC"], data["AirP_hPa"], constant)
    # mols in the heaspace before shaking
    catm_ppm = data[variables[1]]
    ca_mol, _ = ppm_to_mol(cfg, catm_ppm, data["AirT_degC"], data["AirP_hPa"], constant)
    # Estimate Henry's constant
    hcp_molLPa = hcp(
        varname,
        data["Te_degC"],
        data["sal_psu"],
        constant[variables[2]],
        constant[variables[3]],
    )
    # mols in the water after shaking
    caq_mol = hcp_molLPa * ce_pa * (cfg["tot_vol"] - cfg["hs_vol"]) / 1000
    ntotal = ce_mol + caq_mol
    # mols in the water before shaking
    c_molm3 = (ntotal - ca_mol) / ((cfg["tot_vol"] - cfg["hs_vol"]) / 1000 / 1000)
    c_mmolm3 = c_molm3 * 1000
    d13c_permil = dC_equilibrium(cfg, data, variables, constant, caq_mol)
    return c_mmolm3, d13c_permil


def ppm_to_mol(cfg, c_ppm, temp_c, airP_hpa, constant):
    """Transform ppm concentration to mol.

    Parameters
    ----------
    cfg : dict
        Configuration information.
    c_ppm : pandas.Series
        Gas concentration in the headspace in ppm.
    temp_c : pandas.Series
        Temperature of equilibrium.
    airP_hpa : pandas.Series
        Air pressure in hPa
    constant : dict
        Constant information from utils/constant.yml.

    Returns
    -------
    c_mol, c_pa : pd.Series
        Gas concentration in the headspace in mol and Pa.

    """
    c_pa = c_ppm / 10000 * airP_hpa
    c_mol = c_pa * (cfg["hs_vol"]) / 1000 / 1000 / (constant["R"] * (temp_c + 273.15))
    return c_mol, c_pa


def hcp(varname, temp_c, sal_psu, hcp25, dlnHcpd1_T):
    """Calculate Henry's constant corrected by salinity based on XX

    Parameters
    ----------
    varname : string
        Could be ch4 or co2
    temp_c : pandas.Series
        Temperature of equilibrium.
    sal_psu : pandas.Series
        Water salinity in PSU
    hcp25 : float
        Henry's constant at 25 degree C (defined in utils/constant.yml)
    dlnHcpd1_T : float
        parameters to transform Henry's coefficient from 25 to any other temperature (see XX).

    Returns
    -------
    hcp_t : pandas.Series
        Henry's coefficient for temperature of equilibrium.
    """
    hcp_t = hcp25 * np.exp(dlnHcpd1_T * (1 / (temp_c + 273.15) - 1 / 298.15)) / 1000
    if varname == "ch4":
        hcp_t = hcpch4_sal(hcp_t, sal_psu, temp_c)
    elif varname == "co2":
        hcp_t = hcpco2_sal(hcp_t, sal_psu, temp_c)
    return hcp_t


def dC_equilibrium(cfg, data, variables, constant, caq_mol):
    """Calcualte isotopic mass balance for the headspace method.

    Parameters
    ----------
    cfg : dict
        Configuration file.
    data : pandas.DataFrame
        Input data.
    variables : list
        Variables to be used to calculate equilibrium (defined in cfg file).
    constant : dict
        Constant data (defined in utils/constant.yml).
    caq_mol : pandas.Series
        Mols of gas in the water after shaking.

    Returns
    -------
    dc_permil : pandas.Series
        Isotopic signature in the water before shaking.

    """
    ce_ppm = data[variables[0]]
    dce_permil = data[variables[4]]
    ce_mol, _ = ppm_to_mol(cfg, ce_ppm, data["Te_degC"], data["AirP_hPa"], constant)
    e13c_mol, e12c_mol = c13(constant, ce_mol, dce_permil)

    ca_ppm = data[variables[1]]
    dca_permil = data[variables[5]]
    ca_mol, _ = ppm_to_mol(cfg, ca_ppm, data["AirT_degC"], data["AirP_hPa"], constant)
    a13c_mol, a12c_mol = c13(constant, ca_mol, dca_permil)

    aq13c12c = constant["RVDPD"] * (data[variables[4]] / 1000 + 1)
    aq12c_mol = caq_mol / (1 + aq13c12c)
    aq13c_mol = caq_mol - aq12c_mol
    tot13c_mol = e13c_mol + aq13c_mol
    tot12c_mol = e12c_mol + aq12c_mol
    w13c_mol = tot13c_mol - a13c_mol
    w12c_mol = tot12c_mol - a12c_mol
    w13c12c = w13c_mol / w12c_mol
    dc_permil = (w13c12c / constant["RVDPD"] - 1) * 1000
    return dc_permil


def c13(constant, c_mol, dc_permil):
    """Calculates isotopes concentration.

    Parameters
    ----------
    constant : dict
        Constats data (see utils/constat.yml).
    c_mol: pandas.Series
        Air concetation after shaking in mol.
    dc_permil : pandas.Series
        Air isotopic signature after shaking in permille.

    Returns
    -------
    e13c, e12c : pandas.Series, pandas.Series
        e13c is the concentration of C13 in mol after equilibrim.
        e12c is the concentration of C12 in mol after equilibrim.

    """
    e13c12c = constant["RVDPD"] * (dc_permil / 1000 + 1)
    e12c = c_mol / (1 + e13c12c)
    e13c = c_mol - e12c
    return e13c, e12c


def hcpco2_sal(hcp_molLPa, sal_psu, temp_c):
    """Correct Henrys coefficient for CO2 by salinity (Lee et al., 2020).

    Parameters
    ----------
    hcp_molLPa : float
        Henry's coefficient for insitu temperature in molL-1Pa-1
    sal_psu : pandas.Series
        Water salinity in PSU
    temp_c : TODO
        Equilibrium water temperature in deg C.

    Returns
    -------
    hcpsalt_molLPa : float
        Henry's coefficient corrected by salinity
    """
    temp_k = temp_c + 273.15
    ksalt = 0.11572 - 6.0293e-4 * temp_k + 3.5817e-6 * temp_k**2 - 3.772e-9 * temp_k**3
    hcpsalt_molLPa = hcp_molLPa * np.exp(ksalt * sal_psu * 1e-3)
    return hcpsalt_molLPa


def hcpch4_sal(hcp_molLPa, sal_psu, temp_c):
    """Correct Henrys coefficient for CH4 by salinity (Lee et al., 2020).

    Parameters
    ----------
    hcp_molLPa : float
        Henry's coefficient for insitu temperature in molL-1Pa-1
    sal_psu : pandas.Series
        Water salinity in PSU
    temp_c : TODO
        Equilibrium water temperature in deg C.

    Returns
    -------
    hcpsalt_molLPa : float
        Henry's coefficient corrected by salinity

    """
    ksalt = (
        3.38828
        - 0.0318765 * temp_c
        + 1.22003e-4 * temp_c**2
        - 2.31891e-7 * temp_c**3
        + 2.22938e-10 * temp_c**4
        - 8.83764e-14 * temp_c**5
    )
    hcpsalt_molLPa = hcp_molLPa * np.exp(ksalt * sal_psu * 1e-3)

    return hcpsalt_molLPa
