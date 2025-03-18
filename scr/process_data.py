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
    hcp_molm3Pa = hcp(
        varname, data["Te_degC"], data["sal_psu"], constant, cfg["H_coeff"]
    )
    # mols in the water after shaking
    caq_mol = hcp_molm3Pa * ce_pa * (cfg["tot_vol"] - cfg["hs_vol"]) / 1000 / 1000
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


def hcp(varname, temp_c, sal_psu, constant, coeff=1):
    """Calculate Henry's Coefficient

    Parameters
    ----------

    varname : string
        Could be ch4 or co2
    temp_c : pandas.Series
        Temperature of equilibrium.
    sal_psu : pandas.Series
        Water salinity in PSU
    constant : dict
        Constant data (defined in utils/constant.yml).
    coeff : float
        0 Use Sanders 2015 - Only consider correction by temperature (Use this for freshwaters).
        1 Use Wiss 1974 to estimate Henry's Coefficient for CO2 and Wiesenburg & Guinasso 1979 for CH4, considering Salinity and Temperature coerrections (default).

    Returns
    -------
    hcp_t : float
            Henry's coefficient for CH4 or CO2 in mmolm-3Pa-1
    """

    if coeff == 0:
        if varname == "ch4":
            hcp25 = constant["H_CH4_T25"]
            dlnHcpd1_T = constant["dlnH_CH4_d1T"]
        elif varname == "co2":
            hcp25 = constant["H_CO2_T25"]
            dlnHcpd1_T = constant["dlnH_CO2_d1T"]
        hcp_t = hcp25 * np.exp(dlnHcpd1_T * (1 / (temp_c + 273.15) - 1 / 298.15))
        return hcp_t

    elif coeff == 1:
        hcp_t = hcp_sal(varname, sal_psu, temp_c)
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
    dce_permil = data[variables[2]]
    ce_mol, _ = ppm_to_mol(cfg, ce_ppm, data["Te_degC"], data["AirP_hPa"], constant)
    e13c_mol, e12c_mol = c13(constant, ce_mol, dce_permil)

    ca_ppm = data[variables[1]]
    dca_permil = data[variables[3]]
    ca_mol, _ = ppm_to_mol(cfg, ca_ppm, data["AirT_degC"], data["AirP_hPa"], constant)
    a13c_mol, a12c_mol = c13(constant, ca_mol, dca_permil)

    aq13c12c = constant["RVDPD"] * (data[variables[2]] / 1000 + 1)
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


def hcp_sal(var, sal_psu, temp_c):
    """Correct Henrys coefficient for CO2 (Weiss 1974) and CH4 (Wiesenburg and Guinasso 1979) considering salnity and temperature corrections.

    Parameters
    ----------
    var : string.
        Variable name (ch4, co2).
    sal_psu : pandas.Series or float.
        Salinity in PSU.
    temp_c : pandas.Series or float.
        Water temperature in deg C.

    Returns
    -------
    hcpsalt_mmolm3Pa : list
        Henry's coefficient estimated considering salinity and temperature corrections [molm-3Pa-1]
    """
    temp_k = temp_c + 273.15
    if var == "co2":  # folowing Weiss 1974
        A = [-58.0931, 90.5069, 22.2940]
        B = [0.027766, -0.025888, 0.0050578]
        hcpsalt_molm3Pa = (
            np.exp(
                A[0]
                + A[1] * 100 / temp_k
                + A[2] * np.log(temp_k / 100)
                + sal_psu * (B[0] + B[1] * temp_k / 100 + B[2] * (temp_k / 100) ** 2)
            )
            * 1000
            / 101325
        )
    elif var == "ch4":  # folowing Wiesenburg and Guinasso 1979
        A = [-417.5053, 599.8626, 380.3636, -62.0764]
        B = [-0.064236, 0.03498, -0.0052732]
        c_molm3 = (
            np.exp(
                +A[0]
                + A[1] * 100 / temp_k
                + A[2] * np.log(temp_k / 100)
                + A[3] * (temp_k / 100)
                + sal_psu * (B[0] + B[1] * temp_k / 100 + B[2] * (temp_k / 100) ** 2)
            )
            * 1000
            / 1e9
        )
        hcpsalt_molm3Pa = c_molm3 / 101325

    return hcpsalt_molm3Pa


def hcpch4_sal(hcp_molLPa, sal_psu, temp_c):
    """Correct Henrys coefficient for CH4 by salinity (Lee et al., 2020).
    NOT USED

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
    # A1 = -68.8862
    # A2 = 101.4956
    # A3 = 28.7314
    # B1 = -0.076146
    # B2 = 0.043970
    # B3 = -0.0068612

    temp_k = temp_c + 273.15
    # lnb = A1 + A2*(100/temp_k) + A3*np.log(temp_k/100) + sal_psu/100*(B1 + B2*temp_k/100 + (B3*temp_k/100)**2)

    ksalt = (
        3.38828
        - 0.0318765 * temp_k
        + 1.22003e-4 * temp_k**2
        - 2.31891e-7 * temp_k**3
        + 2.22938e-10 * temp_k**4
        - 8.83764e-14 * temp_k**5
    )
    sr_gkg = 35.16504 / 35 * sal_psu
    m = 1000 / 31.4038218 * (sr_gkg / (1 - sr_gkg))
    hcpsalt_molLPa = hcp_molLPa * np.exp(ksalt * m)

    return hcpsalt_molLPa


def hcpco2_sal(hcp_molLPa, sal_psu, temp_c):
    """Correct Henrys coefficient for CO2 by salinity (Lee et al., 2020).
    NOT USED

    Parameters
    ----------
    hcp_molLPa : float
        Henry's coefficient for insitu temperature in molL-1Pa-1
    sal_psu : pandas.Series
        Water salinity in PSU
    temp_c : pandas.Series
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
