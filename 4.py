import streamlit as st
from scipy.integrate import quad       # type: ignore
from scipy.integrate import odeint     # type: ignore
from scipy.optimize import fsolve      # type: ignore
import matplotlib.pyplot as plt        # type: ignore
import numpy as np

# ── Physics helpers ────────────────────────────────────────────────────────────

def cs(H_):
    return 1.005 + 1.884 * H_

def iG_H_TG(H_, TG_):
    return cs(H_) * TG_ + 2502 * H_

def iG(TL_, param):
    kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal = param
    return iG1 + (L * cL / G) * (TL_ - TL1)

def TL(iG_, param):
    kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal = param
    return TL1 + ((iG_ - iG1) * G) / (L * cL)

def iGi_aux(Ti_):
    if Ti_ < 7.6562545:
        return (9.36 + 1.613 * Ti_) / (1 - 0.01265 * Ti_ + 6.0e-5 * Ti_**2)
    elif Ti_ > 50.390373:
        return (2.28 * Ti_ - 16.9) / (1 - 0.015862 * Ti_ + 5.88e-5 * Ti_**2)
    else:
        return (10.42 + 1.37 * Ti_) / (1 - 0.019 * Ti_ + 9.5e-5 * Ti_**2)

def iGi(Ti_):
    if type(Ti_) is list:
        return [iGi_aux(a) for a in Ti_]
    elif type(Ti_) is np.ndarray:
        return np.array([iGi_aux(a) for a in Ti_])
    return iGi_aux(Ti_)

def pendiGi(Ti_):
    h = 1e-6
    return (iGi_aux(Ti_ + h) - iGi_aux(Ti_ - h)) / (2 * h)

def subintv(iG_, param):
    return 1 / (iGi(TL(iG_, param)) - iG_)

def subinth(TL_, param):
    iGsub = iG(TL_, param)
    Ti = fsolve(lambda T: iGsub - iGi(T), TL_)[0]
    return 1 / (TL_ - Ti)

def subintd(iG_, param):
    kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal = param
    pend = -hLa / kYa
    TLsub = TL(iG_, param)
    Ti = fsolve(lambda T: (iG_ - iGi(T)) - pend * (TLsub - T), TLsub)[0]
    return 1 / (iGi(Ti) - iG_)

def z_v(param):
    kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal = param
    iG2 = iG(TL2, param)
    return (G / kYa) * quad(subintv, iG1, iG2, args=param)[0]

def z_v_g_aux(iG_, param):
    kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal = param
    return (G / kYa) * quad(subintv, iG1, iG_, args=param)[0]

def z_v_g(iG_, param):
    if type(iG_) is list:
        return [z_v_g_aux(a, param) for a in iG_]
    elif type(iG_) is np.ndarray:
        return np.array([z_v_g_aux(a, param) for a in iG_])
    return z_v_g_aux(iG_, param)

def z_h(param):
    kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal = param
    return (L * cL / hLa) * quad(subinth, TL1, TL2, args=param)[0]

def z_h_g_aux(TL_, param):
    kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal = param
    return (L * cL / hLa) * quad(subinth, TL1, TL_, args=param)[0]

def z_h_g(TL_, param):
    if type(TL_) is list:
        return [z_h_g_aux(a, param) for a in TL_]
    elif type(TL_) is np.ndarray:
        return np.array([z_h_g_aux(a, param) for a in TL_])
    return z_h_g_aux(TL_, param)

def z_d(param):
    kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal = param
    iG2 = iG(TL2, param)
    return (G / kYa) * quad(subintd, iG1, iG2, args=param)[0]

def z_d_g_aux(iG_, param):
    kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal = param
    return (G / kYa) * quad(subintd, iG1, iG_, args=param)[0]

def z_d_g(iG_, param):
    if type(iG_) is list:
        return [z_d_g_aux(a, param) for a in iG_]
    elif type(iG_) is np.ndarray:
        return np.array([z_d_g_aux(a, param) for a in iG_])
    return z_d_g_aux(iG_, param)

def perfils_H_TG(H_TG_, z_, Hi_, Ti_, param):
    kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal = param
    dHdz  = (Hi_ - H_TG_[0]) / (G / kYa)
    dTGdz = (Ti_ - H_TG_[1]) / (G / kYa)
    return [dHdz, dTGdz]

# ── Streamlit UI ───────────────────────────────────────────────────────────────

st.set_page_config(page_title="Cooling / Dehumidification Tower", layout="wide")
st.title("Water Cooling or Air Dehumidification Tower")

with st.sidebar:
    st.header("Parameters")
    TL2  = st.number_input("TL2 [ºC]",              value=40.0, format="%.4f")
    TL1  = st.number_input("TL1 [ºC]",              value=25.0, format="%.4f")
    TG1  = st.number_input("TG1 [ºC]",              value=20.0, format="%.4f")
    H1   = st.number_input("H1 [kg/kg]",            value=0.01, format="%.6f")
    L    = st.number_input("L/S [kg/(s·m²)]",       value=1.0,  format="%.4f")
    G    = st.number_input("G'/S [kg/(s·m²)]",      value=1.0,  format="%.4f")
    kYa  = st.number_input("kYa [kg/(s·m³)]",       value=0.5,  format="%.4f")
    hLa  = st.number_input("hLa [kJ/(s·ºC·m³)]",   value=0.5,  format="%.4f")
    solve = st.button("Solve", type="primary", use_container_width=True)

result_col, graph_col = st.columns([1, 3])

if solve:
    # ── Input validation ───────────────────────────────────────────────────────
    error = None
    if TL1 >= 70 or TL2 >= 70:
        error = "Water temperature must be below 70 ºC."
    elif TL1 <= 5 or TL2 <= 5:
        error = "Water temperature must be above 5 ºC."
    elif TL2 == TL1:
        error = "TL1 and TL2 can't be the same."
    elif TG1 < -15:
        error = "The temperature of air shouldn't be that low."
    elif TG1 > 90:
        error = "The temperature of air shouldn't be that high."
    elif H1 < 0:
        error = "Humidity can't be negative."
    elif H1 > 1.42:
        error = "Humidity must be lower."
    elif L > 10:
        error = "L/S must be lower."
    elif L < 0.05:
        error = "L/S must be higher."
    elif G > 10:
        error = "G'/S must be lower."
    elif G < 0.05:
        error = "G'/S must be higher."
    elif kYa < 0.05:
        error = "kYa must be higher."
    elif hLa < 0.05:
        error = "hLa must be higher."

    if error:
        st.error(error)
        st.stop()

    # ── Derived quantities ─────────────────────────────────────────────────────
    if kYa / hLa > 100000:
        kYa = float("inf")
    if hLa / kYa > 100000:
        hLa = float("inf")

    cL  = 4.18
    iG1 = iG_H_TG(H1, TG1)
    zreal = None
    param1 = [kYa, hLa, cL, L, G, TL1, TL2, iG1, zreal]
    pend1  = -hLa / kYa

    if iG1 > iGi(TG1):
        st.error("iG1 can't be above the saturation point.")
        st.stop()
    if TL2 > TL1 and iG1 >= iGi(TL1):
        st.error(f"For this humidification iG1 is too high, iG1max = {round(iGi_aux(TL1), 2)} kJ/kg.")
        st.stop()
    if TL2 < TL1 and iG1 <= iGi(TL1):
        st.error(f"For this dehumidification iG1 is too low, iG1min = {round(iGi_aux(TL1), 2)} kJ/kg.")
        st.stop()
    if hLa == float("inf") and kYa == float("inf"):
        st.error("hLa and kYa can't be both infinite.")
        st.stop()

    # ── Intersection check ─────────────────────────────────────────────────────
    TLintersection, infodict, ier, msg = fsolve(
        lambda TLv: iGi(TLv) - iG(TLv, param1), x0=TL1, full_output=True
    )

    if ier == 1 and (TL1 - TLintersection) * (TL2 - TLintersection) < 0:
        if TL2 > TL1:
            pendromax = (iGi(TL2) - iG1) / (TL2 - TL1)
            if pendromax < pendiGi(TL2):
                TLpendiguals = fsolve(
                    lambda TLv: pendiGi(TLv) - (iGi(TLv) - iG1) / (TLv - TL1), TL2
                )[0]
                pendromax = (iGi(TLpendiguals) - iG1) / (TLpendiguals - TL1)
            else:
                pendromax = (iG1 - iGi(TL2)) / (TL1 - TL2)
            Gmin = L * cL / pendromax
            Lmax = pendromax * G / cL
            st.warning(
                f"The operation line intersects the saturation curve. "
                f"(G'/s)min = {round(Gmin, 2)} kg/(h·m²) or (L/S)max = {round(Lmax, 2)} kg/(h·m²)."
            )
            fig_warn, ax_warn = plt.subplots(figsize=(7, 5))
            x = np.linspace(TL1, TL2, 1000)
            ax_warn.plot(x, iGi(x), color="black", label="Saturation curve")
            ax_warn.plot(x, iG(x, param1), color="blue", label="Operation line")
            ax_warn.legend()
            ax_warn.set_xlabel("Temperature [ºC]")
            ax_warn.set_ylabel("iG [kJ/kg]")
            st.pyplot(fig_warn)
            plt.close(fig_warn)
            st.stop()

    # ── Solve z ────────────────────────────────────────────────────────────────
    if pend1 == -float("inf"):
        z = z_v(param1)
    elif pend1 == 0:
        z = z_h(param1)
    else:
        z = z_d(param1)

    with result_col:
        st.subheader("Result")
        st.metric(label="Tower height z [m]", value=f"{round(z, 4)}")

    # ── Tie-line endpoints ─────────────────────────────────────────────────────
    Ti1 = TL1
    Ti5 = TL2

    if pend1 == -float("inf"):
        xru1 = (TL1, TL1);                         yru1 = (iG(TL1, param1),                iGi(TL1))
        xru2 = (TL1+(TL2-TL1)/4,)*2;               yru2 = (iG(TL1+(TL2-TL1)/4, param1),   iGi(TL1+(TL2-TL1)/4))
        xru3 = (TL1+2*(TL2-TL1)/4,)*2;             yru3 = (iG(TL1+2*(TL2-TL1)/4, param1), iGi(TL1+2*(TL2-TL1)/4))
        xru4 = (TL1+3*(TL2-TL1)/4,)*2;             yru4 = (iG(TL1+3*(TL2-TL1)/4, param1), iGi(TL1+3*(TL2-TL1)/4))
        xru5 = (TL2, TL2);                         yru5 = (iG(TL2, param1),                iGi(TL2))
    elif pend1 == 0:
        Ti1 = fsolve(lambda T: iG(TL1, param1) - iGi(T), TL1)[0]
        xru1 = (TL1, Ti1);                         yru1 = (iG(TL1, param1), iGi(Ti1))
        Ti2 = fsolve(lambda T: iG(TL1+(TL2-TL1)/4, param1) - iGi(T), TL1+(TL2-TL1)/4)[0]
        xru2 = (TL1+(TL2-TL1)/4, Ti2);             yru2 = (iG(TL1+(TL2-TL1)/4, param1), iGi(Ti2))
        Ti3 = fsolve(lambda T: iG(TL1+2*(TL2-TL1)/4, param1) - iGi(T), TL1+2*(TL2-TL1)/4)[0]
        xru3 = (TL1+2*(TL2-TL1)/4, Ti3);           yru3 = (iG(TL1+2*(TL2-TL1)/4, param1), iGi(Ti3))
        Ti4 = fsolve(lambda T: iG(TL1+3*(TL2-TL1)/4, param1) - iGi(T), TL1+3*(TL2-TL1)/4)[0]
        xru4 = (TL1+3*(TL2-TL1)/4, Ti4);           yru4 = (iG(TL1+3*(TL2-TL1)/4, param1), iGi(Ti4))
        Ti5 = fsolve(lambda T: iG(TL2, param1) - iGi(T), TL2)[0]
        xru5 = (TL2, Ti5);                         yru5 = (iG(TL2, param1), iGi(Ti5))
    else:
        Ti1 = fsolve(lambda T: (iG1 - iGi(T)) - pend1*(TL1 - T), TL1)[0]
        xru1 = (TL1, Ti1);                         yru1 = (iG(TL1, param1), iGi(Ti1))
        Ti2 = fsolve(lambda T: (iG(TL1+(TL2-TL1)/4, param1) - iGi(T)) - pend1*(TL1+(TL2-TL1)/4 - T), TL1+(TL2-TL1)/4)[0]
        xru2 = (TL1+(TL2-TL1)/4, Ti2);             yru2 = (iG(TL1+(TL2-TL1)/4, param1), iGi(Ti2))
        Ti3 = fsolve(lambda T: (iG(TL1+2*(TL2-TL1)/4, param1) - iGi(T)) - pend1*(TL1+2*(TL2-TL1)/4 - T), TL1+2*(TL2-TL1)/4)[0]
        xru3 = (TL1+2*(TL2-TL1)/4, Ti3);           yru3 = (iG(TL1+2*(TL2-TL1)/4, param1), iGi(Ti3))
        Ti4 = fsolve(lambda T: (iG(TL1+3*(TL2-TL1)/4, param1) - iGi(T)) - pend1*(TL1+3*(TL2-TL1)/4 - T), TL1+3*(TL2-TL1)/4)[0]
        xru4 = (TL1+3*(TL2-TL1)/4, Ti4);           yru4 = (iG(TL1+3*(TL2-TL1)/4, param1), iGi(Ti4))
        Ti5 = fsolve(lambda T: (iG(TL2, param1) - iGi(T)) - pend1*(TL2 - T), TL2)[0]
        xru5 = (TL2, Ti5);                         yru5 = (iG(TL2, param1), iGi(Ti5))

    xeq = np.linspace(TL1, TL2, 1000) if pend1 == -float("inf") else np.linspace(Ti1, Ti5, 1000)
    xro = np.linspace(TL1, TL2, 1000)

    # ── Profile arrays ─────────────────────────────────────────────────────────
    if pend1 == -float("inf"):
        y_iG = np.linspace(iG1, iG(TL2, param1), 1000)
        x_z  = z_v_g(y_iG, param1)
        y_TL = TL(y_iG, param1)
    elif pend1 == 0:
        y_TL = np.linspace(TL1, TL2, 1000)
        x_z  = z_h_g(y_TL, param1)
        y_iG = iG(y_TL, param1)
    else:
        y_iG = np.linspace(iG1, iG(TL2, param1), 1000)
        x_z  = z_d_g(y_iG, param1)
        y_TL = TL(y_iG, param1)

    # ── H / TG profiles ────────────────────────────────────────────────────────
    H_TG0 = [H1, TG1]
    if pend1 == -float("inf"):
        Ti_prof = y_TL
    elif pend1 == 0:
        Ti_prof = [fsolve(lambda T: iGi(T) - a, TL1)[0] for a in y_iG]
    else:
        Ti_prof = [fsolve(lambda T: (a - iGi(T)) - pend1*(b - T), TG1)[0] for a, b in zip(y_iG, y_TL)]

    Hi = [fsolve(lambda H: iGi(a) - iG_H_TG(H, a), H1)[0] for a in Ti_prof]

    if pend1 == 0:
        H_prof  = Hi
        TG_prof = Ti_prof
        x_z0 = (0, 0)
        y_H1_pts  = (H1, H_prof[0])
        y_TG1_pts = (TG1, TG_prof[0])
    else:
        H_prof  = np.empty_like(x_z)
        TG_prof = np.empty_like(x_z)
        H_prof[0]  = H_TG0[0]
        TG_prof[0] = H_TG0[1]
        for i in range(1, 1000):
            zspan = [x_z[i-1], x_z[i]]
            H_TG = odeint(perfils_H_TG, H_TG0, zspan, args=(Hi[i], Ti_prof[i], param1))
            H_prof[i]  = H_TG[1][0]
            TG_prof[i] = H_TG[1][1]
            H_TG0 = H_TG[1]
        x_z0 = (0, 0)
        y_H1_pts  = (H1, H1)
        y_TG1_pts = (TG1, TG1)

    # ── Plots ──────────────────────────────────────────────────────────────────
    with graph_col:
        st.subheader("Graphs")
        fig, ax = plt.subplots(2, 2, figsize=(11, 8))

        # Operation line
        ax[0, 0].plot(xeq, iGi(xeq), color="black", label="Saturation curve")
        ax[0, 0].plot(xro, iG(xro, param1), color="blue", label="Operation line")
        for xr, yr in [(xru1, yru1), (xru2, yru2), (xru3, yru3), (xru4, yru4), (xru5, yru5)]:
            ax[0, 0].plot(xr, yr, color="red", label="Tie line" if xr is xru1 else "")
        ax[0, 0].legend(fontsize="small")
        ax[0, 0].set_title("Operation line")
        ax[0, 0].set_xlabel("Temperature [ºC]")
        ax[0, 0].set_ylabel("iG [kJ/kg]")

        # TL profile
        ax[1, 0].plot(x_z, y_TL, color="red", label="TL")
        ax[1, 0].legend()
        ax[1, 0].set_title("TL profile")
        ax[1, 0].set_xlabel("z [m]")
        ax[1, 0].set_ylabel("TL [ºC]")

        # H profile
        ax[0, 1].plot(x_z, H_prof, color="blue", label="H")
        ax[0, 1].plot(x_z0, y_H1_pts, color="blue")
        ax[0, 1].legend()
        ax[0, 1].set_title("H profile")
        ax[0, 1].set_xlabel("z [m]")
        ax[0, 1].set_ylabel("H [kg/kg]")

        # TG profile
        ax[1, 1].plot(x_z, TG_prof, color="red", label="TG")
        ax[1, 1].plot(x_z0, y_TG1_pts, color="red")
        ax[1, 1].legend()
        ax[1, 1].set_title("TG profile")
        ax[1, 1].set_xlabel("z [m]")
        ax[1, 1].set_ylabel("TG [ºC]")

        plt.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
