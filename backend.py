from ctypes import sizeof
from msilib import type_string
from scipy import signal
import scipy
import numpy as np
import matplotlib.pyplot as plt
from metodosnum import bisection
from matplotlib import patches
from matplotlib.figure import Figure
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D


def get_min_order(filter_name, Wpass, Watt, Gp, Ga):
    if filter_name == "butter":
        ordfunc = signal.buttord
    elif filter_name == "cheby":
        ordfunc = signal.cheb1ord
    elif filter_name == "cheby2":
        ordfunc = signal.cheb2ord
    elif filter_name == "cauer":
        ordfunc = signal.ellipord
    else:
        return
    return ordfunc(Wpass, Watt, Gp, Ga, True)


def get_filter(filter_name, filter_type, N, Wn, Wpass, Watt, Gp, Ga, denorm):

    if filter_name == "butter":
        filterfunc = signal.butter
    elif filter_name == "cheby":
        filterfunc = lambda ord, w3db, type_str, analog: signal.cheby1(
            ord, -Gp, w3db, type_str, analog
        )
    elif filter_name == "cheby2":
        filterfunc = lambda ord, w3db, type_str, analog: signal.cheby2(
            ord, -Ga, w3db, type_str, analog
        )
    elif filter_name == "cauer":
        filterfunc = lambda ord, w3db, type_str, analog: signal.ellip(
            ord, -Gp, -Ga, w3db, type_str, analog
        )
    else:
        return

    def getMagnAtWx(order, W3db, Wx):
        num, den = filterfunc(order, W3db, filter_type, True)
        w, h = signal.freqs(num, den, [Wx])
        res = 20 * np.log10(abs(h[0]))
        return res

    b1, a1 = filterfunc(N, Wn, filter_type, True)

    # Busco dos frecuencias para podes hallar una que clave exacto Gp en Wp
    if filter_type == "lowpass":
        Wp1 = Wn
        while getMagnAtWx(N, Wp1, Wpass) > Gp:
            Wp1 *= 0.99
        Wp2 = Wn * 1.02
        Wn0 = 1.01 * bisection(
            lambda x: getMagnAtWx(N, x, Wpass) - (Gp), [Wp1, Wp2], 0.01
        )
    elif filter_type == "highpass":
        Wp1 = Wn
        while getMagnAtWx(N, Wp1, Wpass) > Gp:
            Wp1 *= 1.01
        Wp2 = Wn * 0.98
        Wn0 = 0.99 * bisection(
            lambda x: getMagnAtWx(N, x, Wpass) - (Gp), [Wp1, Wp2], 0.01
        )
    elif filter_type == "bandpass":
        Wn0 = []
        WaA1 = Wn[0]
        while getMagnAtWx(N, [WaA1, Wn[1]], Wpass[0]) > Gp:
            WaA1 *= 1.01
        WpA2 = Wn[0] * 0.98
        Wn0.append(
            1.01
            * bisection(
                lambda x: getMagnAtWx(N, [x, Wn[1]], Wpass[0]) - (Gp),
                [WaA1, WpA2],
                0.01,
            )
        )

        WpB1 = Wn[1]
        while getMagnAtWx(N, [Wn0[0], WpB1], Wpass[1]) > Gp:
            WpB1 *= 0.99
        WpB2 = Wn[1] * 1.02

        magB1 = getMagnAtWx(N, [Wn0[0], WpB1], Wpass[1])
        magB2 = getMagnAtWx(N, [Wn0[0], WpB2], Wpass[1])

        var = 0.99 * bisection(
            lambda x: getMagnAtWx(N, [Wn0[0], x], Wpass[1]) - (Gp), [WpB2, WpB1], 0.01
        )

        Wn0.append(var)
    elif filter_type == "bandstop":
        Wn0 = []

        W1 = Wn[0]
        W2 = Wn[1]

        if abs(Wpass[1] - Watt[1]) > abs(Wpass[0] - Watt[0]):

            var = scipy.optimize.minimize_scalar(
                lambda delta: (
                    (getMagnAtWx(N, [Wn[0] - delta, Wn[1] + delta], Wpass[0]) - Gp)
                ),
                bounds=(Wn[0] * 0.9, Wn[0] * 1.1),
            )
        else:

            var = scipy.optimize.minimize_scalar(
                lambda delta: (
                    (getMagnAtWx(N, [Wn[0] - delta, Wn[1] + delta], Wpass[1]) - Gp)
                ),
                bounds=(Wn[1] * 0.9, Wn[1] * 1.1),
            )

        Wn0 = [Wn[0] - var.x, Wn[0] + var.x]

    if filter_type == "lowpass":
        Wa1 = Wn
        while getMagnAtWx(N, Wa1, Watt) < Ga:
            Wa1 *= 1.01
        Wa2 = Wn * 0.98
        Wn100 = 0.99 * bisection(
            lambda x: getMagnAtWx(N, x, Watt) - (Ga), [Wa1, Wa2], 0.01
        )
    elif filter_type == "highpass":
        Wa1 = Wn
        while getMagnAtWx(N, Wa1, Watt) < Ga:
            Wa1 *= 0.99
        Wa2 = Wn * 1.02
        Wn100 = 1.01 * bisection(
            lambda x: getMagnAtWx(N, x, Watt) - (Ga), [Wa1, Wa2], 0.01
        )
    elif filter_type == "bandpass":
        Wn100 = []

        W1 = Wn[0]
        W2 = Wn[1]

        if abs(Wpass[1] - Watt[1]) > abs(Wpass[0] - Watt[0]):

            var = scipy.optimize.minimize_scalar(
                lambda delta: (
                    (getMagnAtWx(N, [Wn[0] - delta, Wn[1] + delta], Watt[0]) - Ga)
                ),
                bounds=(Wn[0] * 0.9, Wn[0] * 1.1),
            )
        else:

            var = scipy.optimize.minimize_scalar(
                lambda delta: (
                    (getMagnAtWx(N, [Wn[0] - delta, Wn[1] + delta], Watt[1]) - Ga)
                ),
                bounds=(Wn[1] * 0.9, Wn[1] * 1.1),
            )

        Wn100 = [Wn[0] - var.x, Wn[0] + var.x]

    elif filter_type == "bandstop":
        Wn100 = []
        WaA1 = Wn[0]
        WaA2 = Wn[0]
        while getMagnAtWx(N, [WaA1, Wn[1]], Watt[0]) < Ga:
            WaA1 *= 1.01
        while getMagnAtWx(N, [WaA2, Wn[1]], Watt[0]) > Ga:
            WaA2 *= 0.995

        magA1 = getMagnAtWx(N, [WaA1, Wn[1]], Watt[0])
        magA2 = getMagnAtWx(N, [WaA2, Wn[1]], Watt[0])

        Wn100.append(
            bisection(
                lambda x: getMagnAtWx(N, [x, Wn[1]], Watt[0]) - (Ga),
                [WaA1, WaA2],
                0.001,
            )
        )

        WaB1 = Wn[1]
        WaB2 = Wn[1]
        while getMagnAtWx(N, [Wn100[0], WaB1], Watt[1]) < Ga:
            WaB1 *= 0.99
        while getMagnAtWx(N, [Wn100[0], WaB2], Watt[1]) > Ga:
            WaB1 *= 1.005

        magB1 = getMagnAtWx(N, [Wn100[0], WaB1], Watt[1])
        magB2 = getMagnAtWx(N, [Wn100[0], WaB2], Watt[1])

        var = bisection(
            lambda x: getMagnAtWx(N, [Wn100[0], x], Watt[1]) - (Ga), [WaB2, WaB1], 0.001
        )

        Wn100.append(var)

        while getMagnAtWx(N, [WaA1, Wn[1]], Watt[0]) < Ga:
            WaA1 *= 1.01
        while getMagnAtWx(N, [WaA2, Wn[1]], Watt[0]) > Ga:
            WaA2 *= 0.995

        magA1 = getMagnAtWx(N, [WaA1, Wn[1]], Watt[0])
        magA2 = getMagnAtWx(N, [WaA2, Wn[1]], Watt[0])

        Wn100[0] = bisection(
            lambda x: getMagnAtWx(N, [x, Wn[1]], Watt[0]) - (Ga),
            [WaA1, WaA2],
            0.001,
        )

    # b2, a2 = filterfunc(N, Wn100, filter_type, True)
    # w, h2 = signal.freqs(b2, a2, np.logspace(1, 3, 1000))

    if filter_type == "bandpass" or filter_type == "bandstop":
        Wden = [0, 0]
        for index, element in enumerate(Wden):
            Wden[index] = Wn0[index] ** (1 - denorm) * Wn100[index] ** (denorm)
    else:
        Wden = Wn0 ** (1 - denorm) * Wn100 ** (denorm)

    b3, a3 = filterfunc(N, Wn, filter_type, True)
    # val1 = getMagnAtWx(N, Wden, Wpass[0])

    # b4, a4 = filterfunc(N, Wn0, filter_type, True)
    # w, h4 = signal.freqs(b4, a4, np.logspace(1, 3, 1000))

    # plt.semilogx(w, 20 * np.log10(abs(h4)), label='den 0%')
    # plt.semilogx(w, 20 * np.log10(abs(h2)), label='den 100%')

    return b3, a3


def graph_standard(filter_name, filter_type, N, Wpass, Watt, Gp, Ga, denorm, b, a):

    if filter_type == "lowpass":
        logsp = np.logspace(np.log10(0.01 * Wpass), np.log10(100 * Watt), 10000)
        w, h3 = signal.freqs(b, a, logsp)
    if filter_type == "highpass":
        logsp = np.logspace(np.log10(0.01 * Watt), np.log10(100 * Wpass), 10000)
        w, h3 = signal.freqs(b, a, logsp)
    if filter_type == "bandpass":
        logsp = np.logspace(np.log10(0.01 * Watt[0]), np.log10(100 * Watt[1]), 10000)
        w, h3 = signal.freqs(b, a, logsp)
    if filter_type == "bandstop":
        logsp = np.logspace(np.log10(0.01 * Wpass[0]), np.log10(100 * Wpass[1]), 10000)
        w, h3 = signal.freqs(b, a, logsp)

    plt.semilogx(w, 20 * np.log10(abs(h3)), label="den 50%")

    plt.title(
        "Transfer graph. "
        + filter_name
        + " "
        + filter_type
        + " filter. Order "
        + str(N)
        + "."
    )
    plt.xlabel("Frequency [radians / second]")
    plt.ylabel("Amplitude [dB]")
    plt.grid(which="both", axis="both")
    if filter_type == "lowpass":
        plt.fill(
            [0.01 * Wpass, 0.01 * Wpass, 100 * Watt, 100 * Watt],
            [0, 100 * (-Gp), 100 * (-Gp), 0],
            "0.9",
            lw=0,
        )  # zona prohibida: amplificación
        plt.fill(
            [0.01 * Wpass, 0.01 * Wpass, Wpass, Wpass],
            [100 * Gp, Gp, Gp, 100 * Gp],
            "0.9",
            lw=0,
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Watt, 100 * Watt, 100 * Watt, Watt], [0, 0, Ga, Ga], "0.9", lw=0
        )  # zona prohibida: sobre Ga
        plt.axis([Wpass * 0.1, Watt * 10, Ga - 10, 10])

    if filter_type == "highpass":
        plt.fill(
            [0.01 * Wpass, 0.01 * Wpass, 100 * Watt, 100 * Watt],
            [0, 100 * (-Gp), 100 * (-Gp), 0],
            "0.9",
            lw=0,
        )  # zona prohibida: amplificación
        plt.fill(
            [0.01 * Watt, 0.01 * Watt, Watt, Watt], [0, Ga, Ga, 0], "0.9", lw=0
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Wpass, 100 * Wpass, 100 * Wpass, Wpass],
            [Gp, Gp, 100 * Gp, 100 * Gp],
            "0.9",
            lw=0,
        )  # zona prohibida: sobre Ga
        plt.axis([Watt * 0.1, Wpass * 10, Ga - 10, 10])

    if filter_type == "bandpass":
        plt.fill(
            [1, 1, 100 * Watt[1], 100 * Watt[1]],
            [0, 100 * (-Gp), 100 * (-Gp), 0],
            "0.9",
            lw=0,
        )  # zona prohibida: amplificación
        plt.fill(
            [1, 1, Watt[0], Watt[0]], [0, Ga, Ga, 0], "0.9", lw=0
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Watt[1], Watt[1], 100 * Watt[1], 100 * Watt[1]],
            [0, Ga, Ga, 0],
            "0.9",
            lw=0,
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Wpass[0], Wpass[0], Wpass[1], Wpass[1]],
            [Gp, 100 * Gp, 100 * Gp, Gp],
            "0.9",
            lw=0,
        )  # zona prohibida: sobre Ga
        plt.axis([Watt[0] * 0.1, Watt[1] * 10, Ga - 10, 10])

    if filter_type == "bandstop":
        plt.fill(
            [1, 1, 100 * Wpass[1], 100 * Wpass[1]],
            [0, 100 * (-Gp), 100 * (-Gp), 0],
            "0.9",
            lw=0,
        )  # zona prohibida: amplificación
        plt.fill(
            [1, 1, Wpass[0], Wpass[0]], [Gp, 100 * Gp, 100 * Gp, Gp], "0.9", lw=0
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Wpass[1], Wpass[1], 100 * Wpass[1], 100 * Wpass[1]],
            [Gp, 100 * Gp, 100 * Gp, Gp],
            "0.9",
            lw=0,
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Watt[0], Watt[0], Watt[1], Watt[1]], [0, Ga, Ga, 0], "0.9", lw=0
        )  # zona prohibida: bajo Gp
        plt.axis([Wpass[0] * 0.1, Wpass[1] * 10, Ga - 10, 10])

    # stop

    plt.legend()
    plt.show()


def graph_atte(filter_name, filter_type, N, Wpass, Watt, Gp, Ga, denorm, b, a):

    if filter_type == "lowpass":
        logsp = np.logspace(np.log10(0.01 * Wpass), np.log10(100 * Watt), 10000)
        w, h3 = signal.freqs(b, a, logsp)
    if filter_type == "highpass":
        logsp = np.logspace(np.log10(0.01 * Watt), np.log10(100 * Wpass), 10000)
        w, h3 = signal.freqs(b, a, logsp)
    if filter_type == "bandpass":
        logsp = np.logspace(np.log10(0.01 * Watt[0]), np.log10(100 * Watt[1]), 10000)
        w, h3 = signal.freqs(b, a, logsp)
    if filter_type == "bandstop":
        logsp = np.logspace(np.log10(0.01 * Wpass[0]), np.log10(100 * Wpass[1]), 10000)
        w, h3 = signal.freqs(b, a, logsp)

    Gp *= -1
    Ga *= -1

    plt.semilogx(w, -20 * np.log10(abs(h3)), label="den 50%")

    plt.title(
        "Attenuation graph. "
        + filter_name.capitalize()
        + " "
        + filter_type
        + " filter. Order "
        + str(N)
        + "."
    )
    plt.xlabel("Frequency [radians / second]")
    plt.ylabel("Amplitude [dB]")
    plt.grid(which="both", axis="both")
    if filter_type == "lowpass":
        plt.fill(
            [0.01 * Wpass, 0.01 * Wpass, 100 * Watt, 100 * Watt],
            [0, 100 * (-Gp), 100 * (-Gp), 0],
            "0.9",
            lw=0,
        )  # zona prohibida: amplificación
        plt.fill(
            [0.01 * Wpass, 0.01 * Wpass, Wpass, Wpass],
            [100 * Gp, Gp, Gp, 100 * Gp],
            "0.9",
            lw=0,
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Watt, 100 * Watt, 100 * Watt, Watt], [0, 0, Ga, Ga], "0.9", lw=0
        )  # zona prohibida: sobre Ga
        plt.axis([Wpass * 0.1, Watt * 10, -10, Ga + 10])

    if filter_type == "highpass":
        plt.fill(
            [0.01 * Wpass, 0.01 * Wpass, 100 * Watt, 100 * Watt],
            [0, 100 * (-Gp), 100 * (-Gp), 0],
            "0.9",
            lw=0,
        )  # zona prohibida: amplificación
        plt.fill(
            [0.01 * Watt, 0.01 * Watt, Watt, Watt], [0, Ga, Ga, 0], "0.9", lw=0
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Wpass, 100 * Wpass, 100 * Wpass, Wpass],
            [Gp, Gp, 100 * Gp, 100 * Gp],
            "0.9",
            lw=0,
        )  # zona prohibida: sobre Ga
        plt.axis([Watt * 0.1, Wpass * 10, -10, Ga + 10])

    if filter_type == "bandpass":
        plt.fill(
            [1, 1, 100 * Watt[1], 100 * Watt[1]],
            [0, 100 * (-Gp), 100 * (-Gp), 0],
            "0.9",
            lw=0,
        )  # zona prohibida: amplificación
        plt.fill(
            [1, 1, Watt[0], Watt[0]], [0, Ga, Ga, 0], "0.9", lw=0
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Watt[1], Watt[1], 100 * Watt[1], 100 * Watt[1]],
            [0, Ga, Ga, 0],
            "0.9",
            lw=0,
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Wpass[0], Wpass[0], Wpass[1], Wpass[1]],
            [Gp, 100 * Gp, 100 * Gp, Gp],
            "0.9",
            lw=0,
        )  # zona prohibida: sobre Ga
        plt.axis([Watt[0] * 0.1, Watt[1] * 10, -10, Ga + 10])

    if filter_type == "bandstop":
        plt.fill(
            [1, 1, 100 * Wpass[1], 100 * Wpass[1]],
            [0, 100 * (-Gp), 100 * (-Gp), 0],
            "0.9",
            lw=0,
        )  # zona prohibida: amplificación
        plt.fill(
            [1, 1, Wpass[0], Wpass[0]], [Gp, 100 * Gp, 100 * Gp, Gp], "0.9", lw=0
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Wpass[1], Wpass[1], 100 * Wpass[1], 100 * Wpass[1]],
            [Gp, 100 * Gp, 100 * Gp, Gp],
            "0.9",
            lw=0,
        )  # zona prohibida: bajo Gp
        plt.fill(
            [Watt[0], Watt[0], Watt[1], Watt[1]], [0, Ga, Ga, 0], "0.9", lw=0
        )  # zona prohibida: bajo Gp
        plt.axis([Wpass[0] * 0.1, Wpass[1] * 10, -10, Ga + 10])

    # stop

    plt.legend()
    plt.show()
    return


def graph_filter(
    modality, filter_name, filter_type, N, Wpass, Watt, Gp, Ga, denorm, b, a
):
    if modality == "standard":
        graph_standard(filter_name, filter_type, N, Wpass, Watt, Gp, Ga, denorm, b, a)
    elif modality == "attenuation":
        graph_atte(filter_name, filter_type, N, Wpass, Watt, Gp, Ga, denorm, b, a)
    return


def return_p_z(b, a):
    z, p, g = signal.tf2zpk(b, a)

    zplane(z, p)
    poles = [[abs(x), -x.real / abs(x)] for x in p]

    print("poles:")
    print(len(poles))
    for element in poles:
        xi = element[1]
        if xi < 1e-15:
            xi = 0
        Q = 1 / (2 * xi)
        if Q > 1e14:
            Q = np.inf
        print("wo=" + str(element[0]) + "; xi= " + str(xi) + "; Q= " + str(Q))

    zeros = [[abs(x), -x.real / abs(x)] for x in z]

    print("zeros:")
    print(len(zeros))
    for element in zeros:
        xi = element[1]
        if xi < 1e-15:
            xi = 0

        if xi != 0:
            Q = 1 / (2 * xi)
            if Q > 1e14:
                Q = np.inf
        else:
            Q = np.inf

        print("wo=" + str(element[0]) + "; xi= " + str(xi) + "; Q= " + str(Q))

    return poles, zeros


def zplane(z, p, filename=None):
    """Plot the complex z-plane given a transfer function."""

    # get a figure/plot
    ax = plt.subplot(111)

    # create the unit circle
    # uc = patches.Circle((0,0), radius=1, fill=False,
    # color='black', ls='dashed')
    # ax.add_patch(uc)

    # Plot the zeros and set marker properties
    t1 = plt.plot(z.real, z.imag, "bo", ms=5)
    plt.setp(
        t1,
        markersize=10.0,
        markeredgewidth=1.0,
        markeredgecolor="k",
        markerfacecolor="b",
    )

    # Plot the poles and set marker properties
    t2 = plt.plot(p.real, p.imag, "rx", ms=5)
    plt.setp(
        t2,
        markersize=12.0,
        markeredgewidth=3.0,
        markeredgecolor="r",
        markerfacecolor="r",
    )

    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_position("center")
    ax.spines["right"].set_position("zero")
    ax.spines["top"].set_visible(False)
    ax.set_aspect("equal")
    ax.set_adjustable("datalim")
    ax.set_xlabel("Real", loc="right")
    ax.set_ylabel("Imaginary")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=-45, ha="left")
    ax.grid(which="both", axis="both")

    # set the ticks
    # r = 1.5; plt.axis('scaled'); plt.axis([-r, r, -r, r])
    # ticks = [-1, -.5, .5, 1]; plt.xticks(ticks); plt.yticks(ticks)

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)
