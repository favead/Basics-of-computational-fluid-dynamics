from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np


def g(Cu, beta):
    return (
        np.sqrt(
            (1 + Cu**2. * (np.cos(beta) - 1))**2.
            + np.sin(beta)**2. * Cu ** 2.
        )
    )


def g2(Cu, beta):
    val = (
        np.sqrt(
            (1 - Cu / 3 * (1 - np.cos(beta))**2.)**2.
            + (Cu ** 2. / 9 * np.sin(beta)**2. * (2 - np.cos(beta))**2.)
        )
    )
    return val


def create_dissipation_surface():
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    beta_s = np.linspace(0., 3., num=100)
    Cu_s = np.linspace(0, 2., num=100)
    
    X, Y = np.meshgrid(Cu_s, beta_s)
    g_values = np.array([g2(Cu, beta) for (Cu, beta) in zip(X, Y)])
    surf = ax.plot_surface(X, Y, g_values,
                           cmap=cm.coolwarm,
                           linewidth=0.2,
                           antialiased=False,
                           alpha=0.9)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    fig.set_figheight(8)
    fig.set_figwidth(10)

    ax.set_xlabel('Число куранта, Cu')
    ax.set_ylabel(r'Номер гармоники, $\beta$')
    ax.set_zlabel(r'Модуль коэффициента переноса, $|g|$')

    ax.view_init(20, 35)
    # plt.show()
    plt.savefig("dissipation_2.png", format='png',
                )


if __name__ == "__main__":
    create_dissipation_surface()
