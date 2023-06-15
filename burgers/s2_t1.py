"""Решение задачи Бюргерса №1 с помощью
схемы Леонарда (№ 4)
"""
import math
import copy
import numpy as np
import matplotlib.pyplot as plt
from main import GlobalSolver


class Solver(GlobalSolver):
    """Реализация солвера для решения модельной Бюргерса № 1
    с применением схемы Леонарда (№4)"""

    def _init_scheme_values(self) -> None:
        """Инициализация параметров схемы

        Args:
            None
        Return:
            None
        """
        self.h = self.L / (self.NX - 1)
        self.dt = self.CFL * self.h / self.C

        print(
            f"\tComputed h = {self.h}\n\
            \tComputed dt = {self.dt}\n\
            \tComputed NT = {self.NT}"
        )

        self.x = np.arange(-2 * self.h, (self.NX + 2) * self.h, self.h)
        print(self.x, self.NX, len(self.x), self.x[int(self.NX - 1)])
        return None

    def run_scheme(self) -> None:
        """Алгоритм вычисления решения на схеме

        Args:
            None
        Return:
            None
        """
        for j in range(1, int(self.NT)):
            for i in range(2, len(self.x) - 2):
                if self.v[i] > 0:
                    self.vn[i] = self.v[i] - self.C * self.dt / self.h / 6 * (
                        2 * self.v[i + 1] ** 2. / 2.
                        + 3 * self.v[i] ** 2. / 2.
                        - 6 * self.v[i - 1] ** 2. / 2.
                        + self.v[i - 2] ** 2. / 2.
                    )
                else:
                    self.vn[i] = self.v[i] - self.C * self.dt / self.h / 6 * (
                        - self.v[i + 2] ** 2. / 2.
                        - 3 * self.v[i] ** 2. / 2.
                        + 6 * self.v[i + 1] ** 2. / 2.
                        - 2 * self.v[i - 1] ** 2. / 2.
                    )
            self.init_boundary()
            self.v = copy.deepcopy(self.vn)
        return None

    def init_value(self) -> None:
        """Инициализация начальных значений

        Args:
            None
        Return:
            None
        """
        # Н.У. имеет вид: C0 + C1 * sin(m * pi * x / L)
        self.v = self.C0 + self.C1 * np.sin(self.x * self.m * math.pi / self.L)
        self.vn = np.zeros((len(self.x)))
        return None

    def init_boundary(self) -> None:
        """Инициализация граничных условий

        Args:
            None
        Return:
            None
        """
        self.vn[0] = self.vn[-4]
        self.vn[1] = self.vn[-3]
        self.vn[-2] = self.vn[2]
        self.vn[-1] = self.vn[3]
        return None

    def save_to_file(self) -> None:
        """Запись решения в файл

        Args:
            None
        Return:
            None
        """

        # g = 0.986745
        # fe = -0.28359
        k = self.m * math.pi / self.L
        # beta = k * self.h
        # g = np.sqrt(
        #     (
        #         (1 - self.CFL / 3 * (1 - np.cos(beta)) ** 2.0) ** 2.0
        #         + (
        #             self.CFL**2.0
        #             / 9
        #             * np.sin(beta) ** 2.0
        #             * (2 - np.cos(beta)) ** 2.0
        #         )
        #     )
        # )
        # fe = -np.arctan(
        #     self.CFL
        #     / 3
        #     * np.sin(beta)
        #     * (2 - np.cos(beta))
        #     / (1 - self.CFL / 3 * (1 - np.cos(beta)) ** 2.0)
        # )

        fig, ax = plt.subplots()
        ax.plot(self.x[2:-2], self.v[2:-2], label="numerical")
        ax.plot(
            self.x[2:-2],
            [
                math.sin(k * x - k * self.dt * self.NT * self.C)
                for x in self.x[2:-2]
            ],
            linestyle=':',
            linewidth=1,
            label="analytical",
        )
        # ax.plot(
        #     self.x[2: -2],
        #     [
        #         (g**self.NT) * math.sin(k * x + fe * self.NT)
        #         for x in self.x[2: -2]
        #     ],
        #     label="numer_estimate",
        # )
        ax.set_xlabel('x')
        ax.set_ylabel('u')

        fig.set_figheight(8)
        fig.set_figwidth(10)
        ax.legend()
        plt.savefig("s2_t1_burg.png", format='png')

        with open(self.output_filepath, "w") as f:
            f.write('Variables = "x", "u", "u_exac", "u_num_exac"\n')
            for i, x in enumerate(self.x[2: -2]):
                v_e = math.sin(k * x - k * self.dt * self.NT * self.C)
                # v_num_e = (g**self.NT) * math.sin(k * x + fe * self.NT)
                f.write(f"{x}, {self.v[i + 2]}, {v_e}, \n")  # {v_num_e} \n")
        return None
