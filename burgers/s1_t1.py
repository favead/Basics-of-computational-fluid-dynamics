"""Решение задачи Бюргерса №1 с помощью
явной противопоточной схемой первого порядка (№1)
"""
import math
import copy
import numpy as np
import matplotlib.pyplot as plt
from main import GlobalSolver


class Solver(GlobalSolver):
    """Реализация солвера для решения модельной Бюргерса № 1
    с применением явной противопоточной схемой первого порядка (№1)"""

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

        self.x = np.arange(-self.h, (self.NX + 1) * self.h, self.h)
        return None

    def run_scheme(self) -> None:
        """Алгоритм вычисления решения на схеме

        Args:
            None
        Return:
            None
        """
        for j in range(1, int(self.NT)):
            for i in range(1, len(self.x) - 1):
                if self.v[i] > 0:
                    self.vn[i] = self.v[i] - self.C * self.dt / self.h / 2 * (
                        self.v[i] ** 2.0 / 2 - self.v[i - 1] ** 2 / 2
                    )
                else:
                    self.vn[i] = self.v[i] - self.C * self.dt / self.h / 2 * (
                        self.v[i + 1] ** 2.0 / 2 - self.v[i] ** 2 / 2
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
        self.vn[0] = self.vn[-3]
        self.vn[1] = self.vn[-2]
        self.vn[-1] = self.vn[2]
        return None

    def save_to_file(self) -> None:
        """Запись решения в файл

        Args:
            None
        Return:
            None
        """
        k = self.m * math.pi / self.L
        beta = k * self.h
        g = np.sqrt(
            (1 - self.CFL * (1 - np.cos(beta)) ** 2.0)
            + (self.CFL * np.sin(beta)) ** 2.0
        )
        fe = -np.arctan(
            self.CFL * np.sin(beta)
            / (1 - self.CFL * (1 - np.cos(beta)))
        )

        fig, ax = plt.subplots()

        ax.plot(self.x[1:-1], self.v[1:-1], label="numer")
        ax.plot(
            self.x[1:-1],
            [
                math.sin(k * x - k * self.dt * self.NT * self.C)
                for x in self.x[1:-1]
            ],
            linestyle=':',
            linewidth=1,
            label="analytical",
        )
        ax.plot(
            self.x[1:-1],
            [
                (g**self.NT) * math.sin(k * x + fe * self.NT)
                for x in self.x[1:-1]
            ],
            linestyle='--',
            linewidth=1,
            label="numer_estimate",
        )
        ax.set_xlabel('x')
        ax.set_ylabel('u')

        fig.set_figheight(8)
        fig.set_figwidth(10)

        ax.legend()
        plt.savefig("s1_t1_burg.png", format='png')
        with open(self.output_filepath, "w") as f:
            f.write('Variables = "x", "u", "u_exac", "u_num_exac"\n')
            for i, x in enumerate(self.x[1:-1]):
                v_e = math.sin(k * x - k * self.dt * self.NT * self.C)
                v_num_e = (g**self.NT) * math.sin(k * x + fe * self.NT)
                f.write(f"{x}, {self.v[i + 1]}, {v_e}, {v_num_e} \n")
        return None
