"""Решение модельной задачи конвекции №17
с помощью явной противопоточной схемой первого порядка (№1)"""
import copy
import numpy as np
import matplotlib.pyplot as plt
from main import GlobalSolver


class Solver(GlobalSolver):
    """Реализация солвера для решения модельной задачи конвекции № 17
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

        self.x = np.arange(-self.h, self.L + 2 * self.h, self.h)
        return None

    def init_value(self) -> None:
        """Инициализация начальных значений решения

        Args:
            None
        Return:
            None
        """
        full_dim = len(self.x)
        self.v = np.array([init_func(x) for x in self.x])
        self.vn = np.zeros((full_dim))
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

    def run_scheme(self) -> None:
        """Запуск решения на схеме

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

    def save_to_file(self) -> None:
        """Запись решения в файл

        Args:
            None
        Return:
            None
        """

        fig, ax = plt.subplots()
        ax.plot(self.x[1:-1], self.v[1:-1], label="numerical")
        ax.plot(
            self.x[1:-1],
            [init_func(x) for x in (self.x[1:-1])],
            linestyle=':',
            linewidth=1,
            label="analytical",
        )
        ax.set_xlabel('x')
        ax.set_ylabel('u')

        fig.set_figheight(8)
        fig.set_figwidth(10)
        ax.legend()
        plt.savefig("s1_t2_burg.png", format='png')

        with open(self.output_filepath, "w") as f:
            f.write('Variables = "x", "u",\n')
            for i, x in enumerate(self.x[1:-1]):
                f.write(f"{x}, {self.v[i + 1]}\n")
        return None


def init_func(x: float) -> float:
    if (x < 0.2) or (0.4 <= x < 0.6) or (0.8 <= x):
        return 0.6
    elif 0.6 <= x < 0.8:
        return 0.4
    else:
        return 0.2
