"""Решение задачи конвекции №2 с помощью
схемы Леонарда (№ 4)
"""
import copy
import numpy as np
import matplotlib.pyplot as plt
from main import GlobalSolver


class Solver(GlobalSolver):
    """Реализация солвера для решения модельной задачи конвекции № 2
    с применением схемы Леонарда (№ 4)"""

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

        # Задаём мнимые точки в силу порядка аппроксимации
        self.x = np.arange(-2 * self.h, (self.NX + 2) * self.h, self.h)
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
        fig, ax = plt.subplots()
        ax.plot(self.x[2:-2], self.v[2:-2], label="numerical")
        ax.plot(
            self.x[2:-2],
            [init_func(x) for x in self.x[2:-2]],
            linestyle=':',
            linewidth=1,
            label="analytical",
        )
        ax.set_xlabel('x')
        ax.set_ylabel('u')

        fig.set_figheight(8)
        fig.set_figwidth(10)
        ax.legend()
        plt.savefig("s2_t2_burg.png", format='png')
        with open(self.output_filepath, "w") as f:
            f.write('Variables = "x", "u",\n')
            for i, x in enumerate(self.x[2:-2]):
                f.write(f"{x}, {self.v[i + 2]}\n")
        return None


def init_func(x: float) -> float:
    if (x < 0.2) or (0.4 <= x < 0.6) or (0.8 <= x):
        return 0.6
    elif 0.6 <= x < 0.8:
        return 0.4
    else:
        return 0.2
