""" Решение задачи течения в канале с движущейся крышкой
Схема № 2
"""
import copy
import numpy as np
from main import GlobalSolver


class Solver(GlobalSolver):
    def _init_scheme_values(self) -> None:
        self.h = self.H / (self.NY - 1)
        self.dt = self.VNM * (self.h ** 2.) / self.nu
        self.NT = self.Time / self.dt
        self.y = np.arange(0, self.H + self.h, self.h)
        assert len(self.y) == int(self.NY)
        return None

    def _init_value(self) -> None:
        self.v = np.zeros(int(self.NY))
        self.vn = np.zeros(int(self.NY))

    def init_value(self) -> None:
        self.v = np.zeros(int(self.NY))
        self.vn = np.zeros(int(self.NY))

    def init_boundary(self) -> None:
        self.vn[0] = self.U0
        self.vn[-1] = self.U1

    def run_scheme(self) -> None:
        self.init_boundary()
        self.v = copy.deepcopy(self.vn)
        for j in range(1, int(self.NT) - 1):
            for i in range(1, int(self.NY) - 1):
                self.vn[i] = self.v[i] + self.VNM * (self.v[i+1] - 2.0 * self.v[i] + self.v[i-1]) + self.A * self.dt
            self.init_boundary()
            self.v = copy.deepcopy(self.vn)

    def save_to_file(self) -> None:
        with open(self.output_filepath, "w") as f:
            f.write('variables = "y", "u"')
            for i in range(int(self.NY) - 1):
                f.write(f"{self.y[i]}, {self.v[i]} \n")
