"""Шаблонный Солвер"""
from argparse import ArgumentParser
import numpy as np
import importlib


class GlobalSolver:
    """Класс реализует базовые методы на основе которых
    кастомизируются наследники в зависимости от задачи и типа сетки"""

    def __init__(self, input_filepath: str, output_filepath: str) -> None:
        """Инициализацаия
        Args:
            input_filepath: str - путь до файла с входными данными
            output_filepath: str - путь до файла в который сохраняется
            результат
        """
        self.input_filepath = input_filepath
        self.output_filepath = output_filepath
        self._parse_filedata()
        self._init_scheme_values()
        self._init_value()

    def _parse_filedata(self) -> None:
        """Чтение данных из файла и запись в аргументы экземпляра класса

        Args:
            None
        Return:
            None
        """
        with open(self.input_filepath, "r", encoding="utf-8") as f:
            data = f.read()
        rows = data.split("\n")
        for row in rows:
            val_key = row.split("!")
            val = float(val_key[0].strip())
            key = val_key[1].strip()
            setattr(self, key, val)
            print(f"	Load {key} value from file: {key} = {val}")
        return None

    def _init_scheme_values(self) -> None:
        """Инициализация основных параметров расчётной схемы

        Args:
            None
        Return:
            None
        """
        self.h = self.L / (self.NX - 1)
        self.dt = self.VNM * (self.h**2.0) / self.a
        self.NT = self.Time / self.dt
        print(
            f"	Computed h = {self.h} \n\
                Computed dt = {self.dt} \n\
                Computed NT = {self.NT}"
        )
        self.x = np.arange(0, self.L + self.h, self.h)
        assert len(self.x) == int(self.NX)
        return None

    def _init_value(self) -> None:
        """Инициализация массивов для хранения решения

        Args:
            None
        Return:
            None
        """
        self.v = np.zeros(int(self.NX))
        self.vn = np.zeros(int(self.NX))
        return None

    def init_value(self) -> None:
        """Обязательный для реализации в дочернем классе метод -
        осуществляет присвоение начальных условий для векторов решений"""
        pass

    def init_boundary(self) -> None:
        """Обязательный для реализации в дочернем классе метод -
        осуществляет инициализцаию граничных условий для вектора решения"""
        pass

    def run_scheme(self) -> None:
        """Обязательный для реализации в дочернем классе метод -
        описывает алгоритм вычислений на заданой сетке"""
        pass

    def save_to_file(self) -> None:
        """Обязательный для реализации в дочернем классе метод -
        сохраняет результаты решения задачи в файл"""
        pass

    def solve(self) -> None:
        """Последовательный вызов основных этапов решения задачи

        Args:
            None
        Return:
            None
        """
        self.init_value()
        self.init_boundary()
        self.run_scheme()
        self.save_to_file()
        print(f"	Work is over! \n	Results in {self.output_filepath}")
        return None


def main(input_filepath: str, output_filepath: str, solver_file: str) -> None:
    """Импортирование конкретной реализации солвера и запуск вычислений

    Args:
        input_filepath: str - путь до файла с входными данными
        output_filepath: str - путь до файла в который сохраняется
        результат
        solver_file: str - путь до модуля с конкретной реализацией
        солвера
    Return:
        None"""
    Solver = importlib.import_module(solver_file).Solver
    solver = Solver(input_filepath, output_filepath)
    solver.solve()
    return None


if __name__ == "__main__":
    """Запуск решения из CLI - считывание парамертов, запуск
    основной функции

    Args:
        input_filepath: str - путь до файла с входными данными
        output_filepath: str - путь до файла в который сохраняется
        результат
        solver_file: str - путь до модуля с конкретной реализацией
        солвера
    Return:
        None
    """
    parser = ArgumentParser(prog="Numerical Solver")
    parser.add_argument("input_filepath", type=str)
    parser.add_argument("output_filepath", type=str)
    parser.add_argument("solver_file", type=str)
    args = parser.parse_args()
    main(args.input_filepath, args.output_filepath, args.solver_file)
