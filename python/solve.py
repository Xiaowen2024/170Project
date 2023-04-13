"""Solves an instance.

Modify this file to implement your own solvers.

For usage, run `python3 solve.py --help`.
"""

import argparse
import math
import random
from pathlib import Path
from typing import Callable, Dict, List

import numpy as np

from instance import Instance
from point import Point
from solution import Solution
from file_wrappers import StdinFileWrapper, StdoutFileWrapper


def solve_naive(instance: Instance) -> Solution:

    cities = instance.cities
    grid_length = instance.grid_side_length
    cr = instance.coverage_radius
    pr = instance.penalty_radius
    dic = initialize(instance)

    sol = []
    sol_pen = {}
    for i in range(1000):
        to_select = Point(0, 0)
        dic_sorted = {k: v for k, v in sorted(dic.items(), key=lambda item: item[1], reverse=True)}
        towers = []
        while list(dic_sorted.values())[0] > 0:
            count = 0
            m = 100000000
            r2 = random.randint(100, 500)
            for point in dic_sorted:
                r1 = random.uniform(0, 0.5)
                count += 1
                if dic_sorted[point] == 0:
                    break
                if count == r2:
                    break
                num_overlaps = 0
                for k in towers:
                    if Point.distance_obj(point, k) <= pr:
                        num_overlaps += 1
                penalty = 170 * math.exp(0.17 * num_overlaps)
                covered_cities = dic_sorted[point]
                if (penalty / (covered_cities + r1)) < m:
                    to_select = point
                    m = penalty / (covered_cities + r1)

            towers.append(to_select)
            remove_all(dic_sorted, to_select, cr, grid_length, cities)
            del dic_sorted[to_select]
            dic_sorted = {k: v for k, v in sorted(dic_sorted.items(), key=lambda item: item[1], reverse=True)}

        tmp = Solution(instance=instance,towers=towers)
        sol.append(tmp)
        sol_pen[i] = tmp.penalty()
    sol_pen = {k: v for k, v in sorted(sol_pen.items(), key=lambda item: item[1])}
    return sol[list(sol_pen.keys())[0]]
    # return Solution(
    #     instance=instance,
    #     towers=towers,
    # )


SOLVERS: Dict[str, Callable[[Instance], Solution]] = {
    "naive": solve_naive
}


def initialize(instance: Instance):
    grid_length = instance.grid_side_length
    cr = instance.coverage_radius
    dic = {}
    cities = instance.cities
    count = 0
    for i in range(grid_length):
        for j in range(grid_length):
            p = Point(i, j)
            for c in cities:
                if Point.distance_obj(p, c) <= cr:
                    count = count + 1
            dic[p] = count
            count = 0
    return dic


def remove_one(dic: {}, p: Point, cr: int, grid_length: int):
    x = p.x
    y = p.y
    for i in range(x - cr, x + cr + 1):
        if i < 0 or i >= grid_length:
            continue
        for j in range(y - cr, y + cr + 1):
            if j < 0 or j >= grid_length:
                continue
            p2 = Point(i, j)
            if Point.distance_obj(p, p2) <= cr and p2 in dic:
                if dic[p2] > 0:
                    dic[p2] = dic[p2] - 1


def remove_all(dic: {}, p: Point, cr: int, grid_length: int, cities: List[Point]):
    x = p.x
    y = p.y
    for i in range(x - cr, x + cr + 1):
        if i < 0 or i >= grid_length:
            continue
        for j in range(y - cr, y + cr + 1):
            if j < 0 or j >= grid_length:
                continue
            p2 = Point(i, j)
            if Point.distance_obj(p, p2) <= cr and p2 in dic and p2 in cities:
                remove_one(dic, p2, cr, grid_length)


# You shouldn't need to modify anything below this line.
def infile(args):
    if args.input == "-":
        return StdinFileWrapper()

    return Path(args.input).open("r")


def outfile(args):
    if args.output == "-":
        return StdoutFileWrapper()

    return Path(args.output).open("w")


def main(args):
    with infile(args) as f:
        instance = Instance.parse(f.readlines())
        solver = SOLVERS[args.solver]
        solution = solver(instance)
        assert solution.valid()
        with outfile(args) as g:
            print("# Penalty: ", solution.penalty(), file=g)
            solution.serialize(g)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve a problem instance.")
    parser.add_argument("input", type=str, help="The input instance file to "
                                                "read an instance from. Use - for stdin.")
    parser.add_argument("--solver", required=True, type=str,
                        help="The solver type.", choices=SOLVERS.keys())
    parser.add_argument("output", type=str,
                        help="The output file. Use - for stdout.",
                        default="-")
    main(parser.parse_args())
