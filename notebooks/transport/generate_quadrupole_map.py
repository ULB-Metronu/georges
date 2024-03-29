# Create the map for a combined dipole

import sympy as sy
from generate_transport_maps import generate_transport_matrix_tensor, write_matrix, write_tensor

if __name__ == "__main__":
    exact_dpp = True
    if exact_dpp:
        file = "transport_quadrupole_ex.py"
    else:
        file = "transport_quadrupole.py"

    with open(file, "w") as f:
        # HEADERS
        f.write(
            """import numpy as np
from numpy import sqrt, cos, sin, cosh, sinh
from numba import njit
from numba.typed import List as nList
    \n\n""",
        )

        # MATRIX
        f.write("@njit(cache=True)\n")
        if exact_dpp:
            f.write("""def compute_transport_quadrupole_ex_matrix(element_parameters: nList) -> np.ndarray:\n""")
        else:
            f.write("""def compute_transport_quadrupole_matrix(element_parameters: nList) -> np.ndarray:\n""")
        f.write(
            """
    L: float = element_parameters[0]
    k1: float = element_parameters[1]""",
        )
        if exact_dpp:
            f.write(
                """
    d: float = element_parameters[len(element_parameters) - 1]
    k1 = k1/(1+d)""",
            )
        f.write(
            """
    R = np.zeros((6, 6))
    R[4, 4] = 1
    R[5, 5] = 1 \n""",
        )

        f.write("\tif k1 == 0: \n")
        write_matrix(generate_transport_matrix_tensor(0, 0, 0, 1, 0), f, ntab=2)

        f.write("\tif k1 > 0: \n")
        k1 = sy.symbols("k1", real=True, positive=True, nonzero=True)
        write_matrix(generate_transport_matrix_tensor(0, k1, 0, 1, 1), f, ntab=2)

        f.write("\tif k1 < 0: \n")
        k1 = sy.symbols("k1", real=True, negative=True, nonzero=True)
        write_matrix(generate_transport_matrix_tensor(0, k1, 0, 1, -1), f, ntab=2)

        f.write("\treturn R\n\n\n")

        # TENSOR
        f.write("@njit(cache=True)\n")
        if exact_dpp:
            f.write("""def compute_transport_quadrupole_ex_tensor(element_parameters: nList) -> np.ndarray:\n""")
        else:
            f.write("""def compute_transport_quadrupole_tensor(element_parameters: nList) -> np.ndarray:\n""")
        f.write(
            """
    L: float = element_parameters[0]
    k1: float = element_parameters[1]""",
        )
        if exact_dpp:
            f.write(
                """
    d: float = element_parameters[3]
    k1 = k1/(1+d)""",
            )
        f.write(
            """
    T = np.zeros((6, 6, 6)) \n""",
        )

        f.write("\tif k1 > 0: \n")
        k1 = sy.symbols("k1", real=True, positive=True, nonzero=True)
        write_tensor(generate_transport_matrix_tensor(0, k1, 0, 2, 1), f, ntab=2)

        f.write("\tif k1 < 0: \n")
        k1 = sy.symbols("k1", real=True, negative=True, nonzero=True)
        write_tensor(generate_transport_matrix_tensor(0, k1, 0, 2, -1), f, ntab=2)

        f.write("\treturn T\n")
