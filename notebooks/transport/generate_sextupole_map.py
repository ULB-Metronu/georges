# Create the map for a combined dipole

import sympy as sy
from generate_transport_maps import generate_transport_matrix_tensor, write_matrix, write_tensor

if __name__ == "__main__":
    exact_dpp = True
    if exact_dpp:
        file = "transport_sextupole_ex.py"
    else:
        file = "transport_sextupole.py"

    with open(file, "w") as f:
        # HEADERS
        f.write(
            """import numpy as np
from numba.typed import List as nList
from numba import njit
    \n\n""",
        )

        # MATRIX
        f.write("@njit(cache=True)\n")
        if exact_dpp:
            f.write("""def compute_transport_sextupole_ex_matrix(element_parameters: nList) -> np.ndarray:\n""")
        else:
            f.write("""def compute_transport_sextupole_matrix(element_parameters: nList) -> np.ndarray:\n""")
        f.write(
            """
    L: float = element_parameters[0]
    R = np.zeros((6, 6))
    R[4, 4] = 1
    R[5, 5] = 1 \n""",
        )

        write_matrix(generate_transport_matrix_tensor(0, 0, 0, 1, 0), f, ntab=1)

        f.write("\treturn R\n\n\n")

        # TENSOR
        f.write("@njit(cache=True)\n")
        if exact_dpp:
            f.write("""def compute_transport_sextupole_ex_tensor(element_parameters: nList) -> np.ndarray:\n""")
        else:
            f.write("""def compute_transport_sextupole_tensor(element_parameters: nList) -> np.ndarray:\n""")
        f.write(
            """
    L: float = element_parameters[0]
    k2: float = element_parameters[1]""",
        )
        if exact_dpp:
            f.write(
                """
    d: float = element_parameters[len(element_parameters) - 1]
    k2 = k2/(1+d)""",
            )
        f.write(
            """
    T = np.zeros((6, 6, 6)) \n""",
        )

        k2 = sy.symbols("k2", real=True)
        write_tensor(generate_transport_matrix_tensor(0, 0, k2, 2, 0), f, ntab=1)

        f.write("\treturn T\n")
