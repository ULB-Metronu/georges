# Create the map for a combined dipole

import sympy as sy
from generate_transport_maps import generate_transport_matrix_tensor, write_matrix, write_tensor

if __name__ == "__main__":
    exact_dpp = True
    if exact_dpp:
        file = "transport_combined_dipole_ex.py"
    else:
        file = "transport_combined_dipole.py"

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
            f.write("""def compute_transport_combined_dipole_ex_matrix(element_parameters: nList) -> np.ndarray:\n""")
        else:
            f.write("""def compute_transport_combined_dipole_matrix(element_parameters: nList) -> np.ndarray:\n""")
        f.write(
            """
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    h = alpha / L
    k1: float = element_parameters[2]""",
        )
        if exact_dpp:
            f.write(
                """
    d: float = element_parameters[len(element_parameters) - 1]
    h = h/(1+d)
    k1 = k1/(1+d)
        """,
            )
        f.write(
            """
    R = np.zeros((6, 6))
    R[4, 4] = 1
    R[5, 5] = 1 \n""",
        )

        f.write("\tif h == 0:\n")

        f.write("\t\tif k1 == 0: \n")
        write_matrix(generate_transport_matrix_tensor(0, 0, 0, 1, 0), f)

        f.write("\t\tif k1 > 0: \n")
        k1 = sy.symbols("k1", real=True, positive=True, nonzero=True)
        write_matrix(generate_transport_matrix_tensor(0, k1, 0, 1, 1), f)

        f.write("\t\tif k1 < 0: \n")
        k1 = sy.symbols("k1", real=True, negative=True, nonzero=True)
        write_matrix(generate_transport_matrix_tensor(0, k1, 0, 1, -1), f)

        f.write("\tif h != 0:\n")
        h = sy.symbols("h", real=True, positive=True, nonzero=True)

        f.write("\t\tif k1 == 0: \n")
        write_matrix(generate_transport_matrix_tensor(h, 0, 0, 1, 1), f)

        f.write("\t\tif k1 != 0: \n")
        f.write("\t\t\tif h ** 2 + k1 == 0: \n")
        k1 = sy.symbols("k1", real=True, negative=True)
        write_matrix(generate_transport_matrix_tensor(h, k1, 0, 1, 0), f, ntab=4)

        f.write("\t\t\tif h ** 2 + k1 > 0: \n")
        f.write("\t\t\t\tif k1 > 0: \n")
        k1 = sy.symbols("k1", real=True, positive=True, nonzero=True)
        write_matrix(generate_transport_matrix_tensor(h, k1, 0, 1, 1), f, ntab=5)

        f.write("\t\t\t\tif k1 < 0: \n")
        k1 = sy.symbols("k1", real=True, negative=True, nonzero=True)
        write_matrix(generate_transport_matrix_tensor(h, k1, 0, 1, 1), f, ntab=5)

        f.write("\t\t\tif h ** 2 + k1 < 0: \n")
        k1 = sy.symbols("k1", real=True, negative=True)
        write_matrix(generate_transport_matrix_tensor(h, k1, 0, 1, -1), f, ntab=4)

        f.write("\treturn R\n\n\n")

        # TENSOR
        f.write("@njit(cache=True)\n")
        if exact_dpp:
            f.write("""def compute_transport_combined_dipole_ex_tensor(element_parameters: nList) -> np.ndarray:\n""")
        else:
            f.write("""def compute_transport_combined_dipole_tensor(element_parameters: nList) -> np.ndarray:\n""")
        f.write(
            """
    L: float = element_parameters[0]
    alpha: float = element_parameters[1]
    h = alpha / L
    k1: float = element_parameters[2]""",
        )
        if exact_dpp:
            f.write(
                """
    d: float = element_parameters[len(element_parameters) - 1]
    h = h/(1+d)
    k1 = k1/(1+d)
        """,
            )
        f.write(
            """
    T = np.zeros((6, 6, 6)) \n""",
        )

        f.write("\tif h == 0:\n")

        f.write("\t\tif k1 > 0: \n")
        k1 = sy.symbols("k1", real=True, positive=True, nonzero=True)
        write_tensor(generate_transport_matrix_tensor(0, k1, 0, 2, 0), f)

        f.write("\t\tif k1 < 0: \n")
        k1 = sy.symbols("k1", real=True, negative=True, nonzero=True)
        write_tensor(generate_transport_matrix_tensor(0, k1, 0, 2, 0), f)

        f.write("\tif h != 0:\n")
        f.write("\t\tif k1 == 0: \n")
        write_tensor(generate_transport_matrix_tensor(h, 0, 0, 2, 1), f)

        f.write("\t\tif k1 != 0: \n")
        f.write("\t\t\tif h ** 2 + k1 == 0: \n")
        k1 = sy.symbols("k1", real=True, negative=True, nonzero=True)
        write_tensor(generate_transport_matrix_tensor(h, k1, 0, 2, 0), f, ntab=4)

        f.write("\t\t\tif h ** 2 + k1 > 0: \n")
        f.write("\t\t\t\tif k1 > 0: \n")
        k1 = sy.symbols("k1", real=True, positive=True, nonzero=True)
        write_tensor(generate_transport_matrix_tensor(h, k1, 0, 2, 1), f, ntab=5)

        f.write("\t\t\t\tif k1 < 0: \n")
        k1 = sy.symbols("k1", real=True, negative=True, nonzero=True)
        write_tensor(generate_transport_matrix_tensor(h, k1, 0, 2, 1), f, ntab=5)

        f.write("\t\t\tif h ** 2 + k1 < 0: \n")
        k1 = sy.symbols("k1", real=True, negative=True, nonzero=True)
        write_tensor(generate_transport_matrix_tensor(h, k1, 0, 2, -1), f, ntab=4)

        f.write("\treturn T\n")
