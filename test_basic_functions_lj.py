"""
This script can be used to test your basic functions
that you will submit for the first submission in semester 2.

It won't test the actual values that your code gives, but it will
test that they take the right inputs and outputs

This version is for the Lennard-Jones project.
"""
from particle3D import Particle3D
import numpy as np
import sys


#############
# These three helper functions are to make this
# test program more convenient. You wouldn't
# need anything like them in normal work
def get_arg_names(function):
    """
    Get a list of the names of the arguments a function takes

    Parameters
    ----------
    function: callable
        Any python function

    Returns
    -------
    names: list
        A list of string names
    """
    import inspect

    return inspect.getfullargspec(function).args


def import_by_file_name(filename, module_name):
    """
    Import a python file as a module using the name of the file

    Parameters
    ----------
    filename: str
        The name or path to the file

    module_name: str
        A name to give the module after it is imported

    Returns
    -------
    module: module
        The module, just as if you had used "import""
    """
    import importlib

    spec = importlib.util.spec_from_file_location(module_name, filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def check_float_array_shape(result, shape, shape_text, func_name, filename):
    """
    Perform a series of checks on the type and shape of an array.

    Exceptions are raised for unexpected results.

    Parameters
    ----------
    result: objet
        Object tested to see if it is an array of the right shape and type

    shape: tuple
        Expected shape of the array

    shape_text: str
        Text to use in the error message if the shape is wrong

    func_name: str
        Name of the function to use in the error message

    filename: str
        Name of the file this came from to use in the error message

    Returns
    -------
    None

    """
    if not isinstance(result, np.ndarray):
        raise ValueError(
            f"{func_name} in {filename} was supposed to return an"
            f" array but it instead returned: {type(result)}"
        )

    ndim = len(shape)
    if result.ndim != ndim:
        raise ValueError(
            f"{func_name} in {filename} was supposed to return an {ndim}D array"
            f" but instead it was {result.ndim}D"
        )

    if result.shape != shape:
        raise ValueError(
            f"{func_name} in {filename} was supposed to return an array with shape"
            f" {shape_text} but actually returned {result.shape} for n=4"
        )

    if result.dtype.kind != "f":
        raise ValueError(
            f"{func_name} in {filename} was supposed to return a float"
            f" array but it was of a different type {result.dtype}"
        )

# end of helper functions
#############


def test_compute_separations(module, filename):
    """
    Run tests on the compute_separations function.
    Exceptions are raised if something is wrong.

    Parameters
    ----------
    module: module
        An imported modue

    filename: str
        The name or path to the file

    Returns
    -------
    None

    """
    try:
        module.compute_separations  # import compute_separations
    except:
        raise ValueError(f"compute_separations not found in {filename}")

    x1 = np.array([0.0, 0.0, 1.0])
    x2 = np.array([0.0, 0.0, 2.0])
    x3 = np.array([0.0, 0.0, 3.0])
    x4 = np.array([0.0, 0.0, 4.0])

    v = np.zeros(3)

    try:
        p1 = Particle3D("p1", 1.0, x1, v)
        p2 = Particle3D("p2", 1.0, x2, v)
        p3 = Particle3D("p2", 1.0, x3, v)
        p4 = Particle3D("p2", 1.0, x4, v)
    except:
        raise ValueError("Could not use your Particle3D class to create particles")

    particles = [p1, p2, p3, p4]
    box_size = 10.0

    args = get_arg_names(module.compute_separations)
    if len(args) != 2:
        raise ValueError(
            f"compute_separations function {filename} was supposed to take two"
            " argument, a list of particles and a box size, but yours has these:"
            f" {args}"
        )

    try:
        result = module.compute_separations(particles, box_size)
    except:
        raise ValueError(
            f"compute_separations function in {filename} crashed. Full error message"
            " in the block above this one"
        )
        raise

    if result is None:
        raise ValueError(
            f"compute_separations function in {filename} did not return anything"
        )

    if isinstance(result, tuple):
        raise ValueError(
            f"compute_separations function in {filename} returned more than one thing."
            " You might have a stray comma at the end?"
        )

    check_float_array_shape(
        result, (4, 4, 3), "(n, n, 3)", "compute_separations", filename
    )


def test_compute_forces_potential(module, filename):
    """
    Run tests on the test_compute_forces_potential function.
    Exceptions are raised if something is wrong.

    Parameters
    ----------
    module: module
        An imported modue

    filename: str
        The name or path to the file

    Returns
    -------
    None
    """

    try:
        module.compute_forces_potential
    except:
        raise ValueError(f"compute_forces_potential not found in {filename}")

    x1 = np.array([0.0, 0.0, 1.0])
    x2 = np.array([0.0, 0.0, 2.0])
    x3 = np.array([0.0, 0.0, 3.0])
    x4 = np.array([0.0, 0.0, 4.0])

    v = np.zeros(3)

    try:
        p1 = Particle3D("p1", 1.0, x1, v)
        p2 = Particle3D("p2", 1.0, x2, v)
        p3 = Particle3D("p2", 1.0, x3, v)
        p4 = Particle3D("p2", 1.0, x4, v)
    except:
        raise ValueError("Could not use your Particle3D class to create particles")

    particles = [p1, p2, p3, p4]

    # This is not an actually reasonable set of separations.
    # It is just something of the correct shape
    # (it is ones everywhere except the diagonal)
    seps = np.zeros((4, 4, 3))
    seps[:, :, 0] = (1 - np.eye(4)) * 2

    args = get_arg_names(module.compute_forces_potential)
    if len(args) != 1:
        raise ValueError(
            f"compute_forces_potential function {filename} was supposed to take one"
            f" argument, a 3D array. Yours has: {args}"
        )


    try:
        result = module.compute_forces_potential(seps)
    except:
        raise ValueError(
            f"compute_forces_potential function in {filename} crashed. Full error"
            " message in the block above this one"
        )
        raise

    if result is None:
        raise ValueError(
            f"compute_forces_potential function in {filename} did not return anything"
        )

    if isinstance(result, list):
        raise ValueError(
            f"compute_forces_potential function in {filename} returned a list [...] but"
            " should just return a pair (forces, potential)"
        )

    if not isinstance(result, tuple):
        raise ValueError(
            f"compute_forces_potential function in {filename} was supposed to return a"
            f" pair (forces, potential) but instead returned: {type(result)}"
        )

    if len(result) != 2:
        raise ValueError(
            f"compute_forces_potential function in {filename} was supposed to return a"
            " pair (forces, potential) but instead something with"
            f" {len(result)} objects in"
        )

    forces, potential = result

    check_float_array_shape(
        forces, (4, 3), "(n, 3)", "compute_forces_potential forces", filename
    )

    if not isinstance(potential, float):
        raise ValueError(
            "compute_forces_potential return value potential was supposed to be a"
            f" float, but yours was not, it was {type(potential)}"
        )


def run_tests(filename):
    """
    Run both tests on the given filename

    Parameters
    ----------
    filename: str
        Filename to read python code from
    """
    module = import_by_file_name(filename, "basic_functions")
    test_compute_separations(module, filename)
    test_compute_forces_potential(module, filename)

    # If we get this far with no errors then everything worked.
    print(
        "Everything tested okay, though this does not guarantee your code is right,"
        " just that it has the right kind of outputs and inputs."
    )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Syntax: %run  test_basic_functions_ss  filename.py")
        sys.exit(1)

    run_tests(sys.argv[1])
