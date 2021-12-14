# These are novel tests written for NumPy

import math
from hypothesis import given, strategies as st
from hypothesis.control import assume
from hypothesis.extra import numpy as hynp
import numpy as np
from numpy.core.einsumfunc import einsum
from numpy.testing import assert_array_equal


class TestNovel:
    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=hynp.array_shapes(max_dims=1),
            elements=(st.floats(allow_nan=False)),
        )
    )
    def test_einsum(self, A):
        assert einsum("i,i", A, A, optimize=True) == einsum("i,i", A, A, optimize=False)

    @given(
        hynp.arrays(dtype=np.int64, shape=hynp.array_shapes()),
        hynp.array_shapes(),
    )
    def test_grow_shrink(self, A, new_shape):
        """
        Resize the array to be larger and then to the original size
        """
        assume(np.prod(new_shape) >= A.size)

        old_shape = A.shape
        new_array = np.resize(A, new_shape)
        old_array = np.resize(new_array, old_shape)

        assert_array_equal(A, old_array)

    @given(
        hynp.arrays(
            dtype=st.one_of(
                hynp.floating_dtypes(endianness="?"), hynp.complex_number_dtypes()
            ),
            shape=hynp.array_shapes(),
        )
    )
    def test_rint(self, arr):
        """Test the numpy rint function

        Generate numpy arrays with
            - Big and small endianness
            - 16, 32 and 64 bit sizes
        """

        flat_arr = np.ndarray.flatten(arr)
        flat_rounded = np.ndarray.flatten(np.rint(arr))

        # Number of elements doesn't change
        assert arr.size == len(flat_rounded)

        # Check properties of each element
        for i in range(len(flat_arr)):
            current_arr_value = flat_arr[i]
            current_rounded_value = flat_rounded[i]

            # Check infs and nans stay the same
            if current_arr_value == float("inf"):
                assert current_rounded_value == float("inf")
            elif current_arr_value == float("-inf"):
                assert current_rounded_value == float("-inf")
            elif isinstance(current_arr_value, np.complexfloating):
                assert 1 == 1
            elif math.isnan(current_arr_value):
                assert math.isnan(current_rounded_value)
            elif isinstance(current_arr_value, (np.integer, np.floating)):
                # Check ints do not change and floats become ints
                assert current_rounded_value % 1 == 0

                # Check the rounded values are within 0.5 of the old value
                assert abs(current_arr_value - current_rounded_value) <= 0.5
            else:
                assert current_arr_value == current_rounded_value
