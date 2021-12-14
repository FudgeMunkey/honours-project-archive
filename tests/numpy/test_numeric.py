# These are existing tests from numpy/core/tests/test_numeric.py which have been rewritten

import functools
from hypothesis import given, assume, strategies as st
from hypothesis.extra import numpy as hynp
import numpy as np
import pytest
import warnings

from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
    assert_array_equal,
    assert_equal,
)


class TestResize:
    @given(
        hynp.arrays(dtype=np.float64, shape=hynp.array_shapes()),
        st.integers(min_value=1, max_value=1e6),
        st.integers(min_value=0),
    )
    def test_copies(self, A, copies, axis):
        """
        Rewrite of test_copies
        """
        assume(axis < A.ndim)

        scaler = np.ones(A.ndim, dtype=np.int64)
        scaler[axis] = copies
        new_shape = A.shape * scaler

        out = np.resize(A, new_shape=new_shape)
        assert out.size >= A.size

    @given(
        hynp.arrays(dtype=np.int64, shape=hynp.array_shapes()),
    )
    def test_zeroresize(self, A):
        """Resize the array to (0,)
        Rewrite of test_zeroresize
        """
        Ar = np.resize(A, (0,))
        assert_array_equal(Ar, np.array([]))
        assert_equal(A.dtype, Ar.dtype)

        Ar = np.resize(A, (0, 2))
        assert_equal(Ar.shape, (0, 2))

        Ar = np.resize(A, (2, 0))
        assert_equal(Ar.shape, (2, 0))

    @given(
        hynp.arrays(
            dtype=[("a", np.float64)],
            shape=hynp.array_shapes(),
            elements=st.integers(min_value=0, max_value=0),
        ),
        hynp.array_shapes(),
    )
    def test_reshape_from_zero(self, A, new_shape):
        """Resize an array of zeros
        Rewrite of test_reshape_from_zero
        """
        Ar = np.resize(A, new_shape=new_shape)
        assert_array_equal(Ar, np.zeros(new_shape, Ar.dtype))
        assert_equal(A.dtype, Ar.dtype)

    @given(
        hynp.arrays(dtype=np.int64, shape=hynp.array_shapes()),
        st.tuples(st.integers(max_value=-1), st.integers(max_value=-1)),
    )
    def test_negative_resize(self, A, new_shape):
        """
        Resize the array using negative shapes
        Rewrite of test_negative_resize
        """
        with pytest.raises(ValueError, match=r"negative"):
            np.resize(A, new_shape=new_shape)

    @given(
        hynp.arrays(dtype=np.float64, shape=hynp.array_shapes()), hynp.array_shapes()
    )
    def test_subclass(self, A, new_shape):
        """
        Rewrite of test_subclass
        """

        class MyArray(np.ndarray):
            __array_priority__ = 1.0

        my_arr = np.array(A).view(MyArray)
        assert type(np.resize(my_arr, new_shape=new_shape)) is MyArray


class TestNonarrayArgs:
    # # check that non-array arguments to functions wrap them in arrays

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=hynp.array_shapes(),
            elements=(st.floats(allow_nan=False)),
        ),
        st.floats(allow_nan=False),
        st.floats(allow_nan=False),
    )
    def test_clip(self, A, clip_min, clip_max):
        """
        Rewrite of test_clip
        """
        assume(clip_min < clip_max)

        out = np.clip(A, clip_min, clip_max)

        # Check min/max values in out array
        assert np.amin(out) >= clip_min
        assert np.amax(out) <= clip_max

        # Compare each value (Differential testing? - very slow)
        flat_A = np.ndarray.flatten(A)
        flat_A = [
            fa
            if clip_min <= fa <= clip_max
            else (clip_min if fa < clip_min else clip_max)
            for fa in flat_A
        ]
        flat_out = np.ndarray.flatten(out)

        assert all([a == o for a, o in zip(flat_A, flat_out)])

    @given(
        hynp.arrays(dtype=np.float64, shape=hynp.array_shapes()),
        st.integers(min_value=0),
    )
    def test_count_nonzero(self, A, axis):
        """
        Rewrite of test_count_nonzero
        """
        assume(axis < A.ndim)

        if A.ndim > 1:
            out = np.count_nonzero(A, axis=axis)
            assert len(out) <= A.size
            assert np.sum(out) <= A.size
        else:
            out = np.count_nonzero(A)
            assert out <= A.size

    @given(
        hynp.arrays(dtype=np.float64, shape=hynp.array_shapes()),
        st.integers(min_value=0),
    )
    def test_cumproduct(self, A, axis):
        """
        Rewrite of test_cumproduct
        """
        assume(axis < A.ndim)

        if A.ndim > 1:
            out = np.cumproduct(A, axis=axis)
            assert out.size == A.size
        else:
            out = np.cumproduct(A)
            assert out.size == A.size

    @given(
        hynp.arrays(dtype=np.float64, shape=hynp.array_shapes(min_dims=2, max_dims=2))
    )
    def test_diagonal(self, A):
        """
        Rewrite of test_diagonal
        """

        index = np.eye(A.shape[0], A.shape[1], dtype=np.bool8)
        assert_array_equal(np.diagonal(A), A[index])

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=hynp.array_shapes(),
            elements=st.floats(allow_nan=False, allow_infinity=False),
        ),
    )
    def test_mean(self, A):
        """
        Rewrite of test_mean
        """
        if A.size == 0:
            with warnings.catch_warnings(record=True) as w:
                warnings.filterwarnings("always", "", RuntimeWarning)
                assert np.isnan(np.mean(A))
                assert w[0].category is RuntimeWarning
        else:
            assert np.mean(A) == np.sum(A) / A.size

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=hynp.array_shapes(),
            elements=st.floats(allow_nan=False, allow_infinity=False),
        )
    )
    def test_ptp(self, A):
        """
        Rewrite of test_ptp
        """
        assert np.ptp(A) == np.max(A) - np.min(A)

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=hynp.array_shapes(max_dims=1),
            elements=st.floats(allow_nan=False, allow_infinity=False),
        )
    )
    def test_prod(self, A):
        """
        Rewrite of test_prod
        """
        out = np.prod(A)
        if np.isnan(out):
            return

        assert type(out) == np.float64
        assert out == functools.reduce(lambda a, b: a * b, np.ndarray.flatten(A))

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=hynp.array_shapes(),
        )
    )
    def test_ravel(self, A):
        """
        Rewrite of test_ravel
        """
        assert_array_equal(np.ravel(A), np.ndarray.flatten(A))  # order="C"

        order = ["C", "F", "A", "K"]
        for o in order:
            out = np.ravel(A, o)
            assert type(out) == np.ndarray
            assert out.size == A.size
            assert out.shape == (A.size,)

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=hynp.array_shapes(),
        ),
        st.integers(min_value=0, max_value=1e6),
    )
    def test_repeat(self, A, repeats):
        """
        Rewrite of test_repeat
        Note: Very easily runs out of memory for large `repeats`
        """
        # Size limit of 50M~ otherwise memory errors
        assume(A.size * repeats < 50e6)

        out = np.repeat(A, repeats)

        if repeats == 0:
            assert out.size == 0
        else:
            assert out.size == A.size * repeats

    @given(
        hynp.arrays(dtype=np.float64, shape=hynp.array_shapes()), hynp.array_shapes()
    )
    def test_reshape(self, A, new_shape):
        """
        Rewrite of test_reshape
        """
        assume(np.prod(A.shape) == np.prod(new_shape))

        out = np.reshape(A, new_shape)
        assert out.size == A.size

        # Round Trip
        out = np.reshape(out, A.shape)
        assert_array_equal(out, A)

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=hynp.array_shapes(),
            elements=st.floats(allow_nan=False, allow_infinity=False),
        ),
        st.integers(min_value=0, max_value=15),
    )
    def test_round_1(self, A, decimals):
        """
        Rewrite of test_round
        Note: Very easily runs out of memory for large `repeats`
        """
        out = np.around(A, decimals)

        assert len(out) == len(A)
        out_sum = np.sum(out)
        if np.isnan(out_sum) or np.isinf(out_sum):
            return  # Overflow

        assert abs(np.sum(out)) - abs(np.sum(A)) <= 5 / (10 ** decimals) * len(A)

    @given(
        st.floats(allow_nan=False, allow_infinity=False),
        st.integers(min_value=0, max_value=15),
    )
    def test_round_2(self, V, decimals):
        """
        Rewrite of test_round
        """
        out = np.float64(V)
        out = out.round(decimals)

        assert isinstance(out, np.float64)

    @given(
        hynp.arrays(dtype=np.float64, shape=hynp.array_shapes(max_dims=1)), st.floats()
    )
    def test_searchsorted(self, A, V):
        """
        Rewrite of test_searchsorted
        """
        out = np.searchsorted(A, V)
        assert out <= A.size

    @given(
        hynp.arrays(dtype=np.float64, shape=hynp.array_shapes()),
        st.integers(min_value=0),
    )
    def test_size(self, A, axis):
        """
        Rewrite of test_size
        """
        assume(axis < A.ndim)

        assert np.size(A) == A.size
        assert np.size(A, axis=axis) <= A.shape[axis]

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=hynp.array_shapes(),
            elements=st.floats(allow_nan=False, allow_infinity=False),
        ),
        st.integers(min_value=0),
    )
    def test_std(self, A, axis):
        """
        Rewrite of test_std
        Upper bound: 10.1111/j.1467-9639.2004.00157.x
        """
        assume(axis < A.ndim)

        out = np.std(A)

        if A.size == 0:
            with warnings.catch_warnings(record=True) as w:
                warnings.filterwarnings("always", "", RuntimeWarning)
                assert np.isnan(np.std(A))
                assert w[0].category is RuntimeWarning
        else:
            assert out >= 0
            assert out <= np.ptp(A)

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=hynp.array_shapes(),
            elements=st.floats(allow_nan=False, allow_infinity=False),
        ),
        st.integers(min_value=0),
        st.booleans(),
    )
    def test_sum(self, A, axis, keepdims):
        """
        Rewrite of test_sum
        """
        assume(axis < A.ndim)

        if A.ndim == 1:
            out = np.sum(A)
            assert isinstance(out, np.float64)
            assert out >= np.min(A) * A.size - 1e7
            assert out <= np.max(A) * A.size + 1e7
            assert_almost_equal(
                out, functools.reduce(lambda a, b: a + b, np.ndarray.flatten(A))
            )
        else:
            out = np.sum(A, axis=axis)

            out = np.sum(A, axis=axis, keepdims=keepdims)

    @given(
        hynp.arrays(np.float64, shape=hynp.array_shapes()),
        hynp.arrays(np.int64, shape=hynp.array_shapes(min_side=1)),
    )
    def test_take(self, A, indices):
        """
        Rewrite of test_take
        """
        a_min = np.min(indices)
        assume(abs(a_min) >= 0 and abs(a_min) < A.size)
        assume(np.max(indices) < A.size)

        out = np.take(A, indices=indices)
        assert out.shape == indices.shape

    @given(hynp.arrays(np.float64, shape=hynp.array_shapes()))
    def test_transpose(self, A):
        """
        Rewrite of test_transpose
        """
        out = np.transpose(A)
        assert out.size == A.size
        assert out.shape == A.shape[::-1]

        out = np.transpose(out)
        assert_array_equal(out, A)

    @given(
        hynp.arrays(
            np.float64,
            shape=hynp.array_shapes(),
            elements=st.floats(allow_nan=False, allow_infinity=False),
        ),
        st.integers(min_value=0),
    )
    def test_var(self, A, axis):
        """
        Rewrite of test_var
        """
        assume(axis < A.ndim)

        out = np.var(A)

        if A.size == 0:
            with warnings.catch_warnings(record=True) as w:
                warnings.filterwarnings("always", "", RuntimeWarning)
                assert np.isnan(out)
                assert w[0].category is RuntimeWarning
        else:
            assert out >= 0
            assert out <= np.ptp(A) ** 2

            out_var = np.var(A, axis=axis)
            out_std = np.std(A, axis=axis)

            assert_allclose(out_var, out_std ** 2)
