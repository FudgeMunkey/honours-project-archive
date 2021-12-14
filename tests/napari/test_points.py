from copy import copy
from hypothesis import given, strategies as st
from hypothesis.control import assume
from hypothesis.extra import numpy as hynp
from hypothesis.strategies._internal.core import shared
import numpy as np
from napari.layers import Points
from numpy.testing import assert_almost_equal, assert_array_equal

# Strategies
random_shape_st = hynp.array_shapes(min_side=2, max_side=4, min_dims=2, max_dims=2)
random_point_element_st = st.one_of(
    st.floats(allow_infinity=False, min_value=-1e10, max_value=1e10),
    st.integers(min_value=-1e10, max_value=1e10),
)
random_points_st = hynp.arrays(
    dtype=np.float64, shape=random_shape_st, elements=random_point_element_st
)
random_points_2D_st = hynp.arrays(
    dtype=np.float64,
    shape=hynp.array_shapes(min_side=2, max_side=2, min_dims=2, max_dims=2),
    elements=random_point_element_st,
)
random_opacity_st = st.floats(
    allow_nan=False, allow_infinity=False, min_value=0.0, max_value=1.0
)
random_blending_st = st.sampled_from(["opaque", "translucent", "additive"])
random_symbol_st = st.sampled_from(
    [
        "arrow",
        "clobber",
        "cross",
        "diamond",
        "disc",
        "hbar",
        "ring",
        "square",
        "star",
        "tailed_arrow",
        "triangle_down",
        "triangle_up",
        "vbar",
        "x",
    ]
)


class TestPoints:
    # Tests for the Points layer

    @given(
        hynp.arrays(
            dtype=hynp.unicode_string_dtypes(), shape=hynp.array_shapes(max_dims=1)
        )
    )
    def test_empty_points_with_properties(self, labels):
        """Instantiate empty points layer with properties. Generate an array of unicode labels.
        Rewrite of test_empty_points_with_properties, test_empty_points_with_properties_list
        """

        properties = {
            "label": labels,
            "cont_prop": np.array([0], dtype=float),
        }
        pts = Points(property_choices=properties)
        current_props = {k: v[0] for k, v in properties.items()}
        np.testing.assert_equal(pts.current_properties, current_props)

        assert pts.properties["cont_prop"].dtype == float

        pts.add([10, 10])
        pts.add([10, 20])
        props = {
            "label": np.array([labels[0]] * 2),
            "cont_prop": np.array([0] * 2, dtype=float),
        }

        np.testing.assert_equal(pts.properties, props)

    @given(random_points_st)
    def test_random_points(self, data):
        """Instantiate Points class with random array of 2D points
        Rewrite of test_empty_points, test_random_points, test_empty_points_array, test_negative_points, test_integer_points, test_3D_points, test_4D_points
        """
        pts = Points(data)

        if pts.data.size == 0:
            assert pts.data.shape == (0, 2)
            return

        np.testing.assert_array_equal(pts.data, data)
        assert pts.data.ndim == data.ndim
        assert pts._view_data.ndim == data.ndim
        assert pts.data.size == data.size
        assert len(pts.selected_data) == 0

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=shared(random_shape_st, key="test_changing_points"),
            elements=random_point_element_st,
        ),
        hynp.arrays(
            dtype=np.float64,
            shape=shared(random_shape_st, key="test_changing_points"),
            elements=random_point_element_st,
        ),
    )
    def test_changing_points(self, A, B):
        """Change the points in the Points layer
        Rewrite of test_changing_points
        """

        data_A = np.array(A)
        data_B = np.array(B)

        pts = Points(data_A)
        pts.data = data_B

        assert_array_equal(pts.data, data_B)
        assert pts.ndim == data_B.shape[1]
        assert pts._view_data.ndim == 2
        assert len(pts.data) == len(B)

    @given(
        random_points_st,
        st.lists(elements=st.integers(min_value=0), min_size=1),
        st.lists(elements=st.integers(min_value=0), min_size=1),
    )
    def test_selecting_points(self, data, selected_data, other_selected_data):
        """Test selecting points.
        Rewrite of test_selecting_points
        """
        assume(
            np.max(selected_data) < data.shape[0]
            and np.max(other_selected_data) < data.shape[0]
        )

        selected_data = set(selected_data)
        other_selected_data = set(other_selected_data)

        layer = Points(data)
        layer.mode = "select"
        layer.selected_data = selected_data
        assert layer.selected_data == selected_data

        # test switching to 3D
        layer._slice_dims(ndisplay=3)
        assert layer.selected_data == selected_data

        # select different points while in 3D mode
        layer.selected_data = other_selected_data
        assert layer.selected_data == other_selected_data

        # selection should persist when going back to 2D mode
        layer._slice_dims(ndisplay=2)
        assert layer.selected_data == other_selected_data

        # selection should persist when switching between between select and pan_zoom
        layer.mode = "pan_zoom"
        assert layer.selected_data == other_selected_data
        layer.mode = "select"
        assert layer.selected_data == other_selected_data

        # add mode should clear the selection
        layer.mode = "add"
        assert layer.selected_data == set()

    @given(
        hynp.arrays(
            dtype=np.float64,
            shape=shared(random_shape_st, key="test_adding_points"),
            elements=random_point_element_st,
        ),
        hynp.arrays(
            dtype=np.float64,
            shape=shared(random_shape_st, key="test_adding_points"),
            elements=random_point_element_st,
        ),
    )
    def test_adding_points(self, A, B):
        """Test selecting points.
        Rewrite of test_adding_points, test_adding_points_to_empty
        """
        layer = Points(A)
        assert len(layer.data) == A.shape[0]

        layer.add(B)

        assert len(layer.data) == A.shape[0] + B.shape[0]
        assert_array_equal(layer.data, np.concatenate((A, B)))
        assert layer.selected_data == set(range(A.shape[0], A.shape[0] + B.shape[0]))
        # assert layer.selected_data == {A.shape[0] + B.shape[0] - 1}

        layer.remove_selected()
        assert len(layer.data) == A.shape[0]
        assert_array_equal(layer.data, A)

    @given(random_points_st, st.lists(elements=st.integers(min_value=0), min_size=1))
    def test_removing_selected_points(self, data, delete_indices):
        """Test selecting points.
        Rewrite of test_removing_selected_points
        """
        assume(np.max(delete_indices) < data.shape[0])

        delete_indices = set(delete_indices)

        layer = Points(data)

        layer.remove_selected()
        assert len(layer.data) == data.shape[0]

        layer.selected_data = delete_indices
        layer.remove_selected()
        assert len(layer.data) == data.shape[0] - len(delete_indices)
        assert len(layer.selected_data) == 0

        keep = list(range(data.shape[0]))
        for di in delete_indices:
            keep.remove(di)
        assert_array_equal(layer.data, data[keep])

    @given(
        random_points_st,
        st.lists(elements=st.integers(min_value=0), min_size=1),
        st.integers(min_value=0),
    )
    def test_deleting_selected_value_changes(self, data, delete_indices, value):
        """Test deleting selected points appropriately sets self._value
        Rewrite of test_deleting_selected_value_changes
        """
        assume(np.max(delete_indices) < data.shape[0])

        delete_indices = set(delete_indices)

        layer = Points(data)

        ###
        layer._value = value
        layer.selected_data = delete_indices

        if value in layer.selected_data:
            layer.remove_selected()
            assert layer._value is None
        else:
            layer.remove_selected()
            assert layer._value == value

    @given(
        random_points_2D_st,
        st.lists(elements=st.integers(min_value=0), min_size=1),
        st.tuples(
            st.floats(allow_infinity=False, min_value=-1e10, max_value=1e10),
            st.floats(allow_infinity=False, min_value=-1e10, max_value=1e10),
        ),
    )
    def test_move(self, data, moving_indices, move_offset):
        """Test deleting selected points appropriately sets self._value
        Rewrite of test_move
        """
        assume(np.max(moving_indices) < data.shape[0])

        unmoved = copy(data)
        layer = Points(data)

        move_offset = move_offset[: data.shape[1]]

        layer._move(moving_indices, list(np.zeros(data.shape[1])))
        layer._move(moving_indices, list(move_offset))
        layer._drag_start = None

        assert_almost_equal(
            layer.data[moving_indices], unmoved[moving_indices] + move_offset
        )

        unmoved_indices = list(range(data.shape[0]))
        for mi in moving_indices:
            if mi in unmoved_indices:
                unmoved_indices.remove(mi)
        assert_almost_equal(layer.data[unmoved_indices], unmoved[unmoved_indices])

    @given(random_points_st, st.text(), st.text())
    def test_name(self, data, firstname, secondname):
        """Test setting layer name.
        Rewrite of test_name
        """
        layer = Points(data)
        assert layer.name == "Points"

        layer = Points(data, name=firstname)
        assert layer.name == firstname

        layer.name = secondname
        assert layer.name == secondname

    @given(random_points_st, st.booleans())
    def test_visibility(self, data, visible):
        """Test setting layer visibility.
        Rewrite of test_visibility
        """
        layer = Points(data)
        assert layer.visible is True

        layer = Points(data, visible=visible)
        assert layer.visible is visible

        layer.visible = not visible
        assert layer.visible is not visible

    @given(
        random_points_st,
        random_opacity_st,
        random_opacity_st,
    )
    def test_opacity(self, data, opacity_first, opacity_second):
        """Test setting layer opacity.
        Rewrite of test_opacity
        """
        layer = Points(data)
        assert layer.opacity == 1.0

        layer.opacity = opacity_first
        assert layer.opacity == opacity_first

        layer = Points(data, opacity=opacity_first)
        assert layer.opacity == opacity_first

        layer.opacity = opacity_second
        assert layer.opacity == opacity_second

    @given(random_points_st, random_blending_st, random_blending_st)
    def test_blending(self, data, blending_first, blending_second):
        """Test setting layer blending.
        Rewrite of test_blending
        """
        layer = Points(data)
        assert layer.blending == "translucent"

        layer.blending = blending_first
        assert layer.blending == blending_first

        layer = Points(data, blending=blending_first)
        assert layer.blending == blending_first

        layer.blending = blending_second
        assert layer.blending == blending_second

    @given(random_points_st, random_symbol_st, random_symbol_st)
    def test_symbol(self, data, symbol_first, symbol_second):
        """Test setting symbol.
        Rewrite of test_symbol
        """
        layer = Points(data)
        assert layer.symbol == "disc"

        layer.symbol = symbol_first
        assert layer.symbol == symbol_first

        layer = Points(data, symbol=symbol_first)
        assert layer.symbol == symbol_first

        layer.symbol = symbol_second
        assert layer.symbol == symbol_second

    @given(
        random_points_st,
        st.floats(allow_nan=False, allow_infinity=False),
        st.floats(allow_nan=False, allow_infinity=False),
    )
    def test_edge_width(self, data, edge_width_first, edge_width_second):
        """Test setting edge width.
        Rewrite of test_edge_width
        """
        layer = Points(data)
        assert layer.edge_width == 1

        layer.edge_width = edge_width_first
        assert layer.edge_width == edge_width_first
        assert layer.edge_width >= 0

        layer = Points(data, edge_width=edge_width_first)
        assert layer.edge_width == edge_width_first

        layer.edge_width = edge_width_second
        assert layer.edge_width == edge_width_second

    @given(random_points_st, st.booleans())
    def test_n_dimensional(self, data, n_dimensional):
        """Test setting n_dimensional flag for 2D and 4D data.
        Rewrite of test_n_dimensional
        """
        layer = Points(data)
        assert layer.n_dimensional is False

        layer.n_dimensional = n_dimensional
        assert layer.n_dimensional is n_dimensional

        layer = Points(data, n_dimensional=n_dimensional)
        assert layer.n_dimensional == n_dimensional

        layer.n_dimensional = not n_dimensional
        assert layer.n_dimensional is not n_dimensional

    @given(random_points_st)
    def test_message(self, data):
        """Test converting value and coords to message.
        Rewrite of test_message, test_message_3d
        """
        layer = Points(data)
        msg = layer.get_status(data[0])
        assert type(msg) == str

    # Ignored Tests

    # Couldn't rewrite properly due to weird np.assert_array_equal issues.
    # @given(
    #     hynp.arrays(
    #         dtype=np.float64,
    #         shape=shared(random_shape_st, key="test_size"),
    #         elements=random_point_element_st,
    #     ),
    #     hynp.arrays(
    #         dtype=np.float64,
    #         shape=shared(random_shape_st, key="test_size"),
    #         elements=random_point_element_st,
    #     ),
    #     hynp.arrays(
    #         dtype=np.float64,
    #         shape=shared(random_shape_st, key="test_size"),
    #         elements=random_point_element_st,
    #     ),
    #     st.floats(allow_nan=False, allow_infinity=False),
    #     st.floats(allow_nan=False, allow_infinity=False),
    #     st.lists(elements=st.integers(min_value=0), min_size=1),
    # )
    # def test_size(
    #     self,
    #     data,
    #     new_points_first,
    #     new_points_second,
    #     size_first,
    #     size_second,
    #     new_size_indices,
    # ):
    #     """Test setting size with scalar.
    #     Rewrite of test_size
    #     """
    #     assume(np.max(new_size_indices) < data.shape[0])

    #     default_size = 10

    #     original_indices = list(range(data.shape[0]))
    #     new_indices_first = list(
    #         range(data.shape[0], data.shape[0] + new_points_first.shape[0])
    #     )
    #     new_indices_second = list(
    #         range(
    #             data.shape[0] + new_points_first.shape[0],
    #             data.shape[0] + new_points_first.shape[0] + new_points_second.shape[0],
    #         )
    #     )

    #     layer = Points(data)
    #     assert layer.current_size == default_size
    #     assert layer.size.shape == data.shape

    #     # Add a new point, it should get current size
    #     layer.add(new_points_first)
    #     # TODO: Fails for some reason - problem with nan comparison?
    #     # assert assert_array_equal(
    #     #     layer.size[new_indices_first],
    #     #     [[default_size] * new_points_first.ndim] * new_points_first.shape[0],
    #     # )
    #     assert np.min(layer.size) == default_size and np.max(layer.size) == default_size
    #     assert layer.data.shape[0] == data.shape[0] + new_points_first.shape[0]

    #     layer.current_size = size_first
    #     assert layer.current_size == size_first

    # @given(
    #     random_points_st, st.floats(allow_nan=False, allow_infinity=False, min_value=1)
    # )
    # def test_value(self, data, offset):
    #     """Test getting the value of the data at the current coordinates.
    #     Rewrite of test_value
    #     """
    #     layer = Points(data)
    #     index = layer.get_value(data[0])
    #     assert index == 0
