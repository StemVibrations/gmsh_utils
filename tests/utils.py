from typing import Dict, Any

import numpy.testing as npt
class TestUtils:

    @staticmethod
    def assert_dictionary_almost_equal(expected: Dict[Any, Any], actual: Dict[Any, Any]):
        """
        Checks whether two dictionaries are equal.

        Args:
            expected: Expected dictionary.
            actual: Actual dictionary.

        """

        for k, v in expected.items():
            assert k in actual

            if isinstance(v, dict):
                TestUtils.assert_dictionary_almost_equal(v, actual[k])
            elif isinstance(v, str):
                assert v == actual[k]
            else:
                npt.assert_allclose(v, actual[k])

