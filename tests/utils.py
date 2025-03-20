from typing import Dict, Any
import json

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

    @staticmethod
    def read_mesh_data_from_json(path: str) -> Dict[str, Any]:
        """
        Read mesh data from a JSON file and makes sure all node and element keys are integers

        Args:
            - path (str): Path to the JSON file.

        Returns:
            - Dict[str, Any]: Mesh data.
        """

        with open(path, "r") as file:
            mesh_data = json.load(file)

        # make sure all node and element keys are integers
        mesh_data["nodes"] = {int(k): v for k, v in mesh_data["nodes"].items()}
        mesh_data["elements"] = {k: {int(kk): vv for kk, vv in v.items()} for k, v in mesh_data["elements"].items()}

        return mesh_data
