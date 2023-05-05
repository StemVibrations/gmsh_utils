import numpy.testing as npt
class TestUtils:

    @staticmethod
    def assert_dictionary_almost_equal(expected, actual):

        a=1+1
        for k,v in expected.items():
            assert k in actual

            if isinstance(v, dict):
                TestUtils.assert_dictionary_almost_equal(v, actual[k])
            else:
                npt.assert_allclose(v, actual[k])

