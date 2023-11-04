import pytest


def test_import():
    import yt_libyt
    api_list = ["libytDataset", "libytGrid", "libytHierarchy", "libytFieldInfo", "libytIOHandler"]
    for attr in api_list:
        assert hasattr(yt_libyt, attr) is True, "No class named %s in yt_libyt" % attr


def test_always_pass():
    assert True


@pytest.mark.xfail
def test_expected_failure():
    assert 1 == 2
