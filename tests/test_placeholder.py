import pytest


def test_import():
    import yt_libyt


def test_always_pass():
    assert True


@pytest.mark.xfail
def test_expected_failure():
    assert False
