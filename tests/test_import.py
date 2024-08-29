def test_import():
    import yt_libyt

    api_list = ["libytDataset", "libytGrid", "libytHierarchy", "libytFieldInfo", "libytIOHandler"]
    for attr in api_list:
        assert hasattr(yt_libyt, attr) is True, f"No class named {attr} in yt_libyt"
