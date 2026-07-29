def test_package_import():
    import botas


def test_cli_import():
    from botas.cli.align import main

    assert callable(main)
