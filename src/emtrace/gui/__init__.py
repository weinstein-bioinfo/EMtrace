"""EMtrace desktop GUI (PyQt6).

Public entry point: :func:`emtrace.gui.main_window.launch_gui`.

Nothing is re-exported at the package level so that touching this package
does not force a hard dependency on the optional ``[gui]`` extra (PyQt6).
The CLI imports ``launch_gui`` lazily and falls back to a friendly error
when PyQt6 is missing.
"""
