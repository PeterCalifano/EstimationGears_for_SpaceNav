"""Smoke tests for the EstimationGears Python package."""

from __future__ import annotations

import EstimationGears_for_SpaceNav


class TestPythonSmoke:
    """Verify the source package exposes wrapper availability state."""

    def test_import_exposes_wrapper_status(self) -> None:
        """Import the package and inspect its public wrapper status fields."""

        bHasWrapper_ = EstimationGears_for_SpaceNav.HAS_WRAPPER

        assert isinstance(bHasWrapper_, bool)
        assert hasattr(EstimationGears_for_SpaceNav, "WRAPPER_IMPORT_ERROR")
