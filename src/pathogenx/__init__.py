"""
Top-level module, including resource and optional dependency management.
"""
from typing import Callable, TYPE_CHECKING
from warnings import warn
from functools import wraps
from importlib import import_module
from pathlib import Path

if TYPE_CHECKING:
    from importlib.metadata import PackageMetadata

# Classes --------------------------------------------------------------------------------------------------------------
class Resources:
    """Holds global resources for this package which are generated on demand.

    This class provides a centralized way to access package-level information
    like metadata and to check for the availability of optional dependencies.
    Resources are loaded lazily upon first access.

    Attributes:
        package (str): The name of the package.
        optional_packages (set[str]): A set of available optional packages.
    """
    def __init__(self, *optional_packages: str):
        """Initializes the Resources object.

        Args:
            *optional_packages: A variable number of optional package names to
                check for availability (e.g., 'numpy', 'pandas').
        """
        self.package: str = Path(__file__).parent.name
        self.optional_packages: set[str] = set(filter(self._check_module, optional_packages))
        self._metadata: 'PackageMetadata' = None  # Generated on demand

    @property
    def metadata(self) -> 'PackageMetadata':
        """Package metadata.

        Lazily loads and returns the package's metadata (e.g., version).

        Returns:
            PackageMetadata: The metadata object for the package.
        """
        if self._metadata is None:
            from importlib.metadata import metadata
            self._metadata = metadata(self.package)
        return self._metadata

    @staticmethod
    def _check_module(module_name: str) -> bool:
        """Checks if a module can be imported.

        Args:
            module_name (str): The name of the module to check.

        Returns:
            bool: True if the module can be imported, False otherwise.
        """
        try:
            import_module(module_name)
            return True
        except ImportError:
            return False



class PathogenxWarning(Warning):
    """Base warning class for the pathogenx package.

    This allows users to easily silence all warnings from this package.

    Examples:
        To ignore all warnings from this package:
        >>> import warnings
        >>> from pathogenx import PathogenxWarning
        ... warnings.simplefilter('ignore', PathogenxWarning)
    """

    pass


class DependencyWarning(PathogenxWarning):
    """Warning issued when an optional dependency is not found."""
    pass


def require(*packages: str) -> Callable:
    """A decorator to verify optional dependencies before function execution.

    If any of the specified packages are not installed, this decorator will
    issue a `DependencyWarning` and prevent the decorated function from
    running.

    Args:
        *packages: The names of optional packages required by the function.

    Returns:
        Callable: The wrapped function, which will only execute if all
            dependencies are met.

    Examples:
        >>> from pathogenx import require
        ...
        ... @require('numpy', 'pandas')
        ... def process_data_with_numpy_and_pandas(data):
        ...     # This code will only run if numpy and pandas are installed
        ...     pass
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            if missing_deps := [dep for dep in packages if dep not in RESOURCES.optional_packages]:
                warn(
                    f"Function '{func.__name__}' requires the following missing dependencies: "
                    f"{', '.join(missing_deps)}. Skipping execution.",
                    DependencyWarning
                )
                return None
            return func(*args, **kwargs)
        return wrapper
    return decorator


# Constants ------------------------------------------------------------------------------------------------------------
RESOURCES = Resources('pathogenx.app')
