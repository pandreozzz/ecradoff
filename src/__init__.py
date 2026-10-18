"""Package initializer for :mod:`ecradoff.src`.

Expose immediate subpackages as attributes on first access, e.g.
``ecradoff.src.interp``.
"""

from importlib import import_module as _import_module
from pkgutil import iter_modules as _iter_modules


def __getattr__(name: str):
	module_name = f"{__name__}.{name}"
	try:
		module = _import_module(module_name)
	except ModuleNotFoundError as exc:
		if exc.name == module_name:
			raise AttributeError(f"module '{__name__}' has no attribute '{name}'") from None
		raise

	globals()[name] = module
	return module


def __dir__():
	names = set(globals())
	names.update(module.name for module in _iter_modules(__path__))
	return sorted(names)

