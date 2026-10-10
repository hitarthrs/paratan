"""Paratan viewer platform. Heavy visualization dependencies are imported lazily."""
__all__ = ['build_from_yaml', 'build_simple_mirror_components']


def __getattr__(name):
    if name in __all__:
        from src.paratan.viewer import simple_mirror_meshes
        return getattr(simple_mirror_meshes, name)
    raise AttributeError(name)
