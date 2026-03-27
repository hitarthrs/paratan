"""
MachineComponent — base class for all machine component builders.

Every component in the machine (VV, CC, coils, end cells, and any future
additions) inherits from this and implements the same interface. This is
what lets SimpleMachineBuilder treat all components uniformly.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
import openmc


class MachineComponent(ABC):
    """
    Abstract base class for a machine component builder.

    Subclasses hold their config and material namespace, implement build(),
    and expose consistent outputs so the main builder and costing module
    can query them without knowing the internals.
    """

    @abstractmethod
    def build(self, context) -> list[openmc.Cell]:
        """
        Construct OpenMC cells for this component.

        Parameters
        ----------
        context : BuildContext
            Shared build state. Read what you need; write back what
            downstream components will need.

        Returns
        -------
        list[openmc.Cell]
            All cells belonging to this component.
        """
        ...

    @abstractmethod
    def get_tally_descriptors(self) -> list[dict]:
        """
        Return tally descriptor dicts for this component.

        Each dict has keys: type, location, description, cell,
        cell_tallies, mesh_tallies — matching the TallyBuilder interface.
        """
        ...

    @property
    def outermost_region(self) -> openmc.Region | None:
        """
        The outermost region of this component after build().

        Used by the builder to subtract this component from the room and
        from subsequent components. Return None if not applicable.
        """
        return None

    @property
    def volume_estimate(self) -> float:
        """
        Approximate volume of this component in cm³.

        Used by the costing module. Subclasses should override with a
        proper geometric calculation. Default returns 0.
        """
        return 0.0
