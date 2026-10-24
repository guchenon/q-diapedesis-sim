"""
Modules B2 & B5: sub-voxel MRI blind-window observer and post-task degaussing + retrograde
wall-rolling retrieval package.

Naming note (D11, ELN #003): `PINOSubVoxelObserver` is a Physics-Informed Spectral MAP
Observer - a five-mode regularized Fourier estimator of the cerebral haemodynamic waveform.
There is no neural network in it and no operator is learned. `PhysicsInformedSpectralObserver`
is the accurate alias and both names refer to the same class.
"""

from .pino_observer_retrieval import (
    PINOSubVoxelObserver,
    PhysicsInformedSpectralObserver,
    SubVoxelObserverReport,
    DegaussingRetrogradeRetrieval,
    RetrievalMassBalanceReport,
)

__all__ = [
    "PINOSubVoxelObserver",
    "PhysicsInformedSpectralObserver",
    "SubVoxelObserverReport",
    "DegaussingRetrogradeRetrieval",
    "RetrievalMassBalanceReport",
]
