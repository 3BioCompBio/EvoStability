
# Imports ----------------------------------------------------------------------
from typing import List, Self
import numpy as np

# Single Bin values contains ---------------------------------------------------
class Bin():

    # Constructor --------------------------------------------------------------
    def __init__(self, min_value: float, max_value: float):
        self.min = min_value
        self.max = max_value
        self.center = (min_value + max_value) / 2.0
        self.delta = max_value - min_value
        self.values = []
        assert self.delta > 0, f"ERROR in Bin(): max_value (={max_value}) should be > than min_value (={min_value})."

    # Base properties ----------------------------------------------------------
    def __len__(self) -> int:
        return len(self.values)
    
    def __str__(self) -> str:
        return f"Bin([{self.min:.3f}, {self.max:.3f}], l={len(self)})"
    
    def __getitem__(self, id: int) -> float:
        return self.values[id]

    def __iter__(self):
        return iter(self.values)
    
    # Methods ------------------------------------------------------------------
    def append(self, value: float) -> Self:
        self.values.append(value)
        return self

    def add_list(self, values_list: List[float]) -> Self:
        for value in values_list:
            self.values.append(value)
        return self

    def sort(self) -> Self:
        self.values.sort()
        return self

    def clear(self) -> Self:
        self.values = []
        return self


# Main -------------------------------------------------------------------------
class Bins():

    # Constructor --------------------------------------------------------------
    def __init__(self, x_min: float, x_max: float, n_steps: int):
        
        # Guardians
        assert x_min < x_max, f"ERROR in Bins(): x_min (={x_min}) should be < than x_max (={x_max})."
        assert n_steps >= 1, f"ERROR in Bins(): n_steps (={n_steps}) should be at least 1."

        # Set base properties
        self.x_min = x_min
        self.x_max = x_max
        self.n_steps = n_steps
        self.x_delta = self.x_max - self.x_min
        self.step_size = self.x_delta / n_steps
        self.sorted = False

        # Set dependencies properties
        self._scale_range = self.n_steps / self.x_delta
        self._max_id = self.n_steps-1

        # Init bins
        self.bins = []
        current_values = self.x_min
        for i in range(self.n_steps):
            self.bins.append(Bin(current_values, current_values + self.step_size))
            current_values += self.step_size

    def sort(self) -> Self:
        for bin in self.bins:
            bin.sort()
        self.sorted = True

    # Base properties ----------------------------------------------------------
    def __len__(self) -> int:
        return self.n_steps
    
    def __str__(self) -> str:
        return f"Bins([{self.x_min:.3f}, {self.x_max:.3f}], l={len(self)})"
    
    def __getitem__(self, id: int) -> Bin:
        return self.bins[id]

    def __iter__(self):
        return iter(self.bins)

    def __contains__(self, x_value: float) -> bool:
        return self.x_min <= x_value <= self.x_max
    
    def total(self) -> int:
        return sum([len(bin) for bin in self.bins])
    
    def length(self) -> List[int]:
        return [len(bin) for bin in self]
    
    # Methods ------------------------------------------------------------------
    def add_measure(self, x: float, y: float) -> Self:
        assert x in self, f"ERROR in {self}: x={x} is not in x_range."
        self.sorted = False
        x_scaled = (x - self.x_min) * self._scale_range
        x_id = min(int(x_scaled), self._max_id)
        self.bins[x_id].append(y)
        return self
    
    def add_measures_list(self, x_list: List[float], y_list: List[float]) -> Self:
        for x0, y0 in zip(x_list, y_list):
            self.add_measure(x0, y0)
        return self
    
    def get_bins(self, bin_count_thr: int=1):
        return [bin for bin in self.bins if len(bin) >= bin_count_thr]
    
    def get_center_list(self, bin_count_thr: int=1) -> List[float]:
        bins = self.get_bins(bin_count_thr)
        return np.array([bin.center for bin in bins])
    
    def get_mean_list(self, bin_count_thr: int=1) -> List[float]:
        bins = self.get_bins(bin_count_thr)
        return np.array([np.mean(bin.values) for bin in bins])
    
    def get_std_list(self, bin_count_thr: int=1) -> List[float]:
        bins = self.get_bins(bin_count_thr)
        return np.array([np.std(bin.values) for bin in bins])
    
    def get_percentile_list(self, p: float, bin_count_thr: int=1):
        assert self.sorted, f"ERROR in {self}.get_percentile_list() first Bins should be sorted."
        bins = self.get_bins(bin_count_thr)
        ids = [int(p*len(bin)) for bin in bins]
        return np.array([bin.values[id] for bin, id in zip(bins, ids)])