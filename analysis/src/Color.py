
# Imports ----------------------------------------------------------------------
from typing import Union, Self, Tuple, List
import random
import colorsys
import numpy as np
import matplotlib.pyplot as plt

# Main -------------------------------------------------------------------------
class Color():
    """
    Color Contained Class to generate and modify colors. Values are stored as RGB floats.
        
        
    Usage:

        * Input parameters are RGB values as (1) float in [0.0, 1.0] or (2) integers in [0, 255]
            $ color = Color(0.8, 0.1, 0.3)
        
        * Can output RGB (as float or int), HSL (float), HSV (float)
            $ r, g, b = color.rgb
            $ rint, gint, bint = color.rgb_255
            $ h, s, l = color.hsl()
        
        * Give updated color by variation in R, G, B, H, S, L, V
            $ blue = Color(0.0, 0.0, 1.0)
            $ dark_blue = blue.updated_lightness(-0.3)
        
        * Interpolates colors in N steps from 2 colors
            $ colors_blue_red_range = Color.get_interpolated_range(blue, red, 15)
        
        * Generate R, G, B, H, S, L, V color ranges from an initial colors
            $ red_lightness_range = Color.get_lightness_range(red, 15)
    """

    # Constructor --------------------------------------------------------------
    def __init__(self, red: Union[int, float], green: Union[int, float, None]=None, blue: Union[int, float, None]=None, color_space: Union[str, None]=None):

        # Accept single value input -> red = green = blue
        if green is None or blue is None:
            assert green is None and blue is None, f"ERROR in Color(): Input RGB values should be given with 1 (greyscale) or 3 (RGB) values."
            green = red
            blue = red

        # Convert to RGB if color_space is 'HSL' or 'HSV'
        if color_space is not None:
            if color_space.upper() == "HSL":
                red, green, blue = self.hsl_to_rgb(red, green, blue)
            elif color_space.upper() == "HSV":
                red, green, blue = self.hsv_to_rgb(red, green, blue)
            else:
                raise ValueError(f"ERROR in Color(): if specified, color_space='{color_space}' should be HSL or HSV.")

        # Guardiant by input type
        if isinstance(red, float):
            assert isinstance(green, float) and isinstance(blue, float), f"ERROR in Color(): input values RGB={(red, green, blue)} should all be float in [0.0, 1.0] or int in [0, 255]."
            assert 0.0 <= red <= 1.0 and 0.0 <= green <= 1.0 and 0.0 <= blue <= 1.0, f"ERROR in Color(): input values RGB={(red, green, blue)} should all be float in [0.0, 1.0] or int in [0, 255]."
        elif isinstance(red, int):
            assert isinstance(green, int) and isinstance(blue, int), f"ERROR in Color(): input values RGB={(red, green, blue)} should all be float in [0.0, 1.0] or int in [0, 255]."
            assert 0 <= red <= 255 and 0 <= green <= 255 and 0 <= blue <= 255, f"ERROR in Color(): input values RGB={(red, green, blue)} should all be float in [0.0, 1.0] or int in [0, 255]."
            red = red / 255.0
            green = green / 255.0
            blue = blue / 255.0
        else:
            raise ValueError(f"ERROR in Color(): input values RGB={(red, green, blue)} should all be float in [0.0, 1.0] or int in [0, 255].")

        # Assign
        self._rgb = (red, green, blue)

    def copy(self) -> Self:
        return Color(self.red, self.green, self.blue)

    # Base values --------------------------------------------------------------
    @property
    def red(self) -> float:
        return self._rgb[0]
    
    @property
    def green(self) -> float:
        return self._rgb[1]
    
    @property
    def blue(self) -> float:
        return self._rgb[2]
    
    @property
    def rgb(self) -> Tuple[float, float, float]:
        return (self.red, self.green, self.blue)
    
    @property
    def red_255(self) -> int:
        return Color._get_255(self.red)
    
    @property
    def green_255(self) -> int:
        return Color._get_255(self.green)
    
    @property
    def blue_255(self) -> int:
        return Color._get_255(self.blue)
    
    @property
    def rgb_255(self) -> Tuple[int, int, int]:
        return (self.red_255, self.green_255, self.blue_255)
    
    def __len__(self) -> int:
        return len(self._rgb)
    
    def __getitem__(self, id: int) -> float:
        return self._rgb[id]
    
    def __iter__(self):
        return iter(self._rgb)
    
    def hsl(self) -> Tuple[float, float, float]:
        return Color.rgb_to_hsl(self.red, self.green, self.blue)
    
    def hsv(self) -> Tuple[float, float, float]:
        return Color.rgb_to_hsv(self.red, self.green, self.blue)
    
    def hue(self) -> float:
        return self.hsl()[0]
    
    def saturation(self) -> float:
        return self.hsl()[1]
    
    def lightness(self) -> float:
        return self.hsl()[2]
    
    def value(self) -> float:
        return self.hsv()[2]

    # Base methods -------------------------------------------------------------
    def __str__(self) -> str:
        return f"Color(RGB=[{self.red_255}, {self.green_255}, {self.blue_255}])"
    
    def __eq__(self, other: Self) -> bool:
        return self.rgb_255 == other.rgb_255
    
    def __add__(self, other: Self) -> Self:
        """Add two colors as RGB vectors."""
        return Color(self.red + other.red, self.green + other.green, self.blue + other.blue)
    
    def __mul__(self, l: float) -> Self:
        """Multiply a color as RGB vector by scalar."""
        return Color(self.red*l, self.green*l, self.blue*l)
    
    def __rmul__(self, l: float) -> Self:
        """Multiply a color as RGB vector by scalar."""
        return self.__mul__(l)  # To enable right multiplication by scalar
    
    # Update methods -----------------------------------------------------------
    def updated_red(self, update: float) -> Self:
        r, g, b = self.rgb
        return Color(Color._bound_float(r+update), g, b)

    def updated_green(self, update: float) -> Self:
        r, g, b = self.rgb
        return Color(r, Color._bound_float(g+update), b)
    
    def updated_blue(self, update: float) -> Self:
        r, g, b = self.rgb
        return Color(r, g, Color._bound_float(b+update))
    
    def updated_hue(self, update: float) -> Self:
        h, s, l = self.hsl()
        h = Color._bound_float(h+update)
        r, g, b = Color.hsl_to_rgb(h, s, l)
        return Color(r, g, b)
    
    def updated_saturation(self, update: float) -> Self:
        h, s, l = self.hsl()
        s = Color._bound_float(s+update)
        r, g, b = Color.hsl_to_rgb(h, s, l)
        return Color(r, g, b)
    
    def updated_lightness(self, update: float) -> Self:
        h, s, l = self.hsl()
        l = Color._bound_float(l+update)
        r, g, b = Color.hsl_to_rgb(h, s, l)
        return Color(r, g, b)
    
    def updated_value(self, update: float) -> Self:
        h, s, v = self.hsv()
        v = Color._bound_float(v+update)
        r, g, b = Color.hsv_to_rgb(h, s, v)
        return Color(r, g, b)
    
    # Color generations --------------------------------------------------------

    @classmethod
    def random(cls) -> Self:
        return Color(random.random(), random.random(), random.random())
    
    @classmethod
    def random_list(cls, n: int) -> List[Self]:
        return [cls.random() for _ in range(n)]

    @classmethod
    def interpolate(cls, color1: Self, color2: Self, lambda1: float=0.5, color_space: Union[None, str]=None) -> Self:
        "Interpolate two colors in RGB space at l ratio between color1 and color2."
        assert 0.0 <= lambda1 <= 1.0, f"ERROR in Color().interpolate(): Interpolation by float={lambda1} only works for l in [0.0, 1.0]."

        # Manage HSL and HSV
        if color_space is not None:
            if color_space.upper() == "HSL":
                h1, s1, l1 = color1.hsl()
                h2, s2, l2 = color2.hsl()
                h, s, l = (1.0-lambda1)*h1 + lambda1*h2, (1.0-lambda1)*s1 + lambda1*s2, (1.0-lambda1)*l1 + lambda1*l2
                return Color(h, s, l, color_space=color_space)
            elif color_space.upper() == "HSV":
                h1, s1, v1 = color1.hsl()
                h2, s2, v2 = color2.hsl()
                h, s, v = (1.0-lambda1)*h1 + lambda1*h2, (1.0-lambda1)*s1 + lambda1*s2, (1.0-lambda1)*v1 + lambda1*v2
                return Color(h, s, v, color_space=color_space)
            else:
                raise ValueError(f"ERROR in Color.interpolate(): color_space='{color_space}' is specified should be HSL or HSV.")

        # Base case
        color_arr = (1.0-lambda1)*color1 + lambda1*color2
        return Color(color_arr[0], color_arr[1], color_arr[2])
    
    @classmethod
    def get_interpolated_range(cls, color1: Self, color2: Self, n: int, color_space: Union[None, str]=None) -> List[Self]:
        """Get a list of n interpolated colors between color1 and color2."""

        # Guardians
        assert n > 0, f"ERROR in Color.get_interpolated_range(): number of colors n={n} should be > 1."

        # Edge case: only 1 color in the middle
        if n == 1:
            return [cls.interpolate(color1, color2)]

        # Init
        step = 1.0 / (n-1)

        # Generate
        colors_arr = []
        for i in range(n):
            l = float(i) * step
            new_color = Color.interpolate(color1, color2, l, color_space=color_space)
            colors_arr.append(new_color)
        return colors_arr
    
    @classmethod
    def get_red_range(cls, color: Self, n: int, delta: Union[None, float]=None, direction: int=0) -> List[Self]:

        # Guardians
        assert n > 0, f"ERROR in Color.get_red_range(): number of colors n={n} should be > 1."

        # Edge case: only 1 color -> return input color in a list
        if n == 1:
            return [color.copy()]
        
        # Init
        r, g, b = color.rgb

        # Solve range
        x1, x2 = Color._get_range(r, direction, delta)
        delta = x2 - x1
        step = delta / (n-1)

        # Generate
        colors_arr = []
        for i in range(n):
            lambda1 = float(i) * step
            ri = x1 + lambda1
            new_color = Color(ri, g, b)
            colors_arr.append(new_color)
        return colors_arr
    
    @classmethod
    def get_green_range(cls, color: Self, n: int, delta: Union[None, float]=None, direction: int=0) -> List[Self]:

        # Guardians
        assert n > 0, f"ERROR in Color.get_hue_range(): number of colors n={n} should be > 1."

        # Edge case: only 1 color -> return input color in a list
        if n == 1:
            return [color.copy()]
        
        # Init
        r, g, b = color.rgb

        # Solve range
        x1, x2 = Color._get_range(g, direction, delta)
        delta = x2 - x1
        step = delta / (n-1)

        # Generate
        colors_arr = []
        for i in range(n):
            lambda1 = float(i) * step
            gi = x1 + lambda1
            new_color = Color(r, gi, b)
            colors_arr.append(new_color)
        return colors_arr
    
    @classmethod
    def get_blue_range(cls, color: Self, n: int, delta: Union[None, float]=None, direction: int=0) -> List[Self]:

        # Guardians
        assert n > 0, f"ERROR in Color.get_hue_range(): number of colors n={n} should be > 1."

        # Edge case: only 1 color -> return input color in a list
        if n == 1:
            return [color.copy()]
        
        # Init
        r, g, b = color.rgb

        # Solve range
        x1, x2 = Color._get_range(b, direction, delta)
        delta = x2 - x1
        step = delta / (n-1)

        # Generate
        colors_arr = []
        for i in range(n):
            lambda1 = float(i) * step
            bi = x1 + lambda1
            new_color = Color(r, g, bi)
            colors_arr.append(new_color)
        return colors_arr
    
    @classmethod
    def get_hue_range(cls, color: Self, n: int, delta: Union[None, float]=None, direction: int=0) -> List[Self]:

        # Guardians
        assert n > 0, f"ERROR in Color.get_hue_range(): number of colors n={n} should be > 1."

        # Edge case: only 1 color -> return input color in a list
        if n == 1:
            return [color.copy()]
        
        # Init
        h, s, l = color.hsl()

        # Solve range
        x1, x2 = Color._get_range(h, direction, delta)
        delta = x2 - x1
        step = delta / (n-1)

        # Generate
        colors_arr = []
        for i in range(n):
            lambda1 = float(i) * step
            hi = x1 + lambda1
            ri, gi, bi = Color.hsl_to_rgb(hi, s, l)
            new_color = Color(ri, gi, bi)
            colors_arr.append(new_color)
        return colors_arr
    
    @classmethod
    def get_saturation_range(cls, color: Self, n: int, delta: Union[None, float]=None, direction: int=0) -> List[Self]:

        # Guardians
        assert n > 0, f"ERROR in Color.get_saturation_range(): number of colors n={n} should be > 1."

        # Edge case: only 1 color -> return input color in a list
        if n == 1:
            return [color.copy()]
        
        # Init
        h, s, l = color.hsl()

        # Solve range
        x1, x2 = Color._get_range(s, direction, delta)
        delta = x2 - x1
        step = delta / (n-1)

        # Generate
        colors_arr = []
        for i in range(n):
            lambda1 = float(i) * step
            si = x1 + lambda1
            ri, gi, bi = Color.hsl_to_rgb(h, si, l)
            new_color = Color(ri, gi, bi)
            colors_arr.append(new_color)
        return colors_arr
    
    @classmethod
    def get_lightness_range(cls, color: Self, n: int, delta: Union[None, float]=None, direction: int=0) -> List[Self]:

        # Guardians
        assert n > 0, f"ERROR in Color.get_lightness_range(): number of colors n={n} should be > 1."

        # Edge case: only 1 color -> return input color in a list
        if n == 1:
            return [color.copy()]
        
        # Init
        h, s, l = color.hsl()

        # Solve range
        x1, x2 = Color._get_range(l, direction, delta)
        delta = x2 - x1
        step = delta / (n-1)

        # Generate
        colors_arr = []
        for i in range(n):
            lambda1 = float(i) * step
            li = x1 + lambda1
            ri, gi, bi = Color.hsl_to_rgb(h, s, li)
            new_color = Color(ri, gi, bi)
            colors_arr.append(new_color)
        return colors_arr
    
    @classmethod
    def get_value_range(cls, color: Self, n: int, delta: Union[None, float]=None, direction: int=0) -> List[Self]:

        # Guardians
        assert n > 0, f"ERROR in Color.get_lightness_range(): number of colors n={n} should be > 1."

        # Edge case: only 1 color -> return input color in a list
        if n == 1:
            return [color.copy()]
        
        # Init
        h, s, v = color.hsv()

        # Solve range
        x1, x2 = Color._get_range(v, direction, delta)
        delta = x2 - x1
        step = delta / (n-1)

        # Generate
        colors_arr = []
        for i in range(n):
            lambda1 = float(i) * step
            vi = x1 + lambda1
            ri, gi, bi = Color.hsl_to_rgb(h, s, vi)
            new_color = Color(ri, gi, bi)
            colors_arr.append(new_color)
        return colors_arr
    
    # Convertions --------------------------------------------------------------
    @staticmethod
    def rgb_to_hsl(r: float, g: float, b: float) -> Tuple[float, float, float]:
        h, l, s = colorsys.rgb_to_hls(r, g, b) # Here we invert the order (python uses old standard ordering)
        return (h, s, l)
    
    @staticmethod
    def hsl_to_rgb(h: float, s: float, l: float) -> Tuple[float, float, float]:
        return colorsys.hls_to_rgb(h, l, s) # Here we invert the order (python uses old standard ordering)
    
    @staticmethod
    def rgb_to_hsv(r: float, g: float, b: float) -> Tuple[float, float, float]:
        return colorsys.rgb_to_hsv(r, g, b)
    
    @staticmethod
    def hsv_to_rgb(h: float, s: float, v: float) -> Tuple[float, float, float]:
        return colorsys.hsv_to_rgb(h, s, v)
    
    # Show color or color list -------------------------------------------------
    def show(self) -> None:
        plt.bar(0.0, 1.0, color=self.rgb)
        #plt.scatter([0.0], [0.0], color=self.rgb)
        plt.ylim(0.0, 1.0)
        plt.xlim(-0.1, 0.1)
        plt.xticks([], [])
        plt.yticks([], [])
        plt.show()

    @staticmethod
    def show_list(colors_list: List[Self]) -> None:
        plt.bar([i for i, _ in enumerate(colors_list)], [1.0 for _ in colors_list], color=[c.rgb for c in colors_list])
        plt.ylim(0.0, 1.0)
        plt.xticks([], [])
        plt.yticks([], [])
        plt.show()

    # Dependencies -------------------------------------------------------------
    @staticmethod
    def _get_255(value: float) -> int:
        return round(value * 255.0)
    
    @staticmethod
    def _bound_float(value: float) -> float:
        """Return values bounded in [0, 1]"""
        if 0.0 <= value <= 1.0:
            return value
        else:
            print(f"WARNING in Color._bound_float(): value={value} bounted to interval [0, 1].")
            return max(0.0, min(1.0, value))
    
    @staticmethod
    def _get_range(center: float, direction: int, delta: Union[None, float]):
        """
        Get range (contained in [0, 1]) around center of amplitude delta.
            * If range is too big, it will be cropped.
            * If delta is not specified, maximal possible delta will be taken
            * Direction indicated:
                 0: for symmetric range
                -1: for negative range from the center
                +1: for positive range from the center
        """

        if direction == 0:
            maximal_delta = 2.0 * min(center, 1.0 - center)
            if delta is None: delta = maximal_delta
            if delta > maximal_delta:
                delta = maximal_delta
                print(f"WARNING in Color._get_range(): range of direction={direction} with delta={delta:.3f} around center={center:.3f} goes ourside alowed interval [0, 1] -> range is cropped.")
            if delta == 0.0: print(f"WARNING in Color._get_range(): obtained range in collapsed for direction={direction} around center={center:.3f}.")
            return (center - delta/2.0, center + delta/2.0)
        elif direction == -1:
            maximal_delta = center
            if delta is None: delta = maximal_delta
            if delta > maximal_delta:
                delta = maximal_delta
                print(f"WARNING in Color._get_range(): range of direction={direction} with delta={delta:.3f} around center={center:.3f} goes ourside alowed interval [0, 1] -> range is cropped.")
            if delta == 0.0: print(f"WARNING in Color._get_range(): obtained range in collapsed for direction={direction} around center={center:.3f}.")
            return (center - delta, center)
        elif direction == 1:
            maximal_delta = 1.0 - center
            if delta is None: delta = maximal_delta
            if delta > maximal_delta:
                delta = maximal_delta
                print(f"WARNING in Color._get_range(): range of direction={direction} with delta={delta:.3f} around center={center:.3f} goes ourside alowed interval [0, 1] -> range is cropped.")
            if delta == 0.0: print(f"WARNING in Color._get_range(): obtained range in collapsed for direction={direction} around center={center:.3f}.")
            return (center, center + delta)
        else:
            raise ValueError(f"ERROR in Color()._get_range(): direction={direction} should be 0 (for symmetric), -1 (for negative) or 1 (for poitive).")



# Color Gradient Generator -----------------------------------------------------
class ColorGradientSequential():
    """
    Color Gradient (Sequential: for unsigned values) mapping [0, 1] -> [color1, color2]
        * Values outside range will be bounded to [0, 1]

    Usage:
        cg = ColorGradientSequential(Color(0.8, 0.1, 0.3), Color(0.1, 0.3, 0.9))
        my_color = cg.get_color(0.37)
    """

    # Constructor --------------------------------------------------------------
    def __init__(self, color1: Color, color2: Color):
        self.color1 = color1
        self.color2 = color2

    # Methods ------------------------------------------------------------------
    def __str__(self) -> str:
        return f"ColorGradientSequential: [0.0, 1.0] -> [{self.color1}, {self.color2}]"
    
    def inversed(self) -> Self:
        return ColorGradientSequential(self.color2, self.color1)

    def get_color(self, value: float, print_warnings: bool=True) -> Color:
        if not 0.0 <= value <= 1.0:
            if print_warnings:
                print(f"WARNING in ColorGradientSequential().get_color(): values outsize range [0, 1] -> bounded to interval.")
            value = max(0.0, min(1.0, value))
        return Color.interpolate(self.color1, self.color2, value)
    
    def get_range(self, len_range: int) -> List[Color]:
        if len_range <= 0:
            return []
        if len_range == 1:
            return [self.get_color(0.5)]
        step = 1.0 / (len_range - 1)
        float_range = np.arange(0.0, 1.0 + step, step)
        return [self.get_color(r) for r in float_range]
    
    def show(self) -> None:
        N, SIZE = 200, 200.0
        x = np.random.uniform(0.0, 1.0, N)
        colors = [self.get_color(xi) for xi in x]
        plt.scatter(x, np.zeros(N), color=[c.rgb for c in colors], s=SIZE)
        plt.xticks([], [])
        plt.yticks([], [])
        plt.show()
    
class ColorGradientDivergeant():
    """
    Color Gradient (Divergeant: for signed values) mapping [-1, 0, 1] -> [color1, color2, color3]
        * Values outside range will be bounded to [-1, 1]

    Usage:
        cg = ColorGradientDivergeant(Color(0.8, 0.1, 0.3), Color(1.0, 1.0, 1.0), Color(0.1, 0.3, 0.9))
        my_color = cg.get_color(-0.37)
    """

    # Constructor --------------------------------------------------------------
    def __init__(self, color1: Color, color2: Color, color3: Color):
        self.color1 = color1
        self.color2 = color2
        self.color3 = color3

    # Methods ------------------------------------------------------------------
    def __str__(self) -> str:
        return f"ColorGradientDivergeant: [-1.0, 0.0, +1.0] -> [{self.color1}, {self.color2}, {self.color3}]"
    
    def inversed(self) -> Self:
        return ColorGradientDivergeant(self.color3, self.color2, self.color1)
    
    def get_color(self, value: float, print_warnings: bool=True) -> Color:
        if not -1.0 <= value <= 1.0:
            if print_warnings:
                print(f"WARNING in ColorGradientDivergeant().get_color(): values outsize range [-1, 1] -> bounded to interval.")
            value = max(-1.0, min(1.0, value))
        if value < 0.0:
            return Color.interpolate(self.color1, self.color2, value + 1.0)
        else:
            return Color.interpolate(self.color2, self.color3, value)
        
    def get_range(self, len_range: int) -> List[Color]:
        if len_range <= 0:
            return []
        if len_range == 1:
            return [self.get_color(0.0)]
        step = 2.0 / (len_range - 1)
        float_range = np.arange(-1.0, 1.0 + step, step)
        return [self.get_color(r) for r in float_range]
        
    def show(self) -> None:
        N, SIZE = 200, 200.0
        x = np.random.uniform(-1.0, 1.0, N)
        colors = [self.get_color(xi) for xi in x]
        plt.scatter(x, np.zeros(N), color=[c.rgb for c in colors], s=SIZE)
        plt.xticks([], [])
        plt.yticks([], [])
        plt.show()
    


# COLORS, PALETTES AND SCALES --------------------------------------------------

class COLORS:
    """Colors Library."""

    # Base
    RED = Color(1.0, 0.0, 0.0)
    GREEN = Color(0.0, 1.0, 0.0)
    BLUE = Color(0.0, 0.0, 1.0)
    WHITE = Color(1.0, 1.0, 1.0)
    BLACK = Color(0.0, 0.0, 0.0)
    GREY = Color(0.5)
    LIGHTGREY = Color(0.75)
    DARKGREY = Color(0.25)

    # Custom
    BLUE_NICE_1 = Color(0/255, 95/255, 115/255)
    BLUE_NICE_2 = Color(0/255, 76/255, 146/255)
    ORANGE_NICE_1 = Color(238/255, 155/255, 0/255)
    ORANGE_NICE_2 = Color(251/255, 133/255, 0/255)
    
    BLUE_MILD = Color(0.039, 0.314, 0.631)
    RED_MILD = Color(0.882, 0.000, 0.047)

class PALETTES:
    """Colors Palettes Library."""

    # 2/3 Colors
    COLORS2_MILD_HEATSCALE = [COLORS.RED_MILD, COLORS.BLUE_MILD]
    COLORS3_1 = [COLORS.BLUE_NICE_1, COLORS.ORANGE_NICE_2, Color(115/255,  10/255, 115/255)]
    COLORS3_2 = [COLORS.BLUE_NICE_1, COLORS.ORANGE_NICE_2, Color(142/255, 202/255, 230/255)]

    # Multiple colors
    # ...

class COLOR_GRADIENTS:
    """Colors Gradiends Library."""

    # Sequentials
    SEQUENTIAL_BLUE_ORANGE_1 = ColorGradientSequential(COLORS.BLUE_NICE_1, COLORS.ORANGE_NICE_1)
    SEQUENTIAL_BLUE_ORANGE_2 = ColorGradientSequential(COLORS.BLUE_NICE_2, COLORS.ORANGE_NICE_2)
    SEQUENTIAL_RED_BLUE_MILD = ColorGradientSequential(COLORS.RED_MILD, COLORS.BLUE_MILD)
    SEQUENTIAL_RED_BLUE = ColorGradientSequential(COLORS.RED, COLORS.BLUE)

    # Divergeants
    DIVERGEANT_BLUE_ORANGE_1 = ColorGradientDivergeant(COLORS.BLUE_NICE_1, COLORS.WHITE, COLORS.ORANGE_NICE_1)
    DIVERGEANT_BLUE_ORANGE_2 = ColorGradientDivergeant(COLORS.BLUE_NICE_2, COLORS.WHITE, COLORS.ORANGE_NICE_2)
    DIVERGEANT_RED_BLUE_MILD = ColorGradientDivergeant(COLORS.RED_MILD, COLORS.WHITE, COLORS.BLUE_MILD)
    DIVERGEANT_RED_BLUE = ColorGradientDivergeant(COLORS.RED, COLORS.WHITE, COLORS.BLUE)




# Usage Examples ---------------------------------------------------------------

if __name__ == "__main__":

    # Create and use color
    color = Color(0.8, 0.1, 0.3)
    color_values = color.rgb          # RGB values in float [0.0, 1.0]
    color_values_255 = color.rgb_255  # RGB values in int [0, 255]
    color_hsl = color.hsl()           # HSL vales in float
    print(color, color_values, color_values_255, color_hsl)
    color.show()

    # Saved colors, color modiications and interpolations
    blue = COLORS.BLUE_NICE_1
    blue_dark = blue.updated_lightness(-0.17)             # Update one of the color property   
    blue_dark_desaturated = blue.updated_saturation(-0.9) # Update one of the color property
    blue_palette = [blue, blue_dark, blue_dark_desaturated]
    Color.show_list(blue_palette)                         # Show a list of colors
    
    blue_palette_20 = Color.get_interpolated_range(blue, blue_dark, 10) + Color.get_interpolated_range(blue_dark, blue_dark_desaturated, 10) # Interpolate colors with N steps
    Color.show_list(blue_palette_20)

    # Saved palettes
    my_palette = PALETTES.COLORS3_1 # Palettes is just a list of colors
    Color.show_list(my_palette)
    my_palette_desaturated = [c.updated_lightness(+0.1) for c in my_palette] # Update a color palette
    Color.show_list(my_palette_desaturated)

    # Color gradients: Mapping between a range of floats and a range of colors
    # * Sequential: [0, 1] -> [c1, c2]
    # * Divergeant: [-1, 0, 1] -> [c1, c2, c3]
    cg1 = ColorGradientSequential(COLORS.BLUE_NICE_1, COLORS.ORANGE_NICE_1)
    print(cg1)
    cg1_color_mapping = [cg1.get_color(x) for x in [0.013, 0.21, 0.413, 0.601, 0.899, 1.0]] # Access color gradient value with .get_color()
    print([c.rgb_255 for c in cg1_color_mapping])
    cg1.show()
    
    cg2 = COLOR_GRADIENTS.DIVERGEANT_BLUE_ORANGE_2
    cg2.show()
    
    cg3 = ColorGradientDivergeant(COLORS.RED_MILD, COLORS.WHITE, COLORS.BLUE_MILD)
    cg3.show()