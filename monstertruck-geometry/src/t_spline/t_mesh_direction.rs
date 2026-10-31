use super::TmeshDirection;
use std::fmt;

impl TmeshDirection {
    /// Returns an iterable object which rotates clockwise starting from UP.
    pub fn iter() -> impl Iterator<Item = Self> {
        [Self::Up, Self::Right, Self::Down, Self::Left]
            .iter()
            .copied()
    }

    /// Rotates `self` a quarter of a rotation clockwise.
    #[must_use]
    pub const fn clockwise(self) -> Self {
        match self {
            Self::Up => Self::Right,
            Self::Right => Self::Down,
            Self::Down => Self::Left,
            Self::Left => Self::Up,
        }
    }

    /// Rotates `self` a quarter of a turn anti-clockwise.
    #[must_use]
    pub const fn anti_clockwise(self) -> Self {
        match self {
            Self::Up => Self::Left,
            Self::Left => Self::Down,
            Self::Down => Self::Right,
            Self::Right => Self::Up,
        }
    }

    /// Flips `self`, rotating it half a turn.
    #[must_use]
    pub const fn flip(self) -> Self {
        match self {
            Self::Up => Self::Down,
            Self::Down => Self::Up,
            Self::Left => Self::Right,
            Self::Right => Self::Left,
        }
    }

    /// Returns true if the direction is UP or DOWN
    #[must_use]
    pub const fn vertical(self) -> bool {
        match self {
            Self::Up | Self::Down => true,
            Self::Left | Self::Right => false,
        }
    }

    /// Returns true if the direction is RIGHT or LEFT
    #[must_use]
    pub const fn horizontal(self) -> bool {
        match self {
            Self::Up | Self::Down => false,
            Self::Left | Self::Right => true,
        }
    }

    /// Returns true if a knot interval in the direction is a positive delta.
    #[must_use]
    pub const fn knot_delta_positive(self) -> bool {
        match self {
            Self::Down | Self::Left => false,
            Self::Up | Self::Right => true,
        }
    }

    /// Returns true if a knot interval in the direction is a negative delta.
    #[must_use]
    pub const fn knot_delta_negative(self) -> bool {
        match self {
            Self::Down | Self::Left => true,
            Self::Up | Self::Right => false,
        }
    }

    /// Adds `delta` to an existing set of knot `coords`, taking into account the direction `self`
    /// such that `delta` is correctly added or subtracted from one of `coords`'s members in order
    /// to move `coords` in the direction `self` by a distance `delta`.
    #[must_use]
    pub fn mutate_knot_coordinates(self, coords: (f64, f64), delta: f64) -> (f64, f64) {
        match self {
            Self::Up => (coords.0, coords.1 + delta),
            Self::Down => (coords.0, coords.1 - delta),
            Self::Left => (coords.0 - delta, coords.1),
            Self::Right => (coords.0 + delta, coords.1),
        }
    }
}

impl fmt::Display for TmeshDirection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let dir_string = match self {
            Self::Up => "up",
            Self::Down => "down",
            Self::Left => "left",
            Self::Right => "right",
        };

        write!(f, "{dir_string}")
    }
}
