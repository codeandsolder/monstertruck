use super::TnurccConnection;
use std::fmt;

impl TnurccConnection {
    /// Converts a `usize` to a `TnurccConnection`. Modulo's `i` by `4`, resulting in a lossy
    /// conversion. Thus, the `From` trait has not been used.
    ///
    /// # Panics
    /// Panics if `i % 4 >= 4` (impossible?).
    pub const fn from_usize(i: usize) -> Self {
        let remainder = i % 4;
        if remainder == 0 {
            Self::LeftCw
        } else if remainder == 1 {
            Self::LeftAcw
        } else if remainder == 2 {
            Self::RightCw
        } else {
            Self::RightAcw
        }
    }
}

impl fmt::Display for TnurccConnection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LeftAcw => write!(f, "Anti-clockwise Left"),
            Self::RightAcw => write!(f, "Anti-clockwise Right"),
            Self::LeftCw => write!(f, "Clockwise Left"),
            Self::RightCw => write!(f, "Clockwise Right"),
        }
    }
}
