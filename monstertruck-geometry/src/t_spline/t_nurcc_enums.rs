use super::*;
use std::fmt;

impl TnurccConnection {
    /// Converts a `usize` to a `TnurccConnection`. Modulo's `i` by `4`, resulting in a lossy
    /// conversion. Thus, the `From` trait has not been used.
    ///
    /// # Panics
    /// Panics if `i % 4 >= 4` (impossible?).
    pub fn from_usize(i: usize) -> TnurccConnection {
        let remainder = i % 4;
        if remainder == 0 {
            TnurccConnection::LeftCw
        } else if remainder == 1 {
            TnurccConnection::LeftAcw
        } else if remainder == 2 {
            TnurccConnection::RightCw
        } else {
            TnurccConnection::RightAcw
        }
    }
}

impl fmt::Display for TnurccConnection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TnurccConnection::LeftAcw => write!(f, "Anti-clockwise Left"),
            TnurccConnection::RightAcw => write!(f, "Anti-clockwise Right"),
            TnurccConnection::LeftCw => write!(f, "Clockwise Left"),
            TnurccConnection::RightCw => write!(f, "Clockwise Right"),
        }
    }
}
